#!/usr/bin/env python
"""Design primer + hydrolysis-probe sets with primer3 inside unique sequence.

Template = candidate gene (+/- nothing); positions masked in step 02 (non-unique
18-mers) are passed as excluded regions so oligos and amplicon sit in sequence
found only once in the 7-genome set.
Conditions approximate a probe-based dPCR master mix (primers 800 nM, probe 400 nM).
"""
import os
import numpy as np
import pandas as pd
import primer3
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "..", "work")
GENES_PER_STRAIN = 80      # top-ranked genes to design on per strain
SETS_PER_GENE = 3

COND = dict(mv_conc=50.0, dv_conc=3.8, dntp_conc=0.8)
PRIMER_NM, PROBE_NM = 800.0, 400.0

BASE = {
    "PRIMER_TASK": "generic",
    "PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 1, "PRIMER_PICK_INTERNAL_OLIGO": 1,
    "PRIMER_PRODUCT_SIZE_RANGE": [[70, 150]],
    "PRIMER_OPT_SIZE": 20, "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 25,
    "PRIMER_OPT_TM": 60.0, "PRIMER_MIN_TM": 58.0, "PRIMER_MAX_TM": 62.0, "PRIMER_PAIR_MAX_DIFF_TM": 1.5,
    "PRIMER_MAX_POLY_X": 4, "PRIMER_MAX_END_GC": 3, "PRIMER_MAX_END_STABILITY": 9.0,
    "PRIMER_MAX_SELF_ANY_TH": 35.0, "PRIMER_MAX_SELF_END_TH": 30.0, "PRIMER_PAIR_MAX_COMPL_ANY_TH": 35.0,
    "PRIMER_PAIR_MAX_COMPL_END_TH": 30.0, "PRIMER_MAX_HAIRPIN_TH": 35.0,
    "PRIMER_INTERNAL_OPT_SIZE": 24, "PRIMER_INTERNAL_MIN_SIZE": 18, "PRIMER_INTERNAL_MAX_SIZE": 30,
    "PRIMER_INTERNAL_OPT_TM": 68.0, "PRIMER_INTERNAL_MIN_TM": 66.0, "PRIMER_INTERNAL_MAX_TM": 71.0,
    "PRIMER_INTERNAL_MAX_POLY_X": 4, "PRIMER_INTERNAL_MAX_SELF_ANY_TH": 40.0,
    "PRIMER_INTERNAL_MAX_HAIRPIN_TH": 40.0,
    "PRIMER_SALT_MONOVALENT": COND["mv_conc"], "PRIMER_SALT_DIVALENT": COND["dv_conc"],
    "PRIMER_DNTP_CONC": COND["dntp_conc"], "PRIMER_DNA_CONC": PRIMER_NM,
    "PRIMER_INTERNAL_SALT_MONOVALENT": COND["mv_conc"], "PRIMER_INTERNAL_SALT_DIVALENT": COND["dv_conc"],
    "PRIMER_INTERNAL_DNTP_CONC": COND["dntp_conc"], "PRIMER_INTERNAL_DNA_CONC": PROBE_NM,
    "PRIMER_NUM_RETURN": 10,
}


def gc_limits(genome_gc):
    """Widen the GC window toward the genome composition (AT-rich / GC-rich genomes)."""
    lo, hi = 35.0, 65.0
    if genome_gc < 0.45:
        lo = 25.0
    if genome_gc > 0.62:
        hi = 72.0
    return lo, hi


def orient_probe(p):
    """Probe may sit on either strand: avoid 5' G, prefer more C than G."""
    rc = p.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    opts = [(p, "+"), (rc, "-")]
    opts = [o for o in opts if not o[0].startswith("G")] or opts
    opts.sort(key=lambda o: -(o[0].count("C") - o[0].count("G")))
    return opts[0]


def excluded(mask):
    """Masked (non-unique) positions -> primer3 [start, len] list (0-based)."""
    d = np.diff(np.concatenate(([0], (~mask).astype(np.int8), [0])))
    return [[int(a), int(b - a)] for a, b in zip(np.where(d == 1)[0], np.where(d == -1)[0])]


cand = pd.read_csv(os.path.join(WORK, "candidates_genes.tsv"), sep="\t", index_col=0)
cand = cand[cand.single_copy_nt]
chrom = {}
for f in os.listdir(WORK):
    if f.endswith(".chrom.fna"):
        for rec in SeqIO.parse(os.path.join(WORK, f), "fasta"):
            s, sid = rec.id.split("|", 1)
            chrom[(s, sid)] = str(rec.seq)
genome_gc = {s: gc_fraction("".join(v for (ss, _), v in chrom.items() if ss == s)) for s in cand.strain.unique()}

# rank genes: single-copy core genes (homology 1:1 or named 1x in all genomes) first,
# then unique genes supported by both routes, homology only, text only;
# then longer unique run and lower identity to the closest homolog in the other genomes
def class_rank(r):
    c, t = r["class"] if isinstance(r["class"], str) else "", r.text_class if isinstance(r.text_class, str) else ""
    if c == "B_universal_1to1" and t == "text_core_1x":
        return 0
    if c == "B_universal_1to1" or t == "text_core_1x":
        return 1
    if c == "A_unique" and t == "text_unique":
        return 2
    return 3 if c == "A_unique" else 4


cand["class_rank"] = cand.apply(class_rank, axis=1)
cand = cand.sort_values(["strain", "class_rank", "uniq_run_len", "best_other_pident"],
                        ascending=[True, True, False, True])

out = []
for s, g in cand.groupby("strain"):
    lo, hi = gc_limits(genome_gc[s])
    args = dict(BASE, PRIMER_MIN_GC=lo, PRIMER_MAX_GC=hi, PRIMER_INTERNAL_MIN_GC=lo, PRIMER_INTERNAL_MAX_GC=hi)
    n_genes = 0
    for idx, r in g.iterrows():
        if n_genes >= GENES_PER_STRAIN:
            break
        mask = np.load(os.path.join(WORK, f"mask_{s}_{r.seqid}.npy"))[r.start - 1:r.end]
        tmpl = chrom[(s, r.seqid)][r.start - 1:r.end]
        try:
            res = primer3.bindings.design_primers(
                {"SEQUENCE_ID": f"{s}_{r.locus_tag}", "SEQUENCE_TEMPLATE": tmpl,
                 "SEQUENCE_EXCLUDED_REGION": excluded(mask)}, args)
        except OSError:
            continue
        n = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
        kept, used = 0, set()
        for k in range(n):
            L, R, P = (res[f"PRIMER_LEFT_{k}"], res[f"PRIMER_RIGHT_{k}"], res[f"PRIMER_INTERNAL_{k}"])
            # diversify: skip sets overlapping an already kept amplicon start
            if any(abs(L[0] - u) < 40 for u in used):
                continue
            fseq, rseq = res[f"PRIMER_LEFT_{k}_SEQUENCE"], res[f"PRIMER_RIGHT_{k}_SEQUENCE"]
            pseq, pstrand = orient_probe(res[f"PRIMER_INTERNAL_{k}_SEQUENCE"])
            if pseq.startswith("G"):
                continue
            amp_start, amp_end = L[0], R[0]  # 0-based in template, R[0] = 3'-most base
            amp = tmpl[amp_start:amp_end + 1]
            gstart = r.start + amp_start
            out.append(dict(
                set_id=f"{s}_{r.locus_tag}_{k}", strain=s, route=r.route, class_rank=r.class_rank, locus_tag=r.locus_tag,
                gene=r.gene if isinstance(r.gene, str) else "", product=r["product"], seqid=r.seqid,
                gene_start=r.start, gene_end=r.end, gene_strand=r.strand,
                amp_start=gstart, amp_end=gstart + len(amp) - 1, amp_len=len(amp),
                fwd=fseq, rev=rseq, probe=pseq, probe_strand=pstrand,
                fwd_tm=res[f"PRIMER_LEFT_{k}_TM"], rev_tm=res[f"PRIMER_RIGHT_{k}_TM"],
                probe_tm=res[f"PRIMER_INTERNAL_{k}_TM"],
                fwd_gc=res[f"PRIMER_LEFT_{k}_GC_PERCENT"], rev_gc=res[f"PRIMER_RIGHT_{k}_GC_PERCENT"],
                probe_gc=res[f"PRIMER_INTERNAL_{k}_GC_PERCENT"], amp_gc=100 * gc_fraction(amp),
                penalty=res[f"PRIMER_PAIR_{k}_PENALTY"], best_other_pident=r.best_other_pident,
                amplicon=amp))
            used.add(L[0]); kept += 1
            if kept >= SETS_PER_GENE:
                break
        if kept:
            n_genes += 1
    print(f"{s}: GC {genome_gc[s]:.2f}, primer GC {lo}-{hi}%, designed on {n_genes} genes")

df = pd.DataFrame(out)
df.to_csv(os.path.join(WORK, "designs_raw.tsv"), sep="\t", index=False)
print(df.groupby(["strain", "class_rank"]).size().unstack(fill_value=0))
