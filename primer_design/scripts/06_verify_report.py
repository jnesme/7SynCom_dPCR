#!/usr/bin/env python
"""Verification of the final pool + report files.

1. In-silico PCR of the pooled primers (every primer combination, <= 3 mismatches)
   on all 7 genomes: expect exactly 7 products, one per genome, of the designed length.
2. Knock-out test: for each strain, replace its target amplicon by Ns and confirm the
   product disappears while the 6 others remain (no hidden alternative site).
3. amplicons.fasta: amplicon +/- 20 bp flanks, a template for gBlock standards.
4. oligos.tsv (order sheet) and dimer heatmap PNG.
"""
import itertools, os, subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "..", "work")
RES = os.path.join(HERE, "..", "results")
THREADS = os.environ.get("THREADS", "1")
FLANK = 20
MAX_PROD = 3000


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True)


sel = pd.read_csv(os.path.join(RES, "final_multiplex.tsv"), sep="\t")
genomes = {rec.id: str(rec.seq) for rec in SeqIO.parse(os.path.join(WORK, "all_genomes.fna"), "fasta")}

prim = [(f"{r.strain}_F", r.fwd) for r in sel.itertuples()] + [(f"{r.strain}_R", r.rev) for r in sel.itertuples()]
ptsv = os.path.join(WORK, "final_pairs.tsv")
with open(ptsv, "w") as fh:
    for (na, a), (nb, b) in itertools.combinations_with_replacement(prim, 2):
        fh.write(f"{na}+{nb}\t{a}\t{b}\n")


def pcr(fasta):
    bed = fasta + ".bed"
    sh(f"seqkit amplicon -j {THREADS} -m 3 -p {ptsv} --bed {fasta} > {bed}")
    p = pd.read_csv(bed, sep="\t", header=None, names=["chrom", "start", "end", "name", "score", "strand", "seq"])
    p = p[(p.end - p.start) <= MAX_PROD]
    return p.drop_duplicates(["chrom", "start", "end"])


lines = []
p = pcr(os.path.join(WORK, "all_genomes.fna"))
p["strain"] = p.chrom.str.split("|").str[0]
lines.append(f"Pooled in-silico PCR (all {len(prim)} primers, all combinations, <=3 mismatches): {len(p)} products")
for r in p.itertuples():
    lines.append(f"  {r.chrom}:{r.start + 1}-{r.end} ({r.end - r.start} bp) via {r.name}")
exp = {(f"{r.strain}|{r.seqid}", r.amp_start - 1, r.amp_end) for r in sel.itertuples()}
ok1 = set(zip(p.chrom, p.start, p.end)) == exp
lines.append(f"  -> exactly the 7 designed products: {'PASS' if ok1 else 'FAIL'}")

ok2 = True
for r in sel.itertuples():
    ko = os.path.join(WORK, f"ko_{r.strain}.fna")
    with open(ko, "w") as fh:
        for cid, seq in genomes.items():
            if cid == f"{r.strain}|{r.seqid}":
                seq = seq[:r.amp_start - 1] + "N" * r.amp_len + seq[r.amp_end:]
            fh.write(f">{cid}\n{seq}\n")
    q = pcr(ko)
    q["strain"] = q.chrom.str.split("|").str[0]
    good = (len(q) == 6) and (r.strain not in set(q.strain))
    ok2 &= good
    lines.append(f"Knock-out {r.strain}: {len(q)} products, {r.strain} product absent: {'PASS' if good else 'FAIL'}")
    os.remove(ko); os.remove(ko + ".bed")

with open(os.path.join(RES, "amplicons.fasta"), "w") as fh:
    for r in sel.itertuples():
        g = genomes[f"{r.strain}|{r.seqid}"]
        a, b = max(0, r.amp_start - 1 - FLANK), min(len(g), r.amp_end + FLANK)
        fh.write(f">{r.strain}_{r.locus_tag} {r.species_NCBI} {r.seqid}:{a + 1}-{b} "
                 f"amplicon {r.amp_len} bp + {FLANK} bp flanks\n{g[a:b]}\n")

ol = []
for r in sel.itertuples():
    for role, seq, tm in (("F", r.fwd, r.fwd_tm), ("R", r.rev, r.rev_tm), ("P", r.probe, r.probe_tm)):
        ol.append(dict(name=f"{r.strain}_{r.locus_tag}_{role}", strain=r.strain, species=r.species_NCBI,
                       role={"F": "forward primer", "R": "reverse primer", "P": "hydrolysis probe (5' dye / 3' quencher)"}[role],
                       sequence_5to3=seq, length=len(seq), tm=round(tm, 1),
                       gc=round(100 * (seq.count("G") + seq.count("C")) / len(seq), 1), well_4plus3=r.well_4plus3))
pd.DataFrame(ol).to_csv(os.path.join(RES, "oligos.tsv"), sep="\t", index=False)

m = pd.read_csv(os.path.join(RES, "dimer_matrix.tsv"), sep="\t", index_col=0)
fig, ax = plt.subplots(figsize=(9, 8))
im = ax.imshow(m.values, cmap="magma", vmin=min(-12, m.values.min()), vmax=0)
ax.set_xticks(range(len(m))); ax.set_xticklabels(m.columns, rotation=90, fontsize=7)
ax.set_yticks(range(len(m))); ax.set_yticklabels(m.index, fontsize=7)
fig.colorbar(im, ax=ax, label="heterodimer dG (kcal/mol, 37 C)")
ax.set_title("Final pool: pairwise oligo interactions")
fig.tight_layout(); fig.savefig(os.path.join(RES, "dimer_heatmap.png"), dpi=150)

off_diag = m.values[~np.eye(len(m), dtype=bool)]
lines.append(f"Worst oligo-oligo heterodimer dG in pool (excluding self): {off_diag.min():.2f} kcal/mol")
open(os.path.join(RES, "verification.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
if not (ok1 and ok2):
    raise SystemExit("verification FAILED")
