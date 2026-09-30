#!/usr/bin/env python
"""Verification of the final pool + report files.

1. Exhaustive in-silico PCR of the 14 pooled primers (<= 3 mismatches, every primer
   combination incl. F+F/R+R) on all 7 genomes: expect exactly 7 products, one per genome.
2. Knock-out test: for each strain, drop every primer site overlapping its target
   amplicon (equivalent to replacing the amplicon by Ns) and confirm its product
   disappears while the 6 others remain (no hidden alternative site).
3. Restriction enzymes commonly used to fragment gDNA for dPCR that cut none of the
   7 amplicons (and how often they cut each genome).
4. amplicons.fasta (amplicon +/- 20 bp, template for gBlock standards), oligos.tsv
   (order sheet), dimer heatmap.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq
from insilico_pcr import build_index, find_sites, products

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "..", "work")
RES = os.path.join(HERE, "..", "results")
THREADS = os.environ.get("THREADS", "1")
FLANK = 20
MAX_PROD = 3000
ENZYMES = ["HaeIII", "MseI", "CviQI", "AluI", "EcoRI", "HindIII", "XbaI", "PvuII", "BamHI", "MspI", "Sau3AI", "BsaI"]
allg = os.path.join(WORK, "all_genomes.fna")
idx = os.path.join(WORK, "bt_all")

sel = pd.read_csv(os.path.join(RES, "final_multiplex.tsv"), sep="\t")
genomes = {rec.id: rec.seq for rec in SeqIO.parse(allg, "fasta")}
build_index(allg, idx, THREADS)
sites = find_sites({f"{r.strain}_{role}": getattr(r, role) for r in sel.itertuples() for role in ("fwd", "rev")},
                   idx, WORK, "final", threads=THREADS)
exp = {(f"{r.strain}|{r.seqid}", r.amp_start - 1, r.amp_end) for r in sel.itertuples()}

lines = []
p = products(sites, max_len=MAX_PROD).drop_duplicates(["chrom", "start", "end"])
lines.append(f"Pooled in-silico PCR (14 primers, all combinations, <=3 mismatches, <= {MAX_PROD} bp): {len(p)} products")
for r in p.itertuples():
    lines.append(f"  {r.chrom}:{r.start + 1}-{r.end} ({r.len} bp) {r.left_oligo}+{r.right_oligo} "
                 f"mm {r.left_mm}/{r.right_mm}")
ok1 = set(zip(p.chrom, p.start, p.end)) == exp
lines.append(f"  -> exactly the 7 designed products: {'PASS' if ok1 else 'FAIL'}")

ok2 = True
for r in sel.itertuples():
    hit = (sites.chrom == f"{r.strain}|{r.seqid}") & (sites.end > r.amp_start - 1) & (sites.start < r.amp_end)
    q = products(sites[~hit], max_len=MAX_PROD).drop_duplicates(["chrom", "start", "end"])
    strains_left = set(q.chrom.str.split("|").str[0])
    good = len(q) == 6 and r.strain not in strains_left
    ok2 &= good
    lines.append(f"Knock-out {r.strain}: {len(q)} products remain, {r.strain} product absent: {'PASS' if good else 'FAIL'}")

rb = RestrictionBatch(ENZYMES)
amps = {r.strain: Seq(r.amplicon) for r in sel.itertuples()}
lines.append("Restriction enzymes (gDNA fragmentation) vs amplicons / genome cut frequency:")
for enz in sorted(rb, key=str):
    cut_amps = [s for s, a in amps.items() if enz.search(a)]
    med = np.median([len(enz.search(g, linear=False)) for g in genomes.values()])
    lines.append(f"  {str(enz):8s} site {enz.site:10s} cuts amplicons: {','.join(cut_amps) or 'none':20s} "
                 f"median cuts/replicon {med:,.0f}{'   <- usable' if not cut_amps else ''}")

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
