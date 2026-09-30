#!/usr/bin/env python
"""Sequence-aware route: nucleotide uniqueness across the 7 genomes.

1. Count canonical 18-mers over all replicons of all 7 genomes (plasmids included,
   since they are potential off-targets). A position is "unique" only if every
   18-mer covering it occurs exactly once in the whole 7-genome set, i.e. it is
   neither repeated in its own genome nor present in any other genome.
2. For each eligible CDS from step 01, find the longest unique run inside it.
3. blastn each candidate gene against its own genome: must have a single hit
   (no secondary HSP >= 80 % id over >= 50 bp) -> nucleotide-level single copy.
Also reports unique intergenic windows (class C) as a fallback.
"""
import os, subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "..", "work")
K = 18
MIN_RUN = 150
THREADS = os.environ.get("THREADS", "1")
COMP = str.maketrans("ACGT", "TGCA")


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True)


allg = os.path.join(WORK, "all_genomes.fna")
jf = os.path.join(WORK, f"k{K}.jf")
rep = os.path.join(WORK, f"k{K}_repeated.txt")
if not os.path.exists(rep):
    sh(f"jellyfish count -C -m {K} -s 100M -t {THREADS} -o {jf} {allg}")
    sh(f"jellyfish dump -c -L 2 {jf} | cut -d' ' -f1 > {rep}")
repeated = set(line.strip() for line in open(rep))
print(f"repeated canonical {K}-mers (count>=2 across 7 genomes): {len(repeated):,}")

cds = pd.read_csv(os.path.join(WORK, "cds_classified.tsv"), sep="\t")
txt = pd.read_csv(os.path.join(WORK, "cds_text.tsv"), sep="\t").fillna("")
cds = cds.merge(txt, on="id", how="left")
cds["text_class"] = cds.text_class.fillna("")
cds["class"] = cds["class"].fillna("")
# candidates from either route: homology (01) or GFF text parsing (01a)
cds["eligible"] = (((cds["class"] != "") | (cds.text_class != "")) & (cds.replicon == "Chromosome")
                   & ~cds.pseudo & ~cds.near_mobile & (cds.end - cds.start + 1 >= 300))
cds["route"] = (cds["class"].where(cds["class"] != "", "-") + "/" + cds.text_class.where(cds.text_class != "", "-"))
strains = sorted(cds.strain.unique())

masks = {}  # (strain, seqid) -> bool array, True = unique position
for s in strains:
    for rec in SeqIO.parse(os.path.join(WORK, f"{s}.chrom.fna"), "fasta"):
        seq = str(rec.seq)
        rc = seq.translate(COMP)[::-1]
        n = len(seq)
        bad_kmer = np.zeros(n - K + 1, dtype=bool)
        for i in range(n - K + 1):
            f = seq[i:i + K]
            r = rc[n - i - K:n - i]
            if (f if f < r else r) in repeated or "N" in f:
                bad_kmer[i] = True
        # a position is masked if any k-mer covering it is repeated
        cov = np.convolve(bad_kmer.astype(np.int32), np.ones(K, dtype=np.int32))[:n]
        seqid = rec.id.split("|", 1)[1]
        masks[(s, seqid)] = cov == 0
        np.save(os.path.join(WORK, f"mask_{s}_{seqid}.npy"), masks[(s, seqid)])
        print(f"{s} {seqid}: {masks[(s, seqid)].mean():.1%} unique positions")


def longest_run(arr):
    """(start, length) of the longest True run in a bool array (0-based)."""
    if not arr.any():
        return 0, 0
    d = np.diff(np.concatenate(([0], arr.astype(np.int8), [0])))
    st, en = np.where(d == 1)[0], np.where(d == -1)[0]
    j = np.argmax(en - st)
    return int(st[j]), int(en[j] - st[j])


rows = []
for i, r in cds[cds.eligible].iterrows():
    m = masks[(r.strain, r.seqid)][r.start - 1:r.end]
    st, ln = longest_run(m)
    rows.append((i, m.mean(), ln, r.start + st, r.start + st + ln - 1))
u = pd.DataFrame(rows, columns=["idx", "unique_frac", "uniq_run_len", "uniq_run_start", "uniq_run_end"]).set_index("idx")
cand = cds.join(u, how="inner")
cand = cand[cand.uniq_run_len >= MIN_RUN].copy()

# --- blastn single-copy check at the nucleotide level ---
genomes = {s: {rec.id.split("|", 1)[1]: rec.seq for rec in SeqIO.parse(os.path.join(WORK, f"{s}.all.fna"), "fasta")}
           for s in strains}
cand["n_self_hits"] = 0
for s in strains:
    db = os.path.join(WORK, f"blastdb_{s}")
    if not os.path.exists(db + ".nsq") and not os.path.exists(db + ".00.nsq"):
        sh(f"makeblastdb -in {WORK}/{s}.all.fna -dbtype nucl -parse_seqids -out {db} > /dev/null")
    sub = cand[cand.strain == s]
    q = os.path.join(WORK, f"cand_genes_{s}.fna")
    with open(q, "w") as fh:
        for i, r in sub.iterrows():
            fh.write(f">{i}\n{genomes[s][r.seqid][r.start - 1:r.end]}\n")
    out = os.path.join(WORK, f"cand_genes_{s}.self.tsv")
    sh(f"blastn -query {q} -db {db} -outfmt '6 qseqid sseqid pident length' -evalue 1e-5 "
       f"-num_threads {THREADS} -max_hsps 50 > {out}")
    b = pd.read_csv(out, sep="\t", header=None, names=["q", "s", "pident", "len"])
    b = b[(b.pident >= 80) & (b.len >= 50)]
    cand.loc[sub.index, "n_self_hits"] = b.groupby("q").size().reindex(sub.index, fill_value=0).values
cand["single_copy_nt"] = cand.n_self_hits == 1
cand.to_csv(os.path.join(WORK, "candidates_genes.tsv"), sep="\t")
print("\ncandidate genes with >=%d bp unique run, single-copy at nt level:" % MIN_RUN)
print(cand[cand.single_copy_nt].groupby(["strain", "route"]).size().unstack(fill_value=0))

# --- class C fallback: unique intergenic windows on chromosomes ---
ig = []
for (s, seqid), m in masks.items():
    genic = np.zeros(len(m), dtype=bool)
    for _, r in cds[(cds.strain == s) & (cds.seqid == seqid)].iterrows():
        genic[r.start - 1:r.end] = True
    free = m & ~genic
    d = np.diff(np.concatenate(([0], free.astype(np.int8), [0])))
    for a, b in zip(np.where(d == 1)[0], np.where(d == -1)[0]):
        if b - a >= MIN_RUN:
            ig.append((s, seqid, a + 1, b, b - a))
pd.DataFrame(ig, columns=["strain", "seqid", "start", "end", "len"]).to_csv(
    os.path.join(WORK, "candidates_intergenic.tsv"), sep="\t", index=False)
print(f"\nunique intergenic windows >= {MIN_RUN} bp: {len(ig)}")
