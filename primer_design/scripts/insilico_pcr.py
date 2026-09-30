"""Exhaustive in-silico PCR with bowtie1.

`seqkit amplicon` is not usable here: with mismatches allowed it reports a single
product per primer pair and strand (leftmost forward hit), which can hide the true
on-target product behind a spurious Mb-long one. Instead, every binding site with
<= MM mismatches is enumerated with bowtie1 (-v MM -a, both strands, exhaustive),
and products are formed by pairing any '+' site with any downstream '-' site.
Any oligo can prime from either orientation, so F+F and R+R products are included.
"""
import os, subprocess
import numpy as np
import pandas as pd

MM = 3


def build_index(fasta, prefix, threads="1"):
    if not os.path.exists(prefix + ".1.ebwt"):
        subprocess.run(f"bowtie-build --threads {threads} -q {fasta} {prefix}", shell=True, check=True)


def find_sites(oligos, index, work, tag, mm=MM, threads="1"):
    """oligos: dict name -> seq. Returns one row per binding site:
    name, chrom, strand, start (0-based leftmost on the reference), end, mm, mm3
    (mm3 = mismatches among the 5 3'-terminal bases of the oligo)."""
    fa = os.path.join(work, f"{tag}.oligos.fa")
    out = os.path.join(work, f"{tag}.bowtie.tsv")
    with open(fa, "w") as fh:
        for n, s in oligos.items():
            fh.write(f">{n}\n{s}\n")
    subprocess.run(f"bowtie -f -v {mm} -a -p {threads} --quiet -x {index} {fa} > {out}", shell=True, check=True)
    cols = ["name", "strand", "chrom", "start", "seq", "qual", "n_other", "mism"]
    if os.path.getsize(out) == 0:
        return pd.DataFrame(columns=["name", "chrom", "strand", "start", "end", "mm", "mm3"])
    b = pd.read_csv(out, sep="\t", header=None, names=cols, keep_default_na=False, dtype={"mism": str})
    b["olen"] = b.name.map({n: len(s) for n, s in oligos.items()})
    b["end"] = b.start + b.olen
    # bowtie1 mismatch offsets are relative to the 5' end of the read as given
    offs = b.mism.apply(lambda m: [int(x.split(":")[0]) for x in m.split(",")] if m else [])
    b["mm"] = offs.str.len()
    b["mm3"] = [sum(o >= L - 5 for o in os_) for os_, L in zip(offs, b.olen)]
    return b[["name", "chrom", "strand", "start", "end", "mm", "mm3"]]


def products(sites, max_len=3000, min_len=1):
    """Pair every '+' site with every '-' site downstream on the same sequence."""
    rows = []
    for chrom, g in sites.groupby("chrom"):
        plus = g[g.strand == "+"].sort_values("start")
        minus = g[g.strand == "-"].sort_values("end")
        if plus.empty or minus.empty:
            continue
        mend = minus.end.values
        for p in plus.itertuples():
            lo = np.searchsorted(mend, p.start + min_len, side="left")
            hi = np.searchsorted(mend, p.start + max_len, side="right")
            for m in minus.iloc[lo:hi].itertuples():
                if m.start >= p.start:
                    rows.append((chrom, p.start, m.end, m.end - p.start, p.name, m.name,
                                 p.mm, m.mm, p.mm3, m.mm3))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "len", "left_oligo", "right_oligo",
                                       "left_mm", "right_mm", "left_mm3", "right_mm3"])
