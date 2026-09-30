#!/usr/bin/env python
"""In-silico specificity of every designed set against all 7 genomes (plasmids included).

All binding sites with <= 3 mismatches of every primer and probe are enumerated
exhaustively with bowtie1 (see insilico_pcr.py; seqkit amplicon is not exhaustive).
a) in-silico PCR: pairing the set's own F/R sites (incl. F+F, R+R) must give exactly
   one product <= MAX_PROD bp, at the designed locus.
b) a lone off-target primer site cannot amplify without a partner site (that case is
   caught by a), where 3 mismatches per primer are tolerated), but a near-perfect site
   competes for primer and can seed mispriming. A site is "dangerous" when it has
   <= 1 mismatch, or 2 mismatches none of which is in the 5 3'-terminal bases; sets with
   any dangerous primer site are rejected. Sites with 3 mismatches are the random
   background for 18-25-mers in 36 Mb (~13 per primer) and are only counted
   (n_3mm_sites, used for ranking). Probe off-target mismatches are recorded; probe
   cross-talk with the other final amplicons is checked in step 05.
"""
import os
import pandas as pd
from insilico_pcr import build_index, find_sites, products

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "..", "work")
RES = os.path.join(HERE, "..", "results")
os.makedirs(RES, exist_ok=True)
THREADS = os.environ.get("THREADS", "1")
MAX_PROD = 3000
allg = os.path.join(WORK, "all_genomes.fna")
idx = os.path.join(WORK, "bt_all")

d = pd.read_csv(os.path.join(WORK, "designs_raw.tsv"), sep="\t").drop_duplicates(["fwd", "rev"])
build_index(allg, idx, THREADS)
oligos = {f"{r.set_id}__{role}": getattr(r, role) for r in d.itertuples() for role in ("fwd", "rev", "probe")}
sites = find_sites(oligos, idx, WORK, "designs", threads=THREADS)
sites["set_id"] = sites.name.str.rsplit("__", n=1).str[0]
sites["role"] = sites.name.str.rsplit("__", n=1).str[1]
sites = sites.join(d.set_index("set_id")[["strain", "seqid", "amp_start", "amp_end"]], on="set_id")
# intended site = inside the designed amplicon (amp_start/amp_end are 1-based inclusive)
sites["intended"] = ((sites.chrom == sites.strain + "|" + sites.seqid)
                     & (sites.start >= sites.amp_start - 1) & (sites.end <= sites.amp_end))

rows = []
for sid, g in sites[sites.role != "probe"].groupby("set_id"):
    r = d.set_index("set_id").loc[sid]
    p = products(g, max_len=MAX_PROD)
    on = p[(p.chrom == f"{r.strain}|{r.seqid}") & (p.start == r.amp_start - 1) & (p.end == r.amp_end)]
    off = p.drop(on.index)
    rows.append((sid, len(p), len(on), len(off),
                 ";".join(f"{x.chrom}:{x.start + 1}-{x.end}({x.left_oligo.split('__')[1]}+{x.right_oligo.split('__')[1]})"
                          for x in off.itertuples())))
pcr = pd.DataFrame(rows, columns=["set_id", "n_products", "n_on_target", "n_off_products", "off_products"])

off = sites[~sites.intended].copy()
off["dangerous"] = (off.mm <= 1) | ((off.mm == 2) & (off.mm3 == 0))
agg = off[off.role != "probe"].groupby("set_id").dangerous.sum().rename("n_dangerous_primer_sites")
pmm = off[off.role != "probe"].groupby("set_id").mm.min().rename("primer_min_offtarget_mm")
n3 = off[(off.role != "probe") & (off.mm == 3)].groupby("set_id").size().rename("n_3mm_sites")
prb = off[off.role == "probe"].groupby("set_id").mm.min().rename("probe_min_offtarget_mm")

rep = d.merge(pcr, on="set_id", how="left").join(agg, on="set_id").join(pmm, on="set_id").join(prb, on="set_id").join(n3, on="set_id")
cnt = ["n_products", "n_on_target", "n_off_products", "n_dangerous_primer_sites", "n_3mm_sites"]
rep[cnt] = rep[cnt].fillna(0).astype(int)
# no site with <= 3 mismatches found -> >= 4
rep[["primer_min_offtarget_mm", "probe_min_offtarget_mm"]] = rep[["primer_min_offtarget_mm", "probe_min_offtarget_mm"]].fillna(4).astype(int)
rep["pass"] = (rep.n_products == 1) & (rep.n_on_target == 1) & (rep.n_dangerous_primer_sites == 0)
rep.to_csv(os.path.join(RES, "specificity_report.tsv"), sep="\t", index=False)
print(rep.groupby("strain").agg(designed=("pass", "size"), on_target_ok=("n_on_target", lambda x: (x == 1).sum()),
                                no_off_products=("n_off_products", lambda x: (x == 0).sum()),
                                no_dangerous_sites=("n_dangerous_primer_sites", lambda x: (x == 0).sum()),
                                passed=("pass", "sum")))
