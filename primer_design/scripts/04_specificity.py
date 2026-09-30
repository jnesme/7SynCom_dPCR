#!/usr/bin/env python
"""In-silico specificity of every designed set against all 7 genomes (plasmids included).

a) in-silico PCR (seqkit amplicon, <= MM mismatches per primer, both strands):
   exactly one product <= MAX_PROD bp, at the designed locus.
b) blastn-short of each oligo: an off-target binding site is "dangerous" when it has
   < 4 mismatches in total AND < 2 mismatches in the 5 3'-terminal bases (primers).
   Sets with any dangerous primer site are rejected; probe off-target identity is recorded
   (probe cross-talk with the other final amplicons is checked in step 05).
"""
import os, subprocess
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "..", "work")
RES = os.path.join(HERE, "..", "results")
os.makedirs(RES, exist_ok=True)
THREADS = os.environ.get("THREADS", "1")
MM = 3
MAX_PROD = 3000
allg = os.path.join(WORK, "all_genomes.fna")


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True)


d = pd.read_csv(os.path.join(WORK, "designs_raw.tsv"), sep="\t")
d = d.drop_duplicates(["fwd", "rev"])

# ---- a) in-silico PCR -------------------------------------------------------
ptsv = os.path.join(WORK, "designs_primers.tsv")
d[["set_id", "fwd", "rev"]].to_csv(ptsv, sep="\t", header=False, index=False)
bed = os.path.join(WORK, "designs_insilico_pcr.bed")
sh(f"seqkit amplicon -j {THREADS} -m {MM} -p {ptsv} --bed {allg} > {bed}")
p = pd.read_csv(bed, sep="\t", header=None,
                names=["chrom", "start", "end", "set_id", "score", "strand", "seq"])
p["len"] = p.end - p.start
p = p[p.len <= MAX_PROD]
p["strain"] = p.chrom.str.split("|").str[0]
p["seqid"] = p.chrom.str.split("|").str[1]

pcr = []
for _, r in d.iterrows():
    h = p[p.set_id == r.set_id]
    on = h[(h.strain == r.strain) & (h.seqid == r.seqid) & (h.start + 1 == r.amp_start) & (h.end == r.amp_end)]
    pcr.append((r.set_id, len(h), len(on), len(h) - len(on),
                ";".join(f"{x.chrom}:{x.start + 1}-{x.end}" for x in h.itertuples() if x.Index not in on.index)))
pcr = pd.DataFrame(pcr, columns=["set_id", "n_products", "n_on_target", "n_off_products", "off_products"])

# ---- b) blastn-short per oligo ---------------------------------------------
db = os.path.join(WORK, "blastdb_all")
if not os.path.exists(db + ".nsq"):
    sh(f"makeblastdb -in {allg} -dbtype nucl -out {db} > /dev/null")
olig = os.path.join(WORK, "designs_oligos.fna")
with open(olig, "w") as fh:
    for r in d.itertuples():
        for role in ("fwd", "rev", "probe"):
            fh.write(f">{r.set_id}__{role}\n{getattr(r, role)}\n")
bo = os.path.join(WORK, "designs_oligos.blast.tsv")
sh(f"blastn -task blastn-short -query {olig} -db {db} -word_size 7 -evalue 1000 -dust no "
   f"-max_target_seqs 50 -max_hsps 200 -num_threads {THREADS} "
   f"-outfmt '6 qseqid sseqid qstart qend sstart send qlen qseq sseq' > {bo}")
b = pd.read_csv(bo, sep="\t", header=None,
                names=["q", "s", "qstart", "qend", "sstart", "send", "qlen", "qseq", "sseq"])
b["set_id"] = b.q.str.rsplit("__", n=1).str[0]
b["role"] = b.q.str.rsplit("__", n=1).str[1]


def mism(r):
    """mismatches over the whole oligo (unaligned ends count) and in the 3'-terminal 5 nt."""
    mm_pos = set(range(1, r.qstart)) | set(range(r.qend + 1, r.qlen + 1))
    qi = r.qstart
    for a, c in zip(r.qseq, r.sseq):
        if a == "-":
            mm_pos.add(qi); continue  # insertion in subject: attribute to current position
        if c == "-" or a != c:
            mm_pos.add(qi)
        qi += 1
    return len(mm_pos), sum(1 for x in mm_pos if x > r.qlen - 5)


b[["mm", "mm3"]] = b.apply(lambda r: pd.Series(mism(r)), axis=1)
b["lo"] = b[["sstart", "send"]].min(axis=1)
b["hi"] = b[["sstart", "send"]].max(axis=1)
dd = d.set_index("set_id")
b = b.join(dd[["strain", "seqid", "amp_start", "amp_end"]], on="set_id")
intended = (b.s == b.strain + "|" + b.seqid) & (b.lo >= b.amp_start) & (b.hi <= b.amp_end)
off = b[~intended].copy()
off["dangerous"] = (off.mm < 4) & (off.mm3 < 2)
agg = off[off.role != "probe"].groupby("set_id").dangerous.sum().rename("n_dangerous_primer_sites")
pmin = off[off.role == "probe"].groupby("set_id").mm.min().rename("probe_min_offtarget_mm")
pmm = off[off.role != "probe"].groupby("set_id").mm.min().rename("primer_min_offtarget_mm")

rep = d.merge(pcr, on="set_id").join(agg, on="set_id").join(pmm, on="set_id").join(pmin, on="set_id")
rep["n_dangerous_primer_sites"] = rep.n_dangerous_primer_sites.fillna(0).astype(int)
rep[["primer_min_offtarget_mm", "probe_min_offtarget_mm"]] = rep[["primer_min_offtarget_mm", "probe_min_offtarget_mm"]].fillna(99)
rep["pass"] = (rep.n_products == 1) & (rep.n_on_target == 1) & (rep.n_dangerous_primer_sites == 0)
rep.to_csv(os.path.join(RES, "specificity_report.tsv"), sep="\t", index=False)
print(rep.groupby("strain")["pass"].agg(["size", "sum"]).rename(columns={"size": "designed", "sum": "pass"}))
