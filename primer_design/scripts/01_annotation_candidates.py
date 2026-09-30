#!/usr/bin/env python
"""Annotation route: classify CDSs by homology across the 7 genomes.

Class A  = strain-unique: no homolog in any other genome, no paralog in own genome.
Class B  = universal single-copy: no paralog in own genome and exactly one homolog
           in each of the 6 other genomes (1:1 orthologs, e.g. housekeeping genes).
Mobile/unstable loci (and genes within FLANK bp of them) and plasmid genes are flagged.
"""
import os, re, subprocess, sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "..", "work")
FLANK = 5000
EVALUE = 1e-5
MIN_COV = 0.3  # alignment coverage of the shorter protein

m8 = os.path.join(WORK, "allvsall.m8")
if not os.path.exists(m8):
    subprocess.run(["mmseqs", "easy-search", os.path.join(WORK, "all_proteins.faa"),
                    os.path.join(WORK, "all_proteins.faa"), m8, os.path.join(WORK, "mmseqs_tmp"),
                    "-s", "7.5", "-e", str(EVALUE), "--max-seqs", "1000",
                    "--format-output", "query,target,pident,alnlen,evalue,bits,qlen,tlen",
                    "--threads", os.environ.get("THREADS", "1")], check=True)

h = pd.read_csv(m8, sep="\t", header=None,
                names=["q", "t", "pident", "alnlen", "evalue", "bits", "qlen", "tlen"])
h = h[(h.q != h.t) & (h.evalue <= EVALUE)]
h = h[h.alnlen / h[["qlen", "tlen"]].min(axis=1) >= MIN_COV]
h["qs"] = h.q.str.split("|").str[0]
h["ts"] = h.t.str.split("|").str[0]

cds = pd.read_csv(os.path.join(WORK, "cds_table.tsv"), sep="\t")
cds["id"] = cds.strain + "|" + cds.locus_tag
strains = sorted(cds.strain.unique())

counts = h.groupby(["q", "ts"]).t.nunique().unstack(fill_value=0)
counts = counts.reindex(index=cds.id, columns=strains, fill_value=0)
own = [counts.at[i, s] for i, s in zip(cds.id, cds.strain)]
cds["n_paralogs"] = own
others = counts.values.copy()
for j, s in enumerate(strains):
    others[(cds.strain == s).values, j] = -1  # blank own column
cds["n_genomes_with_homolog"] = (others > 0).sum(axis=1)
cds["max_homologs_other"] = others.max(axis=1)
cds["best_other_pident"] = cds.id.map(h[h.qs != h.ts].groupby("q").pident.max()).fillna(0)

cds["class"] = ""
cds.loc[(cds.n_paralogs == 0) & (cds.n_genomes_with_homolog == 0), "class"] = "A_unique"
cds.loc[(cds.n_paralogs == 0) & (cds.n_genomes_with_homolog == len(strains) - 1)
        & (cds.max_homologs_other == 1), "class"] = "B_universal_1to1"

mobile_re = re.compile(r"transpos|integrase|insertion element|insertion sequence|phage|prophage|"
                       r"site-specific recombinase|resolvase|invertase|conjugal|conjugative|relaxase|"
                       r"mobiliz|toxin|antitoxin|reverse transcriptase|plasmid", re.I)
is_re = re.compile(r"\bIS\d|\bIS[A-Z]")  # IS families, case-sensitive (avoid 'isomerase')
cds["mobile"] = (cds["product"].fillna("").str.contains(mobile_re) | cds["product"].fillna("").str.contains(is_re)) & ~cds["product"].fillna("").str.contains("Holliday")
cds["near_mobile"] = False
for (s, sid), g in cds.groupby(["strain", "seqid"]):
    mob = g[g.mobile]
    for i, r in g.iterrows():
        if ((mob.start - FLANK <= r.end) & (mob.end + FLANK >= r.start)).any():
            cds.at[i, "near_mobile"] = True

cds["eligible"] = ((cds["class"] != "") & (cds.replicon == "Chromosome") & ~cds.pseudo
                   & ~cds.near_mobile & (cds.end - cds.start + 1 >= 300))
cds.to_csv(os.path.join(WORK, "cds_classified.tsv"), sep="\t", index=False)
print(cds.groupby(["strain", "class"]).size().unstack(fill_value=0))
print("\neligible (chromosomal, not mobile-adjacent, >=300 bp):")
print(cds[cds.eligible].groupby(["strain", "class"]).size().unstack(fill_value=0))
