#!/usr/bin/env python
"""Prepare per-strain inputs: chromosome/plasmid FASTA, CDS table, prefixed proteins."""
import csv, glob, gzip, os, sys
from Bio import SeqIO

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
GEN = os.path.join(ROOT, "genomes")
WORK = os.path.join(ROOT, "primer_design", "work")
os.makedirs(WORK, exist_ok=True)

strains = []
with open(os.path.join(GEN, "assemblies.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        strains.append(r)

cds_rows = []
rep_rows = []
all_prot = open(os.path.join(WORK, "all_proteins.faa"), "w")
all_genomes = open(os.path.join(WORK, "all_genomes.fna"), "w")
for s in strains:
    tag = s["strain"]
    d = glob.glob(os.path.join(GEN, s["GenBank_acc"], "refseq_*"))[0]
    pre = glob.glob(os.path.join(d, "*_genomic.fna.gz"))
    fna = [p for p in pre if "cds_from" not in p][0]
    gff = glob.glob(os.path.join(d, "*_genomic.gff.gz"))[0]
    faa = glob.glob(os.path.join(d, "*_protein.faa.gz"))[0]
    rep = glob.glob(os.path.join(d, "*_assembly_report.txt"))[0]

    role = {}
    for line in open(rep):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        role[f[6]] = f[3]  # RefSeq accn -> Chromosome/Plasmid
    chrom = open(os.path.join(WORK, f"{tag}.chrom.fna"), "w")
    full = open(os.path.join(WORK, f"{tag}.all.fna"), "w")
    for rec in SeqIO.parse(gzip.open(fna, "rt"), "fasta"):
        r = role.get(rec.id, "Unknown")
        rep_rows.append([tag, rec.id, r, len(rec.seq)])
        seq = str(rec.seq).upper()
        full.write(f">{tag}|{rec.id}\n{seq}\n")
        all_genomes.write(f">{tag}|{rec.id}\n{seq}\n")
        if r == "Chromosome":
            chrom.write(f">{tag}|{rec.id}\n{seq}\n")
    chrom.close(); full.close()

    prot2cds = {}
    for line in gzip.open(gff, "rt"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        pid = a.get("protein_id")
        locus = a.get("locus_tag", "")
        pseudo = "pseudo" in a and a["pseudo"] == "true"
        row = dict(strain=tag, seqid=f[0], replicon=role.get(f[0], "Unknown"),
                   start=int(f[3]), end=int(f[4]), strand=f[6], locus_tag=locus,
                   protein_id=pid or "", gene=a.get("gene", ""),
                   product=a.get("product", "").replace("%2C", ","), pseudo=pseudo)
        if locus and locus in prot2cds:  # multi-part CDS (same locus_tag): extend
            prev = prot2cds[locus]
            prev["start"] = min(prev["start"], row["start"]); prev["end"] = max(prev["end"], row["end"])
            continue
        cds_rows.append(row)
        if locus:
            prot2cds[locus] = row
    # write proteins keyed by strain|locus_tag (WP_ ids are non-redundant; paralogs may share one)
    seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(gzip.open(faa, "rt"), "fasta")}
    for row in cds_rows:
        if row["strain"] == tag and row["protein_id"] in seqs:
            all_prot.write(f">{tag}|{row['locus_tag']}\n{seqs[row['protein_id']]}\n")
all_prot.close(); all_genomes.close()

import pandas as pd
pd.DataFrame(cds_rows).to_csv(os.path.join(WORK, "cds_table.tsv"), sep="\t", index=False)
pd.DataFrame(rep_rows, columns=["strain", "seqid", "role", "length"]).to_csv(
    os.path.join(WORK, "replicons.tsv"), sep="\t", index=False)
pd.DataFrame(strains).to_csv(os.path.join(WORK, "strains.tsv"), sep="\t", index=False)
print(pd.DataFrame(cds_rows).groupby(["strain", "replicon"]).size())
