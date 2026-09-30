#!/usr/bin/env python
"""Annotation route by pure text parsing of the RefSeq/PGAP GFF features.

Key per CDS = gene symbol (gene=) if present, else the normalised product string.
Uninformative keys (hypothetical / family / domain-containing / DUF ...) are dropped.

text_unique   : key occurs exactly once in the target genome and in no other genome.
text_core_1x  : key occurs exactly once in each of the 7 genomes (named single-copy
                housekeeping genes, e.g. gyrB, rpoB, recA); the target locus then
                needs strain-specific sequence (checked in step 02).
No sequence comparison is done here; step 02 verifies nucleotide uniqueness.
"""
import os, re
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "..", "work")

cds = pd.read_csv(os.path.join(WORK, "cds_table.tsv"), sep="\t")
cds["id"] = cds.strain + "|" + cds.locus_tag
strains = sorted(cds.strain.unique())

generic = re.compile(r"hypothetical|family|domain-containing|DUF\d|uncharacteri[sz]ed|"
                     r"putative|unknown|^protein$|YbaB|-like", re.I)


def key(r):
    if isinstance(r.gene, str) and r.gene:
        return "gene:" + re.sub(r"_\d+$", "", r.gene)   # PGAP adds _1/_2 to duplicated symbols
    p = r["product"] if isinstance(r["product"], str) else ""
    if not p or generic.search(p):
        return ""
    return "product:" + p.lower().strip()


cds["text_key"] = cds.apply(key, axis=1)
# a gene symbol and a product can describe the same gene in another genome where the
# symbol is missing -> also count the product string for named genes
cds["prod_key"] = "product:" + cds["product"].fillna("").str.lower().str.strip()

k = cds[cds.text_key != ""]
by_key = k.groupby(["text_key", "strain"]).size().unstack(fill_value=0).reindex(columns=strains, fill_value=0)
by_prod = cds.groupby(["prod_key", "strain"]).size().unstack(fill_value=0).reindex(columns=strains, fill_value=0)

cls = []
for _, r in cds.iterrows():
    if not r.text_key:
        cls.append(""); continue
    row = by_key.loc[r.text_key]
    prow = by_prod.loc[r.prod_key] if not generic.search(str(r["product"])) else None
    others = [s for s in strains if s != r.strain]
    if row[r.strain] == 1 and (row[others] == 0).all() and (prow is None or (prow[others] == 0).all()):
        cls.append("text_unique")
    elif (row == 1).all():
        cls.append("text_core_1x")
    else:
        cls.append("")
cds["text_class"] = cls
cds[["id", "text_key", "text_class"]].to_csv(os.path.join(WORK, "cds_text.tsv"), sep="\t", index=False)

print(cds.groupby(["strain", "text_class"]).size().unstack(fill_value=0))
print("\nexamples of text_core_1x keys:", ", ".join(sorted(cds[cds.text_class == "text_core_1x"].text_key.unique())[:25]))
