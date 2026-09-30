# 7SynCom_Kolter: multiplex dPCR assays for the maize root SynCom

Goal: design primer + hydrolysis-probe sets that can be multiplexed in **digital PCR**, one per strain. Each set targets a **single-copy chromosomal locus found only in that strain**. They are used to measure the absolute DNA abundance of each member in mixed samples of these 7 strains only (DNA, not expression).

Reference: Niu B, Paulson JN, Zheng X, Kolter R (2017). *Simplified and representative bacterial community of maize roots.* PNAS 114:E2450–E2459. Genomes: BioProject **PRJNA357031**.

## Status (2026-09-30)

| Step | Script | Status |
|---|---|---|
| Genome download | (manual, see below) | done, all md5 OK |
| 00 prep | `00_prep.py` | done |
| 01 homology classes (mmseqs2 all-vs-all) | `01_annotation_candidates.py` | done (LSF job 29536664) |
| 01a GFF text-parsing classes | `01a_gff_text_unique.py` | done |
| 02 nucleotide uniqueness (18-mers) + single-copy BLAST | `02_sequence_unique.py` | done |
| 03 primer/probe design (primer3) | `03_design.py` | done: 1,031 sets |
| 04 in-silico specificity | `04_specificity.py` | **running** (LSF job 29536679) |
| 05 multiplex selection | `05_multiplex.py` | pending |
| 06 verification + report | `06_verify_report.py` | pending |

The final assays do not exist yet, so `primer_design/results/` is still empty.

## Genomes

`genomes/<GCA accession>/` holds the authors' GenBank submission: sequence only, **no annotation**. `genomes/<GCA>/refseq_<GCF>/` holds NCBI's RefSeq version, which adds PGAP annotation (GFF, GBFF, proteins, CDS). The sequence is identical in both. The table below is in `genomes/assemblies.tsv`. The `.gz` files are git-ignored; re-download them from the FTP paths in that table.

| Strain | Name in paper | Current NCBI name | Assembly | Replicons |
|---|---|---|---|---|
| AA1 | *Stenotrophomonas maltophilia* | same | GCA_002025605.1 | chromosome CP018756 |
| AA2 | *Ochrobactrum pituitosum* | *Brucella pituitosa* | GCA_002025625.1 | 3 chromosomes CP018779/80/82 + plasmid pOAAA2 CP018781 |
| AA3 | *Curtobacterium pusillum* | same | GCA_002025645.1 | chromosome CP018783 + plasmid pCPAA3 CP018784 |
| AA4 | *Enterobacter cloacae* | *Enterobacter ludwigii* | GCA_002025685.1 | chromosome CP018785 |
| AA5 | *Chryseobacterium indologenes* | same | GCA_002025665.1 | chromosome CP018786 |
| AA6 | *Herbaspirillum frisingense* | *Herbaspirillum robiniae* | GCA_002025725.1 | chromosome CP018845 |
| AA7 | *Pseudomonas putida* | same | GCA_002025705.1 | chromosome CP018846 |

## Design pipeline (`primer_design/`)

- Environment: conda env `dpcr-design` (`primer_design/envs/dpcr-design.yml`: primer3-py, mmseqs2, BLAST+, seqkit, jellyfish, biopython, pandas).
- Cluster: run `bsub < primer_design/submit_dpcr_design.sh` (LSF, queue `hpc`, 8 cores). Edit `STEPS=` to run a subset of steps; steps reuse outputs that already exist.

1. **00 prep**: per-strain chromosome/plasmid FASTA, a CDS table from the GFF, and proteins keyed `strain|locus_tag`. Only chromosomes are used as targets, because plasmid copy number is not 1.
2. **01 homology route**: mmseqs2 all-vs-all protein search (s 7.5, e ≤ 1e-5, cov ≥ 0.3).
   - *A_unique* = no homolog in the other 6 genomes and no paralog.
   - *B_universal_1to1* = exactly one copy in every genome.
   - Genes that are mobile-element-like (transposase, integrase, phage, TA system, …) are flagged, together with every gene within 5 kb of them.
3. **01a text-parsing route**: GFF text only. The key is the gene symbol, otherwise the product; generic products (hypothetical, "family protein", "domain-containing", …) are ignored.
   - *text_unique* = the key occurs once, in one genome only.
   - *text_core_1x* = the key occurs exactly once in every genome (146 keys, e.g. gyrB, recA, dnaK, atpD).
   - Caveat: this route over-calls uniqueness when naming is uneven. AA4 has 1,076 text-unique genes because Enterobacter genes carry E. coli symbols that the other genomes' annotations lack. Step 02 filters these.
4. **02 sequence-aware route**: canonical 18-mers counted over all replicons of all 7 genomes. A position is kept only if every 18-mer covering it occurs once in the whole set, which gives 75–90 % of each chromosome. A candidate gene from either route needs:
   - a unique stretch ≥ 150 bp;
   - a single BLAST hit in its own genome (≥ 80 % id over ≥ 50 bp).

   6,367 unique intergenic windows are kept as a fallback.
5. **03 design**: primer3 designs within the unique sequence only.
   - Amplicon 70–150 bp. Primers 18–25 nt, Tm 58–62 °C. Probe Tm 66–71 °C, on either strand, with no 5′ G.
   - Buffer: 50 mM K⁺, 3.8 mM Mg²⁺, 800 nM primers, 400 nM probe.
   - Primer GC limits follow each genome's GC, which ranges from 36 % (AA5) to 71 % (AA3).
   - Ranking: core single-copy genes found by both routes first, then genes that both routes call unique, then those found by one route only.
6. **04 specificity**: `seqkit amplicon` (≤ 3 mismatches per primer) on all 7 genomes, plasmids included; a set must give exactly one product, at the designed locus. blastn-short scores every off-target site of each primer; a site with < 4 mismatches in total and < 2 mismatches in the 3′-terminal 5 nt rejects the set.
7. **05 multiplex**: the top 6 genes per strain, all 6⁷ combinations scored with primer3.
   - Objective: best worst-case heterodimer ΔG, then best 3′-anchored ΔG, then smallest Tm spread.
   - Checks on the chosen pool: in-silico PCR with every primer pairing (so no cross-set products), and every probe needs ≥ 5 mismatches to all other amplicons.
   - It also proposes a 4 + 3 split into two wells for instruments with ≤ 5 colours.
8. **06 verification**:
   - Pooled in-silico PCR on the 7 genomes must give exactly 7 products.
   - Knock-out test per strain: mask its target, and its product must disappear.
   - Also writes `amplicons.fasta` (amplicon ± 20 bp, as gBlock standards), `oligos.tsv` (order sheet), and a dimer heatmap.

## Progress so far

- **Candidates after step 02:** every strain keeps single-copy candidate genes, mostly core genes present once in all genomes, with 16–61 text-core genes per strain.
- **Step 03 designs:** 1,031 primer + probe sets over 80 genes per strain.

## Still to do

- Finish steps 04–06 and review the final 7 assays.
- Choose the platform and dyes. The design is generic, and the channel count decides whether the 7 assays fit in one well or need the 4 + 3 split.
- Wet-lab validation:
  - efficiency of each assay alone, and alone vs multiplexed, on the gDNA of each pure strain;
  - a cross-reactivity matrix (each strain's DNA against all 7 assays);
  - calibration against gBlock standards.
