# 7SynCom_Kolter: multiplex dPCR assays for the maize root SynCom

Goal: design primer + hydrolysis-probe sets that can be multiplexed in **digital PCR**, one per strain. Each set targets a **single-copy chromosomal locus found only in that strain**. They are used to measure the absolute DNA abundance of each member in mixed samples of these 7 strains only (DNA, not expression).

Reference: Niu B, Paulson JN, Zheng X, Kolter R (2017). *Simplified and representative bacterial community of maize roots.* PNAS 114:E2450–E2459. Genomes: BioProject **PRJNA357031**.

## Status (2026-09-30): design complete and verified in silico

The seven assays are in `primer_design/results/oligos.tsv` (order sheet) and `final_multiplex.tsv` (all details).

| Strain | Species (NCBI) | Target locus | Product | Route | rel_ori | Amplicon | Well (4+3) |
|---|---|---|---|---|---|---|---|
| AA1 | *S. maltophilia* | BTU49_RS02115 *rplP* | 50S ribosomal protein L16 | core 1:1, both routes | 0.39 | 101 bp | W2 |
| AA2 | *B. pituitosa* | BUE85_RS12620 *chpT* | histidine phosphotransferase ChpT | unique, both routes | 0.64 | 134 bp | W1 |
| AA3 | *C. pusillum* | BUE88_RS14975 | hypothetical protein | unique (homology only) | 0.37 | 143 bp | W1 |
| AA4 | *E. ludwigii* | BUE86_RS06915 | ATP-dependent endonuclease | unique, both routes | 0.65 | 76 bp | W1 |
| AA5 | *C. indologenes* | BUE84_RS01130 *bshC* | BshC | unique, both routes | 0.56 | 107 bp | W1 |
| AA6 | *H. robiniae* | BUQ72_RS22945 | AAA ATPase | unique, both routes | 0.58 | 88 bp | W2 |
| AA7 | *P. putida* | BUQ73_RS07385 | M14 carboxypeptidase | unique, both routes | 0.63 | 135 bp | W2 |

`rel_ori` is the target's distance from the replication origin (dnaA): 0 = origin, 1 = terminus. All seven targets are mid-replichore (0.3–0.7).

### Oligos to order (best set, 21 oligos)

All sequences are 5′→3′. The same list is in `primer_design/results/oligos.tsv`. Tm was calculated at 50 mM K⁺, 3.8 mM Mg²⁺, 800 nM primer and 400 nM probe.

| Name | Strain | Type | Sequence (5′→3′) | nt | Tm °C | GC % |
|---|---|---|---|---|---|---|
| AA1_BTU49_RS02115_F | AA1 | F | `AAGAAGCCCATCGAAGTT` | 18 | 60.4 | 44.4 |
| AA1_BTU49_RS02115_R | AA1 | R | `ACACCCTCGATTTCATAGAT` | 20 | 59.9 | 40.0 |
| AA1_BTU49_RS02115_P | AA1 | probe | `CCCAGTATTCCACGTTGCCCTTAC` | 24 | 68.3 | 54.2 |
| AA2_BUE85_RS12620_F | AA2 | F | `CAATCCACGTTTTGTCATCA` | 20 | 60.7 | 40.0 |
| AA2_BUE85_RS12620_R | AA2 | R | `GAAGCAGTGTGTAATATGGC` | 20 | 60.6 | 45.0 |
| AA2_BUE85_RS12620_P | AA2 | probe | `CGTACCGCCAAAGTTCCTCGAA` | 22 | 68.1 | 54.5 |
| AA3_BUE88_RS14975_F | AA3 | F | `TTGTAGTCCATAACGACCTT` | 20 | 60.0 | 40.0 |
| AA3_BUE88_RS14975_R | AA3 | R | `CTGCATCTCTACAACTCGAT` | 20 | 60.9 | 45.0 |
| AA3_BUE88_RS14975_P | AA3 | probe | `CGGAATCGACGTCAAGACACTCAT` | 24 | 67.7 | 50.0 |
| AA4_BUE86_RS06915_F | AA4 | F | `GGTACATAAAGTGCTCCATG` | 20 | 59.9 | 45.0 |
| AA4_BUE86_RS06915_R | AA4 | R | `CTAAATAATGAGCGCGAAGA` | 20 | 59.7 | 40.0 |
| AA4_BUE86_RS06915_P | AA4 | probe | `AGAGCGCGACCACTTAACCATG` | 22 | 68.2 | 54.5 |
| AA5_BUE84_RS01130_F | AA5 | F | `GGAACGACTTGAAAATCTGT` | 20 | 60.0 | 40.0 |
| AA5_BUE84_RS01130_R | AA5 | R | `CAAGCCACGAATAACCATAA` | 20 | 59.9 | 40.0 |
| AA5_BUE84_RS01130_P | AA5 | probe | `ACGCTAAAATTATATACTCTCTCCTGCCA` | 29 | 66.5 | 37.9 |
| AA6_BUQ72_RS22945_F | AA6 | F | `TCGTCTCTCAAGGATTCTTT` | 20 | 60.0 | 40.0 |
| AA6_BUQ72_RS22945_R | AA6 | R | `TTACCGGAACAAGAAGATCT` | 20 | 60.0 | 40.0 |
| AA6_BUQ72_RS22945_P | AA6 | probe | `AACGTTGACTCTTTCTCCTCCGAC` | 24 | 67.4 | 50.0 |
| AA7_BUQ73_RS07385_F | AA7 | F | `GTATAACCATGCACTGTCTG` | 20 | 59.9 | 45.0 |
| AA7_BUQ73_RS07385_R | AA7 | R | `GGAAATTCCCCATGTGTTTA` | 20 | 59.7 | 40.0 |
| AA7_BUQ73_RS07385_P | AA7 | probe | `TTGTGCAGTGTGACCAAGGACTTC` | 24 | 68.4 | 50.0 |

- **Primers:** standard desalted.
- **Probes:** hydrolysis probes with a 5′ reporter dye and a 3′ quencher. Choose the dyes once the platform and channel layout are set, following the W1/W2 grouping above. Double-quenched probes are advisable because several probes are ≥ 24 nt, and the AA5 probe is 29 nt.
- **Quantification standards:** order the 7 amplicons with ±20 bp flanks as gBlocks, from `primer_design/results/amplicons.fasta`.
- **Backups:** other assays per strain that passed all checks are in `primer_design/results/candidates_pool.tsv`.

In-silico verification (`results/verification.txt`):
- **Pooled in-silico PCR:** all 14 primers, every primer combination, ≤ 3 mismatches, over the 7 full genomes including plasmids. Result: **exactly the 7 designed products** (all perfect matches).
- **Knock-out test:** masking each target removes only that strain's product (7/7 PASS).
- **Primer Tm:** 59.7–60.9 °C.
- **Oligo interactions:** the worst interaction between two different sets is −5.64 kcal/mol (3′-anchored: −4.48). The worst oligo pair in the whole pool is −7.42 kcal/mol, which includes pairs within the same set.
- **Probe cross-talk:** each probe has ≥ 5 mismatches to every other amplicon.
- **Restriction enzymes** for fragmenting gDNA before dPCR that cut none of the 7 amplicons: **EcoRI, HindIII, PvuII, BamHI, BsaI, XbaI**. HaeIII and MspI cut 5 of the amplicons, so avoid them.

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

- Environment: conda env `dpcr-design` (`primer_design/envs/dpcr-design.yml`: primer3-py, mmseqs2, BLAST+, bowtie1, jellyfish, seqkit, biopython, pandas).
- Running time: steps 03–06 take about 7 min on one core; only 01 (mmseqs2) and 02 (k-mer scan) benefit from LSF.
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
4. **02 sequence-aware route**: canonical 18-mers counted over all replicons of all 7 genomes. A position is kept only if every 18-mer covering it occurs once in the whole set, which gives 65–93 % of each chromosome. A candidate gene from either route needs:
   - a unique stretch ≥ 150 bp;
   - a single BLAST hit in its own genome (≥ 80 % id over ≥ 50 bp).

   6,367 unique intergenic windows are kept as a fallback.
5. **03 design**: primer3 designs within the unique sequence only.
   - **Target location:** only the dnaA-bearing chromosome is used. This matters for AA2: its secondary chromosomes are not guaranteed to be present 1:1 with the primary one.
   - **Origin distance:** `rel_ori` = distance from dnaA, from 0 (origin) to 1 (terminus). Mid-replichore targets (0.3–0.7) are preferred, which limits the origin-to-terminus copy-number bias in growing cells.
   - **Gene quota:** per strain, 40 core single-copy genes and 40 strain-unique genes. Rank order: unique by both routes, core by both routes, unique by homology only, core by one route, then text-unique only.
   - Amplicon 70–150 bp. Primers 18–25 nt, Tm 58–62 °C. Probe Tm 66–71 °C, on either strand, with no 5′ G.
   - Buffer: 50 mM K⁺, 3.8 mM Mg²⁺, 800 nM primers, 400 nM probe. Primer GC limits follow each genome's GC.
6. **04 specificity**: `insilico_pcr.py` uses **bowtie1** (`-v 3 -a`) to list *every* binding site with ≤ 3 mismatches of every oligo on all 7 genomes, including plasmids. Any + site paired with a downstream − site counts as a product, including F+F and R+R.
   - A set needs exactly one product (≤ 3 kb), at the designed locus.
   - A lone primer site is *dangerous*, and rejects the set, if it has ≤ 1 mismatch, or 2 mismatches with an intact 3′ pentamer.
   - Sites with 3 mismatches are random background (about 13 per primer in 36 Mb); they are only counted and used for ranking.
   - *Note:* `seqkit amplicon` was used first but was dropped. With mismatches allowed it returns only one (the longest) product per primer pair and strand, which hid true products behind spurious Mb-long ones.
7. **05 multiplex**: the top 6 passing genes per strain, and all 6⁷ = 279,936 combinations.
   - Objective: best worst-case heterodimer ΔG *between sets*, then best 3′-anchored ΔG, then smallest Tm spread. Within-set dimers were already constrained by primer3.
   - Checks: exhaustive in-silico PCR of the pooled primers, and each probe needs ≥ 5 mismatches to every other amplicon.
   - It also proposes a 4 + 3 well split. The step stops with an error if any strain has no passing set.
8. **06 verification**:
   - Pooled in-silico PCR must give exactly 7 products.
   - Knock-out test per strain (drop the primer sites inside the target, as if masked with Ns).
   - Restriction-enzyme compatibility of the amplicons.
   - Writes `amplicons.fasta` (± 20 bp, as gBlock standards), `oligos.tsv` and the dimer heatmap.

## Review follow-up

A separate agent reviewed the workflow. These findings were fixed:
- seqkit's non-exhaustive in-silico PCR, replaced by bowtie1 pairing;
- walltime;
- BLAST HSP truncation, since the BLAST step was removed;
- core-first ranking hid the strain-unique genes, now replaced by a quota;
- AA2 secondary replicons, now dnaA chromosome only;
- the within-set diagonal distorted the multiplex objective;
- a 6-plex could be written silently;
- no ori-position control, now `rel_ori`;
- mismatch thresholds were inconsistent;
- the text route could never be ranked into the designs;
- the empty-product key;
- O(n²) prep.

Still open:
- **No host or contaminant screen.** The 14 primers have not been screened against *Zea mays* (nuclear, chloroplast, mitochondrial) or common reagent contaminants. Do this if the assays will be used on root DNA.
- **Cached files are not tied to their parameters.** `allvsall.m8`, `k18_repeated.txt` and the BLAST/bowtie indexes are reused whenever they exist. After changing `EVALUE`, `MIN_COV` or `K`, delete `primer_design/work/`.
- **AA4 target may sit in an unstable region.** BUE86_RS06915 is an "ATP-dependent endonuclease" (an OLD-family gene; these often sit in defence islands). It is not near an annotated mobile element, but if you want a housekeeping target instead, take the next AA4 set from `results/candidates_pool.tsv`.
- **Dye and amplitude plan.** Assign dyes once the platform is chosen. For a 7-in-1 well, amplitude multiplexing needs different probe concentrations per assay, tuned experimentally.

## Still to do

- Wet-lab validation:
  - efficiency of each assay alone, and alone vs multiplexed, on the gDNA of each pure strain;
  - a cross-reactivity matrix (each strain's DNA against all 7 assays);
  - calibration against gBlock standards (`results/amplicons.fasta`).
- Digest the gDNA before partitioning with one of the compatible enzymes listed above.
