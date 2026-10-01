# R version of the dPCR design pipeline (teaching copy)

This folder is a line-by-line **R / tidyverse** translation of the Python pipeline in `primer_design/scripts/`, made for teaching. The Python scripts are untouched. The R version writes only inside `primer_design/R/`: its own `work/`, `results/` and `figures/` (figure files end in `_R`).

**It reproduces the Python results exactly:** the same 1,078 primer/probe designs, the same specificity table, the same final 7-plex, and the same `verification.txt` results (only the column spacing of its enzyme table differs). See *Validation* below.

## Setup and running

```bash
conda env create -f envs/dpcr-design-r.yml      # R 4.4, tidyverse, Biostrings, circlize, patchwork,
                                                # primer3 (primer3_core + ntthal), mmseqs2, bowtie, blast
cd primer_design/R
bsub < submit_dpcr_design_R.sh                  # LSF, all steps; edit STEPS= to run a subset
# or interactively, one step at a time:
conda activate dpcr-design-r
Rscript scripts/00_prep.R
```

Scripts must be run from `primer_design/R/`, because every path is relative to it. Only step 01 (mmseqs2) really needs the cluster. Every other step runs in 1–5 min on one core.

## Scripts

| Step | Script | What it teaches |
|---|---|---|
| — | `helpers.R` | Shared settings; calling command-line tools from R (`system()`); `ntthal` for ΔG; the bowtie1 in-silico PCR |
| 00 | `00_prep.R` | Reading GFF/FASTA with `readr` + Biostrings; parsing GFF attributes with `stringr` |
| 01 | `01_homology.R` | All-vs-all protein search (mmseqs2), then counting homologs per genome with `dplyr` |
| 01a | `01a_gff_text.R` | The "just parse the annotation" route: unique gene names/products via `count()` + `complete()` |
| 02 | `02_sequence_unique.R` | Vectorised k-mer counting: each 18-mer becomes a number, `duplicated()` finds repeats; `rle()` for unique runs; BLAST single-copy check |
| 03 | `03_design.R` | Driving `primer3_core` through its Boulder-IO text format; ranking with `arrange()` |
| 04 | `04_specificity.R` | Exhaustive in-silico PCR: bowtie1 sites + a `dplyr` join of + and − sites |
| 05 | `05_multiplex.R` | Scoring all 6⁷ = 279,936 combinations with vectorised matrix indexing; picking the 4 + 3 well split |
| 06 | `06_verify_report.R` | Pooled PCR, knock-out test, restriction-site search with `vcountPattern()`, the order sheet |
| 07 | `07_figures.R` | Figures 1–4 with ggplot2 + patchwork; circular genome maps with circlize |

## Outputs

- `results/`: `final_multiplex.tsv`, `oligos.tsv` (order sheet), `amplicons.fasta` (gBlock standards), `specificity_report.tsv`, `candidates_pool.tsv`, `candidates_all.tsv`, `dimer_matrix.tsv`, `dimer_heatmap_R.png`, `verification.txt`
- `figures/` (PNG + PDF):
  - `fig1_workflow_R`: the pipeline, with the number of genes and assays at each step
  - `fig2_chromosome_maps_R`: primary chromosomes (dnaA at the top) showing unique DNA, passing assays and the chosen target
  - `fig3_specificity_R`: an assay × genome matrix of products and closest off-target site, plus the knock-out test
  - `fig4_multiplex_compatibility_R`: the oligo ΔG matrix of the final pool, and the score distribution of all candidate 7-plexes

## Python → R equivalents worth knowing

| Python | R used here | Note |
|---|---|---|
| `primer3.bindings.design_primers` | `primer3_core` binary, Boulder-IO text in/out | Identical results (same primer3 engine) |
| `primer3.calc_heterodimer` / `calc_end_stability` | `ntthal -a ANY` / `-a END1` | Identical ΔG. Called once per pair, because batch mode (`-i`) prints nothing when no structure forms, so its output could not be matched back to the pairs |
| jellyfish k-mer counting | base R: 18-mer → base-4 number, `duplicated()` | Same masks; about 3 GB RAM for all 7 genomes |
| numpy boolean masks, `np.diff` runs | logical vectors, `rle()`, `cumsum()` | |
| `Bio.Restriction` | `Biostrings::vcountPattern` on both strands, replicons treated as circular | |

## Pitfalls met during the translation (useful for teaching)

1. **Masked functions.** Biostrings/IRanges also define `slice()`, `rename()` and `desc()`. `helpers.R` therefore loads Biostrings *before* the tidyverse so that dplyr's versions win.
2. **Empty-vector recycling.** `paste0(character(0), ",", integer(0), collapse = " ")` returns `","`, not `""`. A fully unique gene then sent an illegal excluded region to primer3 and was silently skipped. That changed which genes got designed until it was caught by comparing with Python.
3. **Non-palindromic restriction sites.** BsaI (GGTCTC) must be searched on both strands. Searching one strand undercounted the sites.
4. **Name shadowing inside `tibble()`.** A column named `probe` hides a local variable `probe` for the arguments that follow it.

## Validation against the Python pipeline

| Output | Result |
|---|---|
| CDS table / proteins | 33,170 / 32,804, identical |
| Homology and text classes | identical counts per strain |
| Unique fraction per chromosome, candidate genes | identical (5,412 genes; same unique-run lengths and BLAST counts) |
| Designs | all 1,078 primer + probe sets have identical sequences |
| Specificity | identical pass counts per strain |
| Final 7-plex | same 7 assays, same ΔG (−5.64 / −4.17 kcal/mol), same 4 + 3 split |
| `verification.txt`, `oligos.tsv`, `amplicons.fasta` | same content; `amplicons.fasta` byte-identical, the other two differ only in formatting (column spacing; `40` vs `40.0`) |
