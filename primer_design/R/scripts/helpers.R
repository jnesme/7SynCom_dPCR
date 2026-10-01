# Shared settings and helper functions, sourced by every step.
# All scripts are run from primer_design/R/  (e.g. `Rscript scripts/03_design.R`).

suppressPackageStartupMessages({
  library(Biostrings)   # loaded first so that dplyr's slice/rename/desc/... are not masked
  library(tidyverse)
})
options(dplyr.summarise.inform = FALSE, readr.show_col_types = FALSE)

GENOMES <- "../../genomes"          # downloaded NCBI assemblies
WORK    <- "work"                   # intermediate files
RES     <- "results"                # final tables
FIG     <- "figures"
dir.create(WORK, showWarnings = FALSE)
dir.create(RES, showWarnings = FALSE)
dir.create(FIG, showWarnings = FALSE)
THREADS <- Sys.getenv("THREADS", "1")

# PCR buffer used for every thermodynamic calculation (typical probe dPCR master mix)
BUFFER <- list(mv = 50, dv = 3.8, dntp = 0.8, primer_nM = 800, probe_nM = 400, temp = 37)

# ---------------------------------------------------------------------------
# small utilities
# ---------------------------------------------------------------------------
run <- function(cmd) {
  # run a shell command, stop on error
  status <- system(cmd)
  if (status != 0) stop("command failed: ", cmd)
}

revcomp <- function(x) as.character(reverseComplement(DNAStringSet(x)))

gc_pct <- function(x) 100 * str_count(x, "[GC]") / nchar(x)

# ---------------------------------------------------------------------------
# thermodynamics with primer3's `ntthal` (identical engine to primer3-py)
# ---------------------------------------------------------------------------
# dG (kcal/mol) for many oligo pairs.
#   mode "ANY"  : any duplex between a and b
#   mode "END1" : 3' end of a annealed on b (can the 3' end of a be extended?)
# ntthal is called once per pair (in parallel with xargs): its batch mode (-i) prints
# nothing when no structure forms, so its output could not be matched to the input.
# No structure = dG 0, as in primer3-py.
ntthal_dg <- function(a, b, mode = "ANY", dna_nM = BUFFER$primer_nM) {
  if (length(a) == 0) return(numeric(0))
  inp <- tempfile(fileext = ".txt")
  writeLines(paste(seq_along(a), a, b), inp)
  cmd <- sprintf(paste("xargs -P %s -L 1 sh -c 'echo \"$0 $(ntthal -a %s -mv %s -dv %s -n %s -d %s -t %s -s1 $1 -s2 $2 2>/dev/null",
                       "| grep -o \"dG = [-0-9.]*\" | cut -c6-)\"' < %s"),
                 THREADS, mode, BUFFER$mv, BUFFER$dv, BUFFER$dntp, dna_nM, BUFFER$temp, inp)
  out <- read_delim(I(system(cmd, intern = TRUE)), delim = " ", col_names = c("i", "dg"), col_types = "id")
  dg <- numeric(length(a))
  dg[out$i] <- replace_na(out$dg, 0) / 1000
  dg
}

# ---------------------------------------------------------------------------
# exhaustive in-silico PCR with bowtie1
# ---------------------------------------------------------------------------
# Why not `seqkit amplicon`? With mismatches allowed it reports only one (the
# longest) product per primer pair and strand, which can hide the real product.
# Here every binding site with <= 3 mismatches is listed (bowtie -v 3 -a), and any
# '+' site followed by a '-' site on the same sequence is a potential product.

build_bowtie_index <- function(fasta, prefix) {
  if (!file.exists(paste0(prefix, ".1.ebwt")))
    run(sprintf("bowtie-build --threads %s -q %s %s", THREADS, fasta, prefix))
}

# oligos: named character vector (name = oligo id). Returns one row per site.
find_sites <- function(oligos, index, tag, mm = 3) {
  fa  <- file.path(WORK, paste0(tag, ".oligos.fa"))
  out <- file.path(WORK, paste0(tag, ".bowtie.tsv"))
  writeXStringSet(DNAStringSet(setNames(oligos, names(oligos))), fa)
  run(sprintf("bowtie -f -v %d -a -p %s --quiet -x %s %s > %s", mm, THREADS, index, fa, out))
  if (file.size(out) == 0)
    return(tibble(name = character(), chrom = character(), strand = character(),
                  start = numeric(), end = numeric(), mm = integer(), mm3 = integer()))
  read_tsv(out, col_names = c("name", "strand", "chrom", "start", "seq", "qual", "n_other", "mism"),
           col_types = "ccciccic") |>
    mutate(
      mism = replace_na(mism, ""),
      olen = nchar(oligos[name]),
      end  = start + olen,                         # start is 0-based, end exclusive
      # bowtie1 mismatch offsets count from the 5' end of the oligo as given
      offsets = map(mism, \(m) if (m == "") integer(0) else as.integer(str_extract(str_split_1(m, ","), "^[0-9]+"))),
      mm  = lengths(offsets),
      mm3 = map2_int(offsets, olen, \(o, L) sum(o >= L - 5))   # mismatches in the 3'-terminal 5 nt
    ) |>
    select(name, chrom, strand, start, end, mm, mm3)
}

# Pair each '+' site with every '-' site downstream on the same sequence.
pcr_products <- function(sites, max_len = 3000) {
  plus  <- sites |> filter(strand == "+")
  minus <- sites |> filter(strand == "-")
  inner_join(plus, minus, by = "chrom", suffix = c("_l", "_r"), relationship = "many-to-many") |>
    filter(start_r >= start_l, end_r - start_l <= max_len) |>
    transmute(chrom, start = start_l, end = end_r, len = end_r - start_l,
              left_oligo = name_l, right_oligo = name_r, left_mm = mm_l, right_mm = mm_r)
}
