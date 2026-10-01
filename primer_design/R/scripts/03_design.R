# Step 03 - design primer + probe sets with primer3 (primer3_core, Boulder-IO format)
#
# * Targets only on the chromosome carrying dnaA (primary chromosome; the secondary
#   chromosomes of AA2 need not be present 1:1 with it).
# * rel_ori = distance from dnaA (~ origin of replication): 0 = origin, 1 = terminus.
#   Mid-replichore genes (0.3-0.7) are preferred: in growing cells, genes near the
#   origin are present in more copies than genes near the terminus.
# * Per strain, 40 "core" genes (single copy in every genome) and 40 "unique" genes
#   (no homolog in the other genomes) are designed; steps 04-05 decide.
# * Non-unique positions (step 02) are excluded, so primers and probe sit in DNA
#   found only once in the 7 genomes.
source("scripts/helpers.R")

GENES_PER_GROUP <- 40
SETS_PER_GENE   <- 3

primer3_settings <- c(
  PRIMER_TASK = "generic",
  PRIMER_PICK_LEFT_PRIMER = 1, PRIMER_PICK_RIGHT_PRIMER = 1, PRIMER_PICK_INTERNAL_OLIGO = 1,
  PRIMER_PRODUCT_SIZE_RANGE = "70-150",
  PRIMER_OPT_SIZE = 20, PRIMER_MIN_SIZE = 18, PRIMER_MAX_SIZE = 25,
  PRIMER_OPT_TM = 60, PRIMER_MIN_TM = 58, PRIMER_MAX_TM = 62, PRIMER_PAIR_MAX_DIFF_TM = 1.5,
  PRIMER_MAX_POLY_X = 4, PRIMER_MAX_END_GC = 3, PRIMER_MAX_END_STABILITY = 9,
  PRIMER_MAX_SELF_ANY_TH = 35, PRIMER_MAX_SELF_END_TH = 30, PRIMER_PAIR_MAX_COMPL_ANY_TH = 35,
  PRIMER_PAIR_MAX_COMPL_END_TH = 30, PRIMER_MAX_HAIRPIN_TH = 35,
  PRIMER_INTERNAL_OPT_SIZE = 24, PRIMER_INTERNAL_MIN_SIZE = 18, PRIMER_INTERNAL_MAX_SIZE = 30,
  PRIMER_INTERNAL_OPT_TM = 68, PRIMER_INTERNAL_MIN_TM = 66, PRIMER_INTERNAL_MAX_TM = 71,
  PRIMER_INTERNAL_MAX_POLY_X = 4, PRIMER_INTERNAL_MAX_SELF_ANY_TH = 40, PRIMER_INTERNAL_MAX_HAIRPIN_TH = 40,
  PRIMER_SALT_MONOVALENT = BUFFER$mv, PRIMER_SALT_DIVALENT = BUFFER$dv,
  PRIMER_DNTP_CONC = BUFFER$dntp, PRIMER_DNA_CONC = BUFFER$primer_nM,
  PRIMER_INTERNAL_SALT_MONOVALENT = BUFFER$mv, PRIMER_INTERNAL_SALT_DIVALENT = BUFFER$dv,
  PRIMER_INTERNAL_DNTP_CONC = BUFFER$dntp, PRIMER_INTERNAL_DNA_CONC = BUFFER$probe_nM,
  PRIMER_NUM_RETURN = 10
)

# primer GC window follows the genome composition (36 % GC Chryseobacterium ... 71 % Curtobacterium)
gc_limits <- function(genome_gc) c(min = if (genome_gc < 0.45) 25 else 35, max = if (genome_gc > 0.62) 72 else 65)

# non-unique stretches of a gene -> "start,length start,length ..." (0-based, primer3 format)
excluded_regions <- function(unique) {
  r <- rle(!unique)
  if (!any(r$values)) return("")          # gene entirely unique: nothing to exclude
  ends <- cumsum(r$lengths); starts <- ends - r$lengths
  paste0(starts[r$values], ",", r$lengths[r$values], collapse = " ")
}

# a probe can be on either strand: avoid a 5' G (quenches the dye), prefer more C than G
orient_probe <- function(p) {
  opts <- tibble(seq = c(p, revcomp(p)), strand = c("+", "-")) |>
    mutate(no5G = !str_starts(seq, "G"), cg = str_count(seq, "C") - str_count(seq, "G"))
  if (any(opts$no5G)) opts <- filter(opts, no5G)
  opts |> arrange(desc(cg)) |> slice(1)
}

# ---- candidate genes, ranked ------------------------------------------------------------
masks <- readRDS(file.path(WORK, "unique_masks.rds"))
chrom <- map(set_names(list.files(WORK, "\\.chrom\\.fna$", full.names = TRUE)), readDNAStringSet) |>
  unname() |> do.call(what = c)
genes <- read_tsv(file.path(WORK, "cds_table.tsv"))
dnaA <- genes |> filter(gene == "dnaA") |> distinct(strain, .keep_all = TRUE) |>
  transmute(strain, ori_seqid = seqid, ori_pos = start)

cand <- read_tsv(file.path(WORK, "candidates_genes.tsv")) |>
  filter(single_copy_nt) |>
  inner_join(dnaA, by = "strain") |>
  filter(seqid == ori_seqid) |>
  mutate(chrom_len = width(chrom[chrom_id]),
         dist = abs((start + end) / 2 - ori_pos),
         rel_ori = pmin(dist, chrom_len - dist) / (chrom_len / 2),
         mid_replichore = between(rel_ori, 0.3, 0.7),
         class_rank = case_when(
           class == "A_unique" & text_class == "text_unique"             ~ 0,  # unique by both routes
           class == "B_universal_1to1" & text_class == "text_core_1x"    ~ 1,  # core by both routes
           class == "A_unique"                                           ~ 2,
           class == "B_universal_1to1" | text_class == "text_core_1x"    ~ 3,
           TRUE                                                          ~ 4), # text_unique only
         group = if_else(class_rank %in% c(1, 3), "core", "unique")) |>
  arrange(strain, desc(mid_replichore), class_rank, desc(uniq_run_len), best_other_pident)

genome_gc <- tibble(chrom_id = names(chrom), gc = letterFrequency(chrom, "GC", as.prob = TRUE)[, 1],
                    len = width(chrom)) |>
  mutate(strain = word(chrom_id, 1, sep = fixed("|"))) |>
  group_by(strain) |> summarise(gc = weighted.mean(gc, len))

# ---- run primer3 on all candidates of a strain in one batch -------------------------------
run_primer3 <- function(genes, gc_lim) {
  settings <- c(primer3_settings, PRIMER_MIN_GC = gc_lim[["min"]], PRIMER_MAX_GC = gc_lim[["max"]],
                PRIMER_INTERNAL_MIN_GC = gc_lim[["min"]], PRIMER_INTERNAL_MAX_GC = gc_lim[["max"]])
  settings_txt <- paste0(names(settings), "=", settings, collapse = "\n")
  records <- pmap_chr(genes, \(id, template, excluded, ...)
    paste0("SEQUENCE_ID=", id, "\nSEQUENCE_TEMPLATE=", template,
           if (excluded != "") paste0("\nSEQUENCE_EXCLUDED_REGION=", excluded) else "",
           "\n", settings_txt, "\n="))
  inp <- tempfile(fileext = ".bio")
  writeLines(records, inp)
  out <- system(paste("primer3_core <", inp), intern = TRUE)
  # parse Boulder-IO: records end with "=", lines are KEY=VALUE
  rec <- cumsum(c(0, head(out == "=", -1)))
  tibble(rec = rec, line = out) |>
    filter(line != "=") |>
    separate_wider_delim(line, "=", names = c("key", "value"), too_many = "merge") |>
    group_by(rec) |>
    summarise(res = list(set_names(value, key)))
}

designs <- list()
for (s in unique(cand$strain)) {
  lim <- gc_limits(genome_gc$gc[genome_gc$strain == s])
  g <- cand |> filter(strain == s) |>
    mutate(template = as.character(subseq(chrom[chrom_id], start, end)),
           excluded = pmap_chr(list(chrom_id, start, end), \(cid, a, b) excluded_regions(masks[[cid]][a:b])),
           id = paste0(s, "_", locus_tag))

  # genes are sent to primer3 in chunks of 100 (in rank order) until both quotas are full
  n_genes <- c(core = 0, unique = 0)
  res <- list()
  for (i in seq_len(nrow(g))) {
    if (all(n_genes >= GENES_PER_GROUP)) break
    r <- g[i, ]
    if (n_genes[[r$group]] >= GENES_PER_GROUP) next
    if (i > length(res)) {
      chunk <- i:min(nrow(g), i + 99)
      res[chunk] <- run_primer3(g[chunk, ] |> select(id, template, excluded), lim)$res
    }
    p <- res[[i]]
    if (!is.na(p["PRIMER_ERROR"])) next
    n <- as.integer(p["PRIMER_PAIR_NUM_RETURNED"])
    kept <- 0; used <- integer(0)
    for (k in seq_len(n) - 1) {
      key <- function(what) p[[sprintf(what, k)]]
      left_start <- as.integer(str_split_1(key("PRIMER_LEFT_%d"), ",")[1])
      right_end  <- as.integer(str_split_1(key("PRIMER_RIGHT_%d"), ",")[1])   # 3'-most base, 0-based
      if (any(abs(left_start - used) < 40)) next                              # diversify amplicons
      pr <- orient_probe(key("PRIMER_INTERNAL_%d_SEQUENCE"))
      if (str_starts(pr$seq, "G")) next
      amp <- substr(r$template, left_start + 1, right_end + 1)
      designs[[length(designs) + 1]] <- tibble(
        set_id = paste0(r$id, "_", k), strain = s, route = r$route, class_rank = r$class_rank,
        group = r$group, rel_ori = round(r$rel_ori, 3), locus_tag = r$locus_tag,
        gene = replace_na(r$gene, ""), product = r$product, seqid = r$seqid,
        gene_start = r$start, gene_end = r$end, gene_strand = r$strand,
        amp_start = r$start + left_start, amp_end = r$start + left_start + nchar(amp) - 1, amp_len = nchar(amp),
        fwd = key("PRIMER_LEFT_%d_SEQUENCE"), rev = key("PRIMER_RIGHT_%d_SEQUENCE"),
        probe = pr$seq, probe_strand = pr$strand,
        fwd_tm = as.numeric(key("PRIMER_LEFT_%d_TM")), rev_tm = as.numeric(key("PRIMER_RIGHT_%d_TM")),
        probe_tm = as.numeric(key("PRIMER_INTERNAL_%d_TM")),
        fwd_gc = as.numeric(key("PRIMER_LEFT_%d_GC_PERCENT")), rev_gc = as.numeric(key("PRIMER_RIGHT_%d_GC_PERCENT")),
        probe_gc = as.numeric(key("PRIMER_INTERNAL_%d_GC_PERCENT")), amp_gc = gc_pct(amp),
        penalty = as.numeric(key("PRIMER_PAIR_%d_PENALTY")), best_other_pident = r$best_other_pident,
        amplicon = amp)
      used <- c(used, left_start); kept <- kept + 1
      if (kept >= SETS_PER_GENE) break
    }
    if (kept > 0) n_genes[[r$group]] <- n_genes[[r$group]] + 1
  }
  cat(sprintf("%s: GC %.2f, primer GC %g-%g %%, designed on %d core + %d unique genes\n",
              s, genome_gc$gc[genome_gc$strain == s], lim[["min"]], lim[["max"]], n_genes[["core"]], n_genes[["unique"]]))
}
designs <- bind_rows(designs)
write_tsv(designs, file.path(WORK, "designs_raw.tsv"))
designs |> count(strain, class_rank) |> pivot_wider(names_from = class_rank, values_from = n, values_fill = 0) |> print()
