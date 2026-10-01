# Step 02 - sequence-aware route: which parts of each chromosome are unique DNA?
#
# 1. Every 18-mer of all replicons of all 7 genomes (plasmids included) is turned into
#    a number (A=0, C=1, G=2, T=3 in base 4; 4^18 fits exactly in a double). The two
#    strands are merged by keeping the smaller of k-mer and reverse complement
#    ("canonical" k-mer). A k-mer seen more than once anywhere is "repeated".
#    A position is unique only if no repeated 18-mer covers it: it is neither
#    duplicated in its own genome nor present in any other genome.
# 2. For each candidate gene (from step 01 or 01a), find the longest unique stretch.
# 3. BLAST each gene against its own genome: one hit only = single copy at the DNA level.
# Unique intergenic windows are also written as a fallback.
source("scripts/helpers.R")

K <- 18
MIN_RUN <- 150

# ---- 1. canonical k-mer codes, vectorised --------------------------------------------
base_code <- rep(NA_real_, 256)
base_code[utf8ToInt("A")] <- 0; base_code[utf8ToInt("C")] <- 1
base_code[utf8ToInt("G")] <- 2; base_code[utf8ToInt("T")] <- 3

canonical_kmers <- function(dna) {
  x <- base_code[as.integer(charToRaw(as.character(dna)))]
  n <- length(x); m <- n - K + 1
  fwd <- numeric(m); rev <- numeric(m)
  for (j in 0:(K - 1)) {
    fwd <- fwd * 4 + x[(1 + j):(m + j)]                  # k-mer read 5'->3'
    rev <- rev * 4 + (3 - x[(K - j):(m + K - 1 - j)])     # its reverse complement
  }
  pmin(fwd, rev)                                          # NA if the k-mer contains N
}

genomes <- readDNAStringSet(file.path(WORK, "all_genomes.fna"))
cat("encoding", length(genomes), "replicons ...\n")
codes <- lapply(genomes, canonical_kmers)
all_codes <- unlist(codes, use.names = FALSE)
repeated_all <- duplicated(all_codes) | duplicated(all_codes, fromLast = TRUE) | is.na(all_codes)
repeated <- split(repeated_all, rep(names(codes), lengths(codes)))
rm(all_codes, repeated_all); invisible(gc())

# position p is masked if any of the k-mers starting in [p-K+1, p] is repeated
position_unique <- function(bad_kmer) {
  m <- length(bad_kmer); n <- m + K - 1
  cs <- c(0, cumsum(bad_kmer))
  p <- seq_len(n)
  (cs[pmin(p, m) + 1] - cs[pmax(1, p - K + 1)]) == 0
}

reps <- read_tsv(file.path(WORK, "replicons.tsv"))
chrom_ids <- reps |> filter(role == "Chromosome") |> mutate(id = paste0(strain, "|", seqid)) |> pull(id)
masks <- map(set_names(chrom_ids), \(id) position_unique(repeated[[id]]))
saveRDS(masks, file.path(WORK, "unique_masks.rds"))
tibble(id = names(masks), pct_unique = map_dbl(masks, mean) * 100) |>
  mutate(pct_unique = round(pct_unique, 1)) |> print()

# ---- 2. longest unique run inside each candidate gene -----------------------------------
longest_run <- function(v) {
  r <- rle(v)
  if (!any(r$values)) return(c(start = NA, len = 0))
  ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  k <- which(r$values)[which.max(r$lengths[r$values])]
  c(start = starts[k], len = r$lengths[k])
}

cds <- read_tsv(file.path(WORK, "cds_classified.tsv")) |>
  left_join(read_tsv(file.path(WORK, "cds_text.tsv")), by = "id") |>
  mutate(class = replace_na(class, ""), text_class = replace_na(text_class, ""),
         route = paste0(if_else(class == "", "-", class), "/", if_else(text_class == "", "-", text_class)),
         eligible = (class != "" | text_class != "") & replicon == "Chromosome" & !pseudo &
                    !near_mobile & (end - start + 1 >= 300))

cand <- cds |> filter(eligible) |>
  mutate(chrom_id = paste0(strain, "|", seqid),
         run = pmap(list(chrom_id, start, end), \(cid, s, e) longest_run(masks[[cid]][s:e])),
         unique_frac  = pmap_dbl(list(chrom_id, start, end), \(cid, s, e) mean(masks[[cid]][s:e])),
         uniq_run_len = map_dbl(run, "len"),
         uniq_run_start = start + map_dbl(run, "start") - 1) |>
  select(-run) |>
  filter(uniq_run_len >= MIN_RUN)

# ---- 3. single copy at the DNA level (BLAST against own genome) ----------------------------
cand$n_self_hits <- 0L
for (s in unique(cand$strain)) {
  db <- file.path(WORK, paste0("blastdb_", s))
  if (!file.exists(paste0(db, ".nsq")))
    run(sprintf("makeblastdb -in %s/%s.all.fna -dbtype nucl -out %s > /dev/null", WORK, s, db))
  sub <- which(cand$strain == s)
  g <- readDNAStringSet(file.path(WORK, paste0(s, ".all.fna")))
  q <- subseq(g[cand$chrom_id[sub]], cand$start[sub], cand$end[sub])
  names(q) <- cand$id[sub]
  qf <- file.path(WORK, paste0("cand_genes_", s, ".fna")); writeXStringSet(q, qf)
  out <- file.path(WORK, paste0("cand_genes_", s, ".self.tsv"))
  run(sprintf("blastn -query %s -db %s -outfmt '6 qseqid sseqid pident length' -evalue 1e-5 -max_hsps 50 -num_threads %s > %s",
              qf, db, THREADS, out))
  hits <- read_tsv(out, col_names = c("q", "s", "pident", "len")) |>
    filter(pident >= 80, len >= 50) |> count(q)
  cand$n_self_hits[sub] <- coalesce(hits$n[match(cand$id[sub], hits$q)], 0L)
}
cand <- cand |> mutate(single_copy_nt = n_self_hits == 1)
write_tsv(cand, file.path(WORK, "candidates_genes.tsv"))
cat("\ncandidate genes with >=", MIN_RUN, "bp unique run and single copy:\n")
cand |> filter(single_copy_nt) |> count(strain, route) |>
  pivot_wider(names_from = strain, values_from = n, values_fill = 0) |> print(n = Inf)

# ---- fallback: unique intergenic windows ------------------------------------------------------
ig <- imap_dfr(masks, \(m, cid) {
  genic <- logical(length(m))
  g <- cds |> filter(paste0(strain, "|", seqid) == cid)
  for (i in seq_len(nrow(g))) genic[g$start[i]:g$end[i]] <- TRUE
  genic <- genic[seq_along(m)]      # a gene annotated across the sequence end would extend the vector
  r <- rle(m & !genic)
  ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  keep <- r$values & r$lengths >= MIN_RUN
  tibble(chrom_id = cid, start = starts[keep], end = ends[keep], len = r$lengths[keep])
})
write_tsv(ig, file.path(WORK, "candidates_intergenic.tsv"))
cat("\nunique intergenic windows >=", MIN_RUN, "bp:", nrow(ig), "\n")
