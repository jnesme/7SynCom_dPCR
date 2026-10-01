# Step 05 - choose one set per strain that works in a multiplex
#
# Pool: the best passing set per gene (mid-replichore first, then route class, fewest
# 3-mismatch background sites, primer3 penalty), top 6 genes per strain.
# Between two sets of different strains we compute (ntthal, 37 C, PCR buffer):
#   ANY  = strongest duplex among all 9 oligo pairs (primers and probes)
#   END  = strongest duplex that anchors a primer's 3' end (probes are 3'-blocked)
# All 6^7 combinations are scored: best worst-ANY, then best worst-END, then smallest
# primer Tm spread, then preferred route classes. The best ones are then checked for
# cross-set PCR products and for probes binding another set's amplicon.
source("scripts/helpers.R")

N_PER_STRAIN <- 6
MAX_PROD     <- 3000
PROBE_MIN_MM <- 5          # a probe needs >= 5 mismatches to every other amplicon
index <- file.path(WORK, "bt_all")

report <- read_tsv(file.path(RES, "specificity_report.tsv"))
ok <- report |> filter(pass) |>
  mutate(mid = between(rel_ori, 0.3, 0.7)) |>
  arrange(strain, desc(mid), class_rank, n_3mm_sites, penalty)
missing <- setdiff(unique(report$strain), ok$strain)
if (length(missing) > 0) stop("no specificity-passing set for ", paste(missing, collapse = ", "))

pool <- ok |>
  distinct(strain, locus_tag, .keep_all = TRUE) |>        # best set per gene
  group_by(strain) |> slice_head(n = N_PER_STRAIN) |> ungroup() |>
  mutate(idx = row_number())
strains <- sort(unique(pool$strain))
print(count(pool, strain, mid))

# ---- pairwise interactions between sets ----------------------------------------------------
oligo <- pool |> select(idx, strain, fwd, rev, probe) |>
  pivot_longer(c(fwd, rev, probe), names_to = "role", values_to = "seq")
pairs <- cross_join(oligo, oligo, suffix = c("_a", "_b")) |>
  filter(strain_a != strain_b)

any_pairs <- pairs |> filter(idx_a < idx_b) |>
  mutate(dg = ntthal_dg(seq_a, seq_b, "ANY"))
end_pairs <- pairs |> filter(role_a != "probe") |>             # 3' end of a primer (a) on b
  mutate(dg = ntthal_dg(seq_a, seq_b, "END1"))

n <- nrow(pool)
ANY <- matrix(0, n, n); END <- matrix(0, n, n)      # 0 = neutral (same strain / diagonal)
a <- any_pairs |> group_by(idx_a, idx_b) |> summarise(dg = min(dg))
ANY[cbind(a$idx_a, a$idx_b)] <- a$dg; ANY[cbind(a$idx_b, a$idx_a)] <- a$dg
e <- end_pairs |> mutate(i = pmin(idx_a, idx_b), j = pmax(idx_a, idx_b)) |>
  group_by(i, j) |> summarise(dg = min(dg))
END[cbind(e$i, e$j)] <- e$dg; END[cbind(e$j, e$i)] <- e$dg

# ---- score every 7-way combination (vectorised) -------------------------------------------------
combos <- as.matrix(expand.grid(map(set_names(strains), \(s) pool$idx[pool$strain == s])))
worst <- function(M) {
  w <- rep(0, nrow(combos))
  for (p in combn(ncol(combos), 2, simplify = FALSE))
    w <- pmin(w, M[cbind(combos[, p[1]], combos[, p[2]])])
  w
}
tm_max <- pmax(pool$fwd_tm, pool$rev_tm); tm_min <- pmin(pool$fwd_tm, pool$rev_tm)
scores <- tibble(
  combo = seq_len(nrow(combos)),
  worst_any = worst(ANY), worst_end = worst(END),
  tm_spread = apply(matrix(tm_max[combos], ncol = ncol(combos)), 1, max) -
              apply(matrix(tm_min[combos], ncol = ncol(combos)), 1, min),
  class_sum = rowSums(matrix(pool$class_rank[combos], ncol = ncol(combos)))) |>
  arrange(desc(round(worst_any, 1)), desc(round(worst_end, 1)), round(tm_spread, 1), class_sum)
write_tsv(scores |> select(-combo), file.path(WORK, "combination_scores.tsv"))
cat(sprintf("%s combinations scored; best worst-dG any = %.2f, end = %.2f kcal/mol\n",
            format(nrow(scores), big.mark = ","), scores$worst_any[1], scores$worst_end[1]))

# ---- cross-reactivity check of the best combinations -----------------------------------------------
build_bowtie_index(file.path(WORK, "all_genomes.fna"), index)
primers <- oligo |> filter(role != "probe") |> mutate(name = paste0(idx, "__", role))
psites <- find_sites(set_names(primers$seq, primers$name), index, "pool") |>
  mutate(idx = as.integer(word(name, 1, sep = "__")))

probe_hits_amplicon <- function(probe, amplicon) {
  # does the probe match either strand of the amplicon with < PROBE_MIN_MM mismatches?
  amp <- DNAString(amplicon)
  length(matchPattern(probe, amp, max.mismatch = PROBE_MIN_MM - 1)) > 0 ||
    length(matchPattern(probe, reverseComplement(amp), max.mismatch = PROBE_MIN_MM - 1)) > 0
}

cross_check <- function(ids) {
  sel <- pool |> filter(idx %in% ids)
  prod <- pcr_products(psites |> filter(idx %in% ids), MAX_PROD) |> distinct(chrom, start, end)
  expected <- paste(paste0(sel$strain, "|", sel$seqid), sel$amp_start - 1, sel$amp_end)
  found <- paste(prod$chrom, prod$start, prod$end)
  probe_hits <- 0
  for (i in seq_len(nrow(sel))) for (j in seq_len(nrow(sel)))
    if (i != j && probe_hits_amplicon(sel$probe[i], sel$amplicon[j])) probe_hits <- probe_hits + 1
  list(ok = all(expected %in% found) && all(found %in% expected) && probe_hits == 0,
       n_extra = sum(!found %in% expected), probe_hits = probe_hits)
}

chosen <- NULL
for (k in seq_len(min(500, nrow(scores)))) {
  ids <- combos[scores$combo[k], ]
  chk <- cross_check(ids)
  if (chk$ok) { chosen <- k; break }
  cat(sprintf("combination #%d rejected: %d extra products, %d probe hits\n", k, chk$n_extra, chk$probe_hits))
}
if (is.null(chosen)) stop("no combination passed the cross-reactivity check")
ids <- combos[scores$combo[chosen], ]
cat(sprintf("\nchosen: worst between-set dG %.2f (3' end %.2f) kcal/mol, primer Tm spread %.2f C\n",
            scores$worst_any[chosen], scores$worst_end[chosen], scores$tm_spread[chosen]))

# ---- 4 + 3 split into two wells (for instruments with <= 5 colours) -----------------------------
splits <- combn(ids, 4, simplify = FALSE)
split_score <- map_dbl(splits, \(w1) {
  w2 <- setdiff(ids, w1)
  min(ANY[w1, w1], ANY[w2, w2])
})
w1 <- splits[[which.max(split_score)]]

species <- read_tsv(file.path(WORK, "strains.tsv")) |> select(strain, species_NCBI = NCBI_species)
final <- pool |> filter(idx %in% ids) |>
  mutate(well_4plus3 = if_else(idx %in% w1, "W1", "W2")) |>
  left_join(species, by = "strain") |>
  relocate(species_NCBI, .after = set_id) |>
  arrange(strain)
write_tsv(final |> select(-idx), file.path(RES, "final_multiplex.tsv"))
write_tsv(pool |> select(-idx), file.path(RES, "candidates_pool.tsv"))
write_tsv(ok, file.path(RES, "candidates_all.tsv"))

# oligo x oligo dG matrix of the final pool (diagonal = homodimer)
fo <- final |> select(strain, fwd, rev, probe) |>
  pivot_longer(-strain, names_to = "role", values_to = "seq") |>
  mutate(name = paste0(strain, "_", role))
grid <- expand_grid(a = seq_len(nrow(fo)), b = seq_len(nrow(fo)))
grid$dg <- ntthal_dg(fo$seq[grid$a], fo$seq[grid$b], "ANY")
dimer <- matrix(grid$dg, nrow(fo), byrow = TRUE, dimnames = list(fo$name, fo$name))
write_tsv(as_tibble(round(dimer, 2), rownames = "oligo"), file.path(RES, "dimer_matrix.tsv"))

final |> select(strain, species_NCBI, locus_tag, gene, product, rel_ori, amp_len, fwd, rev, probe, well_4plus3) |>
  print(width = Inf)
