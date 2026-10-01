# Step 07 - figures (PNG + PDF in figures/)
#   Fig 1  workflow, with the number of genes / assays left after each step
#   Fig 2  circular maps of the primary chromosomes: unique DNA, passing assays, chosen target
#   Fig 3  in-silico specificity: assay x genome matrix + knock-out test
#   Fig 4  multiplex compatibility: oligo dG matrix + score of all candidate 7-plexes
# Colours: one blue ramp for magnitudes, orange only for the chosen assays, greys otherwise.
source("scripts/helpers.R")
suppressPackageStartupMessages({ library(circlize); library(patchwork) })

INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#8a8983"; GRID <- "#e4e3df"; PANEL <- "#f5f4f1"
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; GREY <- "#b9b8b3"
RAMP <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
theme_set(theme_minimal(base_size = 8) +
            theme(text = element_text(colour = INK), axis.text = element_text(colour = INK2),
                  panel.grid = element_blank(), plot.title = element_text(size = 9, face = "bold"),
                  plot.title.position = "plot"))
save_fig <- function(p, name, w, h) {
  for (ext in c("png", "pdf")) ggsave(file.path(FIG, paste0(name, "_R.", ext)), p, width = w, height = h,
                                      dpi = 300, bg = "white", device = if (ext == "pdf") cairo_pdf else NULL)
  cat("wrote", paste0(name, "_R"), "\n")
}

# ---- data ----------------------------------------------------------------------------------
strains <- read_tsv(file.path(WORK, "strains.tsv"))
reps    <- read_tsv(file.path(WORK, "replicons.tsv"))
cds     <- read_tsv(file.path(WORK, "cds_classified.tsv")) |>
  left_join(read_tsv(file.path(WORK, "cds_text.tsv")), by = "id") |>
  mutate(class = replace_na(class, ""), text_class = replace_na(text_class, ""))
cand    <- read_tsv(file.path(WORK, "candidates_genes.tsv"))
designs <- read_tsv(file.path(WORK, "designs_raw.tsv"))
spec    <- read_tsv(file.path(RES, "specificity_report.tsv"))
pool    <- read_tsv(file.path(RES, "candidates_pool.tsv"))
final   <- read_tsv(file.path(RES, "final_multiplex.tsv")) |> arrange(strain)
scores  <- read_tsv(file.path(WORK, "combination_scores.tsv"))
masks   <- readRDS(file.path(WORK, "unique_masks.rds"))
dnaA    <- cds |> filter(gene == "dnaA") |> distinct(strain, .keep_all = TRUE) |>
  transmute(strain, chrom_id = paste0(strain, "|", seqid), ori = start)
S <- final$strain
pct_unique <- map_dbl(set_names(dnaA$chrom_id, dnaA$strain), \(id) 100 * mean(masks[[id]]))

# =============================================================================================
# Figure 1 - workflow
# =============================================================================================
eligible <- cds |> filter((class != "" | text_class != ""), replicon == "Chromosome", !pseudo, !near_mobile,
                          end - start + 1 >= 300)
cand_primary <- cand |> filter(single_copy_nt, chrom_id %in% dnaA$chrom_id)
n_combos <- nrow(scores)
fmt <- \(x) format(x, big.mark = ",")

steps <- tribble(
  ~y,   ~title,                                                            ~detail,
  10,   "7 complete genomes (PRJNA357031), RefSeq annotation",
        sprintf("%s CDS  ·  %d chromosomes, %d plasmids", fmt(nrow(cds)), sum(reps$role == "Chromosome"), sum(reps$role == "Plasmid")),
  6.2,  "Eligible candidate genes",
        sprintf("%s genes  ·  route A or B, chromosomal, not pseudo, > 5 kb from mobile elements, >= 300 bp", fmt(nrow(eligible))),
  5.0,  "Sequence-aware filter: unique DNA + single copy",
        sprintf("%s genes on the dnaA chromosome  ·  >= 150 bp of unique 18-mers, one BLAST hit in own genome", fmt(nrow(cand_primary))),
  3.8,  "Primer / probe design (primer3)",
        sprintf("%s assays on %d genes  ·  40 core + 40 unique genes per strain, 70-150 bp amplicons", fmt(nrow(designs)), n_distinct(designs$locus_tag)),
  2.6,  "Exhaustive in-silico PCR (bowtie1, <= 3 mismatches, 7 genomes + plasmids)",
        sprintf("%s assays pass  ·  one on-target product, no near-perfect off-target primer site", fmt(sum(spec$pass))),
  1.4,  "Multiplex selection",
        sprintf("%d assays (top 6 genes per strain)  ->  %s 7-plex combinations scored on dimer dG", nrow(pool), fmt(n_combos)),
  0.2,  "Final 7-plex, verified in silico",
        sprintf("pooled PCR gives exactly 7 products  ·  knock-out test 7/7  ·  worst between-assay dG %.2f kcal/mol", scores$worst_any[1])
) |> mutate(final = row_number() == n())

routes <- tibble(
  x = c(0, 3.4, 6.8),
  title = c("A · GFF text parsing", "B · Protein homology (mmseqs2)", "C · 18-mer uniqueness"),
  l1 = c(sprintf("unique name/product: %s", fmt(sum(cds$text_class == "text_unique"))),
         sprintf("no homolog in other 6: %s", fmt(sum(cds$class == "A_unique"))),
         "unique fraction of the"),
  l2 = c(sprintf("named once in all 7: %s", fmt(sum(cds$text_class == "text_core_1x"))),
         sprintf("1:1 in all 7 genomes: %s", fmt(sum(cds$class == "B_universal_1to1"))),
         sprintf("primary chromosome: %.0f-%.0f %%", min(pct_unique), max(pct_unique))))
BW <- 10; BH <- 0.8; RW <- 3.2; RY <- 7.5; RH <- 1.5

p1 <- ggplot() +
  # route boxes
  geom_rect(data = routes, aes(xmin = x, xmax = x + RW, ymin = RY, ymax = RY + RH), fill = PANEL, colour = GRID) +
  geom_text(data = routes, aes(x = x + 0.12, y = RY + RH - 0.3, label = title), hjust = 0, size = 2.7, fontface = "bold") +
  geom_text(data = routes, aes(x = x + 0.12, y = RY + RH - 0.75, label = l1), hjust = 0, size = 2.4, colour = INK2) +
  geom_text(data = routes, aes(x = x + 0.12, y = RY + RH - 1.1, label = l2), hjust = 0, size = 2.4, colour = INK2) +
  # step boxes
  geom_rect(data = steps, aes(xmin = 0, xmax = BW, ymin = y, ymax = y + BH, colour = final, linewidth = final), fill = "white") +
  scale_colour_manual(values = c(`FALSE` = MUTED, `TRUE` = ORANGE), guide = "none") +
  scale_linewidth_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.8), guide = "none") +
  geom_text(data = steps, aes(x = 0.15, y = y + 0.56, label = title), hjust = 0, size = 2.8, fontface = "bold") +
  geom_text(data = steps, aes(x = 0.15, y = y + 0.22, label = detail), hjust = 0, size = 2.35, colour = INK2) +
  # arrows: genomes -> routes, routes A/B -> eligible, route C -> sequence filter, then down the chain
  annotate("segment", x = routes$x + RW / 2, xend = routes$x + RW / 2, y = 10, yend = RY + RH + 0.05,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"), colour = MUTED, linewidth = 0.3) +
  annotate("segment", x = routes$x[1:2] + RW / 2, xend = routes$x[1:2] + RW / 2, y = RY, yend = 6.2 + BH + 0.05,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"), colour = MUTED, linewidth = 0.3) +
  annotate("segment", x = routes$x[3] + RW / 2, xend = routes$x[3] + RW / 2, y = RY, yend = 5.0 + BH + 0.05,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"), colour = MUTED, linewidth = 0.3, linetype = "22") +
  annotate("segment", x = 1.7, xend = 1.7, y = head(steps$y[-1], -1), yend = steps$y[-(1:2)] + BH + 0.05,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"), colour = MUTED, linewidth = 0.3) +
  coord_cartesian(xlim = c(-0.05, BW + 0.05), ylim = c(0.05, 10.9), expand = FALSE) +
  theme_void()
save_fig(p1, "fig1_workflow", 7.2, 6.4)

# =============================================================================================
# Figure 2 - chromosome maps (circlize, base graphics)
# =============================================================================================
BIN <- 10000
passing <- spec |> filter(pass)
seq_col <- colorRamp2(seq(0, 1, length.out = length(RAMP)), RAMP)

draw_chromosome <- function(s) {
  o <- dnaA |> filter(strain == s)
  m <- masks[[o$chrom_id]]; L <- length(m)
  rot <- \(pos) (pos - o$ori) %% L                        # dnaA at position 0 = top
  f <- final |> filter(strain == s)
  circos.clear()
  circos.par(start.degree = 90, clock.wise = TRUE, gap.after = 0, cell.padding = c(0, 0, 0, 0),
             track.margin = c(0.004, 0.004), points.overflow.warning = FALSE)
  circos.initialize("chr", xlim = c(0, L))
  # track 1: fraction of unique positions per 10 kb bin
  nb <- L %/% BIN
  frac <- colMeans(matrix(m[seq_len(nb * BIN)], nrow = BIN))
  circos.track(ylim = c(0, 1), track.height = 0.12, bg.border = NA, panel.fun = function(x, y) {
    st <- rot((seq_len(nb) - 1) * BIN + 1)
    circos.rect(st, 0, pmin(st + BIN, L), 1, col = seq_col(frac), border = NA)
  })
  # track 2: mid-replichore band + passing assays + chosen target
  circos.track(ylim = c(0, 1), track.height = 0.22, bg.border = NA, panel.fun = function(x, y) {
    for (b in list(c(0.15, 0.35), c(0.65, 0.85))) circos.rect(b[1] * L, 0, b[2] * L, 1, col = "#f2f1ee", border = NA)
    pa <- passing |> filter(strain == s, paste0(strain, "|", seqid) == o$chrom_id)
    if (nrow(pa)) circos.segments(rot(pa$amp_start), 0.25, rot(pa$amp_start), 0.9, col = GREY, lwd = 0.5)
    circos.segments(rot(f$amp_start), 0, rot(f$amp_start), 1, col = ORANGE, lwd = 2.2)
    circos.points(rot(f$amp_start), 1, pch = 21, bg = ORANGE, col = "white", cex = 0.9)
  })
  sp <- strains$NCBI_species[strains$strain == s]
  gene <- if (!is.na(f$gene) && f$gene != "") f$gene else str_remove(f$locus_tag, "^.*_")
  text(0, 1.06, "ori", cex = 0.5, col = INK)            # upright labels outside the ring
  text(0, -1.06, "ter", cex = 0.5, col = INK2)
  text(0, 0.16, s, font = 2, cex = 0.9)
  text(0, 0.0, bquote(italic(.(paste0(substr(sp, 1, 1), ". ", word(sp, 2))))), cex = 0.62, col = INK2)
  text(0, -0.15, sprintf("%.2f Mb · %.0f%% unique", L / 1e6, pct_unique[[s]]), cex = 0.5, col = MUTED)
  mtext(sprintf("target: %s · rel_ori %.2f", gene, f$rel_ori), side = 1, line = -0.2, cex = 0.5, col = ORANGE)
  circos.clear()
}

draw_legend <- function() {
  plot.new(); plot.window(c(0, 1), c(0, 1))
  text(0.02, 0.97, "Outer ring: unique fraction\n(18-mers, 10 kb bins)", adj = c(0, 1), cex = 0.62)
  xs <- seq(0.05, 0.85, length.out = 60)
  rect(xs, 0.72, xs + diff(xs)[1], 0.78, col = seq_col(seq(0, 1, length.out = 60)), border = NA)
  text(c(0.05, 0.87), 0.68, c("0", "1"), cex = 0.55, col = INK2)
  segments(0.06, 0.52, 0.06, 0.62, col = GREY, lwd = 1); text(0.13, 0.57, "assay passing specificity", adj = 0, cex = 0.6)
  segments(0.06, 0.37, 0.06, 0.47, col = ORANGE, lwd = 2.2); points(0.06, 0.47, pch = 21, bg = ORANGE, col = "white")
  text(0.13, 0.42, "chosen target", adj = 0, cex = 0.6)
  rect(0.02, 0.2, 0.1, 0.3, col = "#f2f1ee", border = NA); text(0.13, 0.25, "mid-replichore (rel_ori 0.3-0.7)", adj = 0, cex = 0.6)
  text(0.02, 0.08, "dnaA (origin) at the top, clockwise", adj = 0, cex = 0.6)
}

for (ext in c("png", "pdf")) {
  f <- file.path(FIG, paste0("fig2_chromosome_maps_R.", ext))
  if (ext == "png") png(f, width = 7.2, height = 4.2, units = "in", res = 300) else cairo_pdf(f, width = 7.2, height = 4.2)
  par(mfrow = c(2, 4), mar = c(0.8, 0.2, 0.2, 0.2), bg = "white")
  for (s in S) draw_chromosome(s)
  draw_legend()
  dev.off()
}
cat("wrote fig2_chromosome_maps_R\n")

# =============================================================================================
# Figure 3 - specificity
# =============================================================================================
index <- file.path(WORK, "bt_all")
primers <- final |> select(strain, fwd, rev) |> pivot_longer(-strain, names_to = "role", values_to = "seq") |>
  mutate(name = paste0(strain, "_", role))
sites <- find_sites(set_names(primers$seq, primers$name), index, "fig3") |>
  mutate(assay = word(name, 1, sep = "_"), genome = word(chrom, 1, sep = fixed("|")))

best_mm <- sites |> group_by(assay, genome) |> summarise(mm = min(mm)) |>
  right_join(expand_grid(assay = S, genome = S), by = c("assay", "genome")) |>
  left_join(final |> select(assay = strain, amp_len), by = "assay") |>
  mutate(on = assay == genome,
         label = case_when(on ~ paste0(amp_len, "\nbp"), is.na(mm) ~ ">=4", TRUE ~ as.character(mm)))

pA <- ggplot(best_mm, aes(genome, assay)) +
  geom_tile(aes(fill = on), width = 0.92, height = 0.92) +
  geom_text(aes(label = label, colour = on, fontface = if_else(on, "bold", "plain")), size = 2.4) +
  scale_fill_manual(values = c(`TRUE` = BLUE, `FALSE` = PANEL), guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK2), guide = "none") +
  scale_y_discrete(limits = rev, labels = \(x) paste(x, "assay")) +
  coord_equal() +
  labs(x = "genome (all replicons)", y = NULL, title = "A  Products and closest off-target primer site")

expected_key <- paste(paste0(final$strain, "|", final$seqid), final$amp_start - 1, final$amp_end)
ko <- map_dfr(S, \(masked) {
  r <- final |> filter(strain == masked)
  inside <- sites$chrom == paste0(r$strain, "|", r$seqid) & sites$end > r$amp_start - 1 & sites$start < r$amp_end
  prod <- pcr_products(sites[!inside, ], 3000)
  found <- paste(prod$chrom, prod$start, prod$end)
  tibble(masked = masked, assay = S, formed = expected_key %in% found)
})
pB <- ggplot(ko, aes(assay, masked)) +
  geom_tile(data = filter(ko, formed), fill = BLUE, width = 0.92, height = 0.92) +
  geom_tile(data = filter(ko, !formed), fill = "white", colour = ORANGE, linewidth = 0.6, width = 0.9, height = 0.9) +
  geom_text(data = filter(ko, !formed), label = "lost", colour = ORANGE, fontface = "bold", size = 2.3) +
  scale_y_discrete(limits = rev, labels = \(x) paste(x, "target masked")) +
  coord_equal() +
  labs(x = "assay", y = NULL, title = "B  Knock-out test (all 14 primers pooled)")

p3 <- (pA | pB) + plot_annotation(
  caption = paste0("A: blue = in-silico PCR product of the designed length; grey cells give the fewest mismatches of any primer site\n",
                   "in that genome (bowtie1, exhaustive; >=4 = no site with <=3). No off-target product with <=3 mismatches per primer.\n",
                   "B: each row masks one target amplicon; blue = product still formed; only the masked target is lost."),
  theme = theme(plot.caption = element_text(hjust = 0, size = 6.5, colour = INK2)))
save_fig(p3, "fig3_specificity", 7.2, 3.8)

# =============================================================================================
# Figure 4 - multiplex compatibility
# =============================================================================================
dimer <- read_tsv(file.path(RES, "dimer_matrix.tsv"))
well <- set_names(final$well_4plus3, final$strain)
ord <- final |> arrange(well_4plus3, strain) |> pull(strain)
lvl <- as.vector(t(outer(ord, c("fwd", "rev", "probe"), paste, sep = "_")))
dl <- dimer |> pivot_longer(-oligo, names_to = "partner", values_to = "dg") |>
  mutate(oligo = factor(oligo, lvl), partner = factor(partner, lvl),
         between = word(oligo, 1, sep = "_") != word(partner, 1, sep = "_"))
worst_pair <- dl |> filter(between) |> slice_min(dg, n = 1, with_ties = FALSE)
nice <- \(x) str_replace_all(x, c("_fwd" = " F", "_rev" = " R", "_probe" = " P"))
w1_n <- 3 * sum(final$well_4plus3 == "W1")

p4a <- ggplot(dl, aes(partner, oligo, fill = pmax(dg, -10))) +
  geom_tile() +
  geom_tile(data = worst_pair, fill = NA, colour = ORANGE, linewidth = 0.8) +
  geom_vline(xintercept = w1_n + 0.5, colour = INK, linewidth = 0.6) +
  geom_hline(yintercept = length(lvl) - w1_n + 0.5, colour = INK, linewidth = 0.6) +
  scale_fill_gradientn(colours = rev(RAMP), limits = c(-10, 0), breaks = c(-10, -8, -6, -4, -2, 0),
                       labels = c("<= -10", "-8", "-6", "-4", "-2", "0"), name = "heterodimer dG\n(kcal/mol, 37 C)") +
  scale_x_discrete(labels = nice) + scale_y_discrete(limits = rev, labels = nice) +
  coord_equal() +
  labs(x = NULL, y = NULL, title = "A  Oligo interactions in the final pool",
       subtitle = sprintf("black lines: wells W1 | W2; diagonal blocks: oligos of the same assay\norange: strongest between-assay pair (%s x %s, %.2f kcal/mol)",
                          nice(worst_pair$oligo), nice(worst_pair$partner), worst_pair$dg)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
        axis.text.y = element_text(size = 5.5), plot.subtitle = element_text(size = 6.5, colour = INK2),
        legend.title = element_text(size = 6.5), legend.text = element_text(size = 6))

chosen_val <- scores$worst_any[1]
p4b <- ggplot(scores, aes(worst_any)) +
  geom_histogram(binwidth = 0.25, boundary = 0, fill = GREY, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = chosen_val, colour = ORANGE, linewidth = 0.8) +
  annotate("text", x = chosen_val, y = Inf, label = sprintf("chosen\n%.2f ", chosen_val), hjust = 1.1, vjust = 1.2,
           colour = ORANGE, size = 2.4) +
  scale_y_continuous(labels = \(v) if_else(v == 0, "0", paste0(v / 1000, "k")), expand = expansion(c(0, 0.05))) +
  scale_x_continuous(expand = expansion(add = c(0.3, 0.4))) +
  labs(x = "worst between-assay dG (kcal/mol)", y = "7-plex combinations",
       # the add-on set (results folder with candidates_pool_wide.tsv) was chosen from a much
       # wider pool than the one enumerated here, so the title must not claim "all" candidates
       title = if (file.exists(file.path(RES, "candidates_pool_wide.tsv")))
         sprintf("B  Reduced pool only\n    %s 7-plexes\n    (6 genes per strain)", fmt(nrow(scores)))
       else sprintf("B  All %s candidate 7-plexes", fmt(nrow(scores)))) +
  theme(panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3), axis.line.x = element_line(colour = MUTED))

p4 <- p4a + p4b + plot_layout(widths = c(1.35, 1))
save_fig(p4, "fig4_multiplex_compatibility", 7.2, 4.0)
