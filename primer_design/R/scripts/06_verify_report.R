# Step 06 - verify the final pool and write the report files
#  1. in-silico PCR with all 14 primers together: exactly 7 products, one per genome
#  2. knock-out test: remove each target (drop primer sites inside it); only that
#     strain's product must disappear
#  3. restriction enzymes to fragment gDNA before dPCR that cut none of the amplicons
#  4. amplicons.fasta (gBlock standards), oligos.tsv (order sheet), dimer heatmap
source("scripts/helpers.R")

MAX_PROD <- 3000
FLANK    <- 20
ENZYMES  <- c(AluI = "AGCT", BamHI = "GGATCC", BsaI = "GGTCTC", CviQI = "GTAC", EcoRI = "GAATTC",
              HaeIII = "GGCC", HindIII = "AAGCTT", MseI = "TTAA", MspI = "CCGG", PvuII = "CAGCTG",
              Sau3AI = "GATC", XbaI = "TCTAGA")

final <- read_tsv(file.path(RES, "final_multiplex.tsv"))
genomes <- readDNAStringSet(file.path(WORK, "all_genomes.fna"))
index <- file.path(WORK, "bt_all")
build_bowtie_index(file.path(WORK, "all_genomes.fna"), index)

primers <- final |> select(strain, fwd, rev) |>
  pivot_longer(-strain, names_to = "role", values_to = "seq") |>
  mutate(name = paste0(strain, "_", role))
sites <- find_sites(set_names(primers$seq, primers$name), index, "final")
expected <- paste(paste0(final$strain, "|", final$seqid), final$amp_start - 1, final$amp_end)

lines <- character()
say <- function(...) lines <<- c(lines, sprintf(...))

# ---- 1. pooled PCR --------------------------------------------------------------------
prod <- pcr_products(sites, MAX_PROD) |> distinct(chrom, start, end, .keep_all = TRUE)
say("Pooled in-silico PCR (14 primers, all combinations, <=3 mismatches, <= %d bp): %d products", MAX_PROD, nrow(prod))
for (i in seq_len(nrow(prod)))
  say("  %s:%d-%d (%d bp) %s+%s mm %d/%d", prod$chrom[i], prod$start[i] + 1, prod$end[i], prod$len[i],
      prod$left_oligo[i], prod$right_oligo[i], prod$left_mm[i], prod$right_mm[i])
found <- paste(prod$chrom, prod$start, prod$end)
ok1 <- setequal(found, expected)
say("  -> exactly the 7 designed products: %s", if (ok1) "PASS" else "FAIL")

# ---- 2. knock-out test ---------------------------------------------------------------------
ok2 <- TRUE
for (i in seq_len(nrow(final))) {
  r <- final[i, ]
  inside <- sites$chrom == paste0(r$strain, "|", r$seqid) & sites$end > r$amp_start - 1 & sites$start < r$amp_end
  q <- pcr_products(sites[!inside, ], MAX_PROD) |> distinct(chrom, start, end)
  good <- nrow(q) == 6 && !r$strain %in% word(q$chrom, 1, sep = fixed("|"))
  ok2 <- ok2 && good
  say("Knock-out %s: %d products remain, %s product absent: %s", r$strain, nrow(q), r$strain, if (good) "PASS" else "FAIL")
}

# ---- 3. restriction enzymes ----------------------------------------------------------------------
amps <- DNAStringSet(set_names(final$amplicon, final$strain))
say("Restriction enzymes (gDNA fragmentation) vs amplicons / genome cut frequency:")
for (enz in sort(names(ENZYMES))) {
  site <- ENZYMES[[enz]]
  # count on both strands (a non-palindromic site such as BsaI GGTCTC also occurs as GAGACC)
  both <- \(x) vcountPattern(site, x) + if (revcomp(site) != site) vcountPattern(revcomp(site), x) else 0
  cuts_amp <- names(amps)[both(amps) > 0]
  # replicons are circular: append the first bases so a site spanning the origin is counted
  med <- median(both(xscat(genomes, subseq(genomes, 1, nchar(site) - 1))))
  say("  %-8s site %-8s cuts amplicons: %-20s median cuts/replicon %s%s", enz, site,
      if (length(cuts_amp)) paste(cuts_amp, collapse = ",") else "none",
      format(med, big.mark = ","), if (length(cuts_amp) == 0) "   <- usable" else "")
}

# ---- 4. report files ----------------------------------------------------------------------------
chrom_id <- paste0(final$strain, "|", final$seqid)
from <- pmax(1, final$amp_start - FLANK)
to   <- pmin(width(genomes[chrom_id]), final$amp_end + FLANK)
std <- subseq(genomes[chrom_id], from, to)
names(std) <- sprintf("%s_%s %s %s:%d-%d amplicon %d bp + %d bp flanks", final$strain, final$locus_tag,
                      final$species_NCBI, final$seqid, from, to, final$amp_len, FLANK)
writeXStringSet(std, file.path(RES, "amplicons.fasta"), width = 1e4)

oligos <- final |>
  select(strain, locus_tag, species = species_NCBI, well_4plus3, fwd, rev, probe, fwd_tm, rev_tm, probe_tm) |>
  pivot_longer(c(fwd, rev, probe), names_to = "type", values_to = "sequence_5to3") |>
  mutate(tm = round(case_when(type == "fwd" ~ fwd_tm, type == "rev" ~ rev_tm, TRUE ~ probe_tm), 1),
         name = paste0(strain, "_", locus_tag, "_", c(fwd = "F", rev = "R", probe = "P")[type]),
         role = c(fwd = "forward primer", rev = "reverse primer", probe = "hydrolysis probe (5' dye / 3' quencher)")[type],
         length = nchar(sequence_5to3), gc = round(gc_pct(sequence_5to3), 1)) |>
  select(name, strain, species, role, sequence_5to3, length, tm, gc, well_4plus3)
write_tsv(oligos, file.path(RES, "oligos.tsv"))

dimer <- read_tsv(file.path(RES, "dimer_matrix.tsv"))
dl <- dimer |> pivot_longer(-oligo, names_to = "partner", values_to = "dg") |>
  mutate(oligo = factor(oligo, levels = dimer$oligo), partner = factor(partner, levels = dimer$oligo))
p <- ggplot(dl, aes(partner, oligo, fill = dg)) +
  geom_tile() +
  scale_fill_gradient(low = "#0d366b", high = "#cde2fb", limits = c(min(-12, min(dl$dg)), 0),
                      name = "heterodimer dG\n(kcal/mol, 37 C)") +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL, title = "Final pool: pairwise oligo interactions") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid = element_blank())
ggsave(file.path(RES, "dimer_heatmap_R.png"), p, width = 7, height = 6, dpi = 150, bg = "white")

off_diag <- dl |> filter(oligo != partner)
say("Worst oligo-oligo heterodimer dG in pool (excluding self): %.2f kcal/mol", min(off_diag$dg))
writeLines(lines, file.path(RES, "verification.txt"))
cat(lines, sep = "\n")
if (!(ok1 && ok2)) stop("verification FAILED")
