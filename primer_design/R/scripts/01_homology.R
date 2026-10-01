# Step 01 - homology route: which genes have homologs in the other genomes?
#
# All proteins of the 7 genomes are searched against each other (mmseqs2).
#   A_unique         : no homolog in any other genome and no paralog in its own genome
#   B_universal_1to1 : no paralog, and exactly one homolog in each of the 6 other genomes
# Genes that look mobile (transposase, phage, integrase, toxin-antitoxin ...) are flagged,
# together with every gene within 5 kb of them.
source("scripts/helpers.R")

EVALUE  <- 1e-5
MIN_COV <- 0.3       # alignment length / length of the shorter protein
FLANK   <- 5000

m8 <- file.path(WORK, "allvsall.m8")
if (!file.exists(m8)) {
  prot <- file.path(WORK, "all_proteins.faa")
  run(sprintf("mmseqs easy-search %s %s %s %s/mmseqs_tmp -s 7.5 -e %g --max-seqs 1000 --threads %s -v 1 \\
              --format-output query,target,pident,alnlen,evalue,bits,qlen,tlen",
              prot, prot, m8, WORK, EVALUE, THREADS))
}

hits <- read_tsv(m8, col_names = c("q", "t", "pident", "alnlen", "evalue", "bits", "qlen", "tlen")) |>
  filter(q != t, evalue <= EVALUE, alnlen / pmin(qlen, tlen) >= MIN_COV) |>
  mutate(q_strain = word(q, 1, sep = fixed("|")), t_strain = word(t, 1, sep = fixed("|")))

cds <- read_tsv(file.path(WORK, "cds_table.tsv"))
all_strains <- sort(unique(cds$strain))

# number of distinct hit proteins per query and per target genome
per_genome <- hits |>
  distinct(q, t, q_strain, t_strain) |>
  count(q, q_strain, t_strain, name = "n_hits")

summary_hits <- per_genome |>
  group_by(q) |>
  summarise(n_paralogs             = sum(n_hits[t_strain == q_strain]),
            n_genomes_with_homolog = sum(t_strain != q_strain),
            max_homologs_other     = max(c(0, n_hits[t_strain != q_strain])))

best_other <- hits |> filter(q_strain != t_strain) |> group_by(q) |> summarise(best_other_pident = max(pident))

cds <- cds |>
  left_join(summary_hits, by = c("id" = "q")) |>
  left_join(best_other, by = c("id" = "q")) |>
  mutate(across(c(n_paralogs, n_genomes_with_homolog, max_homologs_other, best_other_pident), \(x) replace_na(x, 0)),
         class = case_when(
           n_paralogs == 0 & n_genomes_with_homolog == 0 ~ "A_unique",
           n_paralogs == 0 & n_genomes_with_homolog == length(all_strains) - 1 & max_homologs_other == 1 ~ "B_universal_1to1",
           TRUE ~ ""))

# ---- mobile elements -----------------------------------------------------------------
mobile_re <- regex(paste("transpos", "integrase", "insertion element", "insertion sequence", "phage",
                         "prophage", "site-specific recombinase", "resolvase", "invertase", "conjugal",
                         "conjugative", "relaxase", "mobiliz", "toxin", "antitoxin",
                         "reverse transcriptase", "plasmid", sep = "|"), ignore_case = TRUE)
is_re <- "\\bIS\\d|\\bIS[A-Z]"            # IS families, case-sensitive (so not "isomerase")
cds <- cds |>
  mutate(mobile = (str_detect(product, mobile_re) | str_detect(product, is_re)) &
                  !str_detect(product, "Holliday"))   # RuvC/RuvX are core genes

# genes within FLANK bp of a mobile gene (same replicon)
mob <- cds |> filter(mobile) |> select(strain, seqid, m_start = start, m_end = end)
near <- cds |>
  select(id, strain, seqid, start, end) |>
  inner_join(mob, by = c("strain", "seqid"), relationship = "many-to-many") |>
  filter(m_start - FLANK <= end, m_end + FLANK >= start) |>
  distinct(id)
cds <- cds |> mutate(near_mobile = id %in% near$id)

write_tsv(cds, file.path(WORK, "cds_classified.tsv"))
cds |> count(strain, class = if_else(class == "", "(none)", class)) |> pivot_wider(names_from = class, values_from = n) |> print()
cat("\neligible (chromosomal, not near mobile, not pseudo, >= 300 bp):\n")
cds |> filter(class != "", replicon == "Chromosome", !pseudo, !near_mobile, end - start + 1 >= 300) |>
  count(strain, class) |> pivot_wider(names_from = class, values_from = n) |> print()
