# Step 00 - prepare per-strain inputs from the NCBI RefSeq files
#   * chromosome / plasmid FASTA per strain + one FASTA with everything
#   * a CDS table parsed from the GFF (one row per gene)
#   * all proteins, renamed "strain|locus_tag"
source("scripts/helpers.R")

strains <- read_tsv(file.path(GENOMES, "assemblies.tsv")) |>
  mutate(dir = map_chr(GenBank_acc, \(a) list.files(file.path(GENOMES, a), pattern = "^refseq_", full.names = TRUE)))

# ---- replicon roles (chromosome vs plasmid) from the assembly report -------------
read_report <- function(dir) {
  f <- list.files(dir, pattern = "_assembly_report.txt$", full.names = TRUE)
  read_tsv(f, comment = "#", col_names = FALSE) |>
    transmute(seqid = X7, role = X4)          # col 7 = RefSeq accession, col 4 = role
}

# ---- genes from the GFF ------------------------------------------------------------
gff_attr <- function(x, key) str_match(x, paste0("(?:^|;)", key, "=([^;]*)"))[, 2]

read_cds <- function(dir, strain) {
  f <- list.files(dir, pattern = "_genomic.gff.gz$", full.names = TRUE)
  read_tsv(f, comment = "#", col_names = c("seqid", "source", "type", "start", "end",
                                           "score", "strand", "phase", "attributes"),
           col_types = "cccddcccc") |>
    filter(type == "CDS") |>
    mutate(strain     = strain,
           locus_tag  = gff_attr(attributes, "locus_tag"),
           protein_id = gff_attr(attributes, "protein_id"),
           gene       = gff_attr(attributes, "gene"),
           product    = URLdecode(replace_na(gff_attr(attributes, "product"), "")),
           pseudo     = replace_na(gff_attr(attributes, "pseudo") == "true", FALSE)) |>
    # a CDS split in several parts shares one locus_tag -> merge into one gene
    group_by(strain, seqid, locus_tag, strand, protein_id, gene, product, pseudo) |>
    summarise(start = min(start), end = max(end)) |>
    ungroup()
}

cds_all <- list(); rep_all <- list()
all_genomes <- DNAStringSet(); all_prot <- AAStringSet()
for (i in seq_len(nrow(strains))) {
  s   <- strains$strain[i]
  dir <- strains$dir[i]
  role <- read_report(dir)

  fna <- list.files(dir, pattern = "_genomic.fna.gz$", full.names = TRUE) |> str_subset("cds_from", negate = TRUE)
  g <- readDNAStringSet(fna)
  names(g) <- word(names(g), 1)
  rep_all[[s]] <- tibble(strain = s, seqid = names(g), length = width(g)) |> left_join(role, by = "seqid")

  names(g) <- paste0(s, "|", names(g))
  writeXStringSet(g, file.path(WORK, paste0(s, ".all.fna")))
  chrom_ids <- paste0(s, "|", role$seqid[role$role == "Chromosome"])
  writeXStringSet(g[names(g) %in% chrom_ids], file.path(WORK, paste0(s, ".chrom.fna")))
  all_genomes <- c(all_genomes, g)

  cds <- read_cds(dir, s) |> left_join(role |> dplyr::rename(replicon = role), by = "seqid")
  cds_all[[s]] <- cds

  # proteins, keyed by locus_tag (identical WP_ ids can be shared by paralogs)
  faa <- readAAStringSet(list.files(dir, pattern = "_protein.faa.gz$", full.names = TRUE))
  names(faa) <- word(names(faa), 1)
  have <- cds |> filter(!is.na(protein_id), protein_id %in% names(faa))
  p <- faa[have$protein_id]
  names(p) <- paste0(s, "|", have$locus_tag)
  all_prot <- c(all_prot, p)
}

writeXStringSet(all_genomes, file.path(WORK, "all_genomes.fna"))
writeXStringSet(all_prot, file.path(WORK, "all_proteins.faa"))
cds <- bind_rows(cds_all) |> mutate(id = paste0(strain, "|", locus_tag))
write_tsv(cds, file.path(WORK, "cds_table.tsv"))
write_tsv(bind_rows(rep_all), file.path(WORK, "replicons.tsv"))
write_tsv(strains |> select(-dir), file.path(WORK, "strains.tsv"))

cds |> count(strain, replicon) |> pivot_wider(names_from = replicon, values_from = n) |> print()
cat("proteins written:", length(all_prot), "\n")
