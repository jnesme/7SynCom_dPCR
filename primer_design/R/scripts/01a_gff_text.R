# Step 01a - text route: find unique genes by reading the GFF annotation only.
#
# Each gene gets a text key: its gene symbol (gene=) if present, otherwise its product.
# Uninformative products (hypothetical, "... family protein", DUF ...) get no key.
#   text_unique  : the key occurs once in this genome and in no other genome
#   text_core_1x : the key occurs exactly once in every genome (e.g. gyrB, recA, dnaK)
# No sequence is compared here; step 02 checks the DNA.
source("scripts/helpers.R")

generic <- regex("hypothetical|family|domain-containing|DUF\\d|uncharacteri[sz]ed|putative|unknown|^protein$|YbaB|-like",
                 ignore_case = TRUE)

cds <- read_tsv(file.path(WORK, "cds_table.tsv")) |>
  mutate(product = replace_na(product, ""),
         # PGAP adds _1, _2 to duplicated gene symbols -> strip it
         text_key = case_when(
           !is.na(gene) & gene != ""                 ~ paste0("gene:", str_remove(gene, "_\\d+$")),
           product == "" | str_detect(product, generic) ~ "",
           TRUE                                      ~ paste0("product:", str_to_lower(str_trim(product)))),
         # the product string is also checked, in case the same gene has no symbol elsewhere
         prod_key = if_else(product != "", paste0("product:", str_to_lower(str_trim(product))), ""))
all_strains <- sort(unique(cds$strain))
n_strains <- length(all_strains)

# how many times does each key occur in each genome?
key_counts <- cds |> filter(text_key != "") |> count(text_key, strain) |>
  complete(text_key, strain = all_strains, fill = list(n = 0))
prod_counts <- cds |> filter(prod_key != "") |> count(prod_key, strain) |>
  complete(prod_key, strain = all_strains, fill = list(n = 0))

key_summary <- key_counts |> group_by(text_key) |>
  summarise(n_genomes = sum(n > 0), all_once = all(n == 1), total = sum(n))
prod_summary <- prod_counts |> group_by(prod_key) |> summarise(prod_genomes = sum(n > 0))

cds <- cds |>
  left_join(key_summary, by = "text_key") |>
  left_join(prod_summary, by = "prod_key") |>
  mutate(prod_comparable = prod_key != "" & !str_detect(product, generic),
         text_class = case_when(
           text_key == "" ~ "",
           # once, in one genome only (and the product string is not found elsewhere either)
           total == 1 & n_genomes == 1 & (!prod_comparable | prod_genomes == 1) ~ "text_unique",
           all_once & n_genomes == n_strains ~ "text_core_1x",
           TRUE ~ ""))

write_tsv(cds |> select(id, text_key, text_class), file.path(WORK, "cds_text.tsv"))
cds |> count(strain, text_class = if_else(text_class == "", "(none)", text_class)) |> pivot_wider(names_from = text_class, values_from = n) |> print()
cat("\nexamples of text_core_1x keys:",
    paste(head(sort(unique(cds$text_key[cds$text_class == "text_core_1x"])), 25), collapse = ", "), "\n")
