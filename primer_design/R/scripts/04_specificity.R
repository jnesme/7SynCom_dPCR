# Step 04 - in-silico specificity of every designed set against all 7 genomes
#
# bowtie1 lists every binding site (<= 3 mismatches) of every primer and probe on
# all replicons, plasmids included.
#  a) in-silico PCR: pairing the set's own primer sites (also F+F and R+R) must give
#     exactly one product <= 3 kb, at the designed position.
#  b) a lone off-target primer site cannot amplify on its own (that is covered by a),
#     but a near-perfect site competes for primer and can misprime. It is "dangerous"
#     with <= 1 mismatch, or 2 mismatches both outside the 5 3'-terminal bases.
#     Sites with 3 mismatches are the random background (~13 per primer in 36 Mb):
#     only counted (n_3mm_sites) and used for ranking.
source("scripts/helpers.R")

MAX_PROD <- 3000
index <- file.path(WORK, "bt_all")
build_bowtie_index(file.path(WORK, "all_genomes.fna"), index)

d <- read_tsv(file.path(WORK, "designs_raw.tsv")) |> distinct(fwd, rev, .keep_all = TRUE)

oligos <- d |>
  select(set_id, fwd, rev, probe) |>
  pivot_longer(-set_id, names_to = "role", values_to = "seq") |>
  mutate(name = paste0(set_id, "__", role))
sites <- find_sites(set_names(oligos$seq, oligos$name), index, "designs") |>
  separate_wider_delim(name, "__", names = c("set_id", "role"), cols_remove = FALSE) |>
  left_join(d |> select(set_id, strain, seqid, amp_start, amp_end), by = "set_id") |>
  # intended site = inside the designed amplicon (amp_start/amp_end are 1-based)
  mutate(intended = chrom == paste0(strain, "|", seqid) & start >= amp_start - 1 & end <= amp_end)

# ---- a) in-silico PCR per set -------------------------------------------------------------
products <- sites |>
  filter(role != "probe") |>
  group_by(set_id) |>
  group_modify(\(s, key) pcr_products(s, MAX_PROD)) |>
  ungroup() |>
  left_join(d |> select(set_id, strain, seqid, amp_start, amp_end), by = "set_id") |>
  mutate(on_target = chrom == paste0(strain, "|", seqid) & start == amp_start - 1 & end == amp_end)

pcr <- products |>
  group_by(set_id) |>
  summarise(n_products = n(), n_on_target = sum(on_target), n_off_products = sum(!on_target),
            off_products = paste0(chrom[!on_target], ":", start[!on_target] + 1, "-", end[!on_target], collapse = ";"))

# ---- b) off-target binding sites -------------------------------------------------------------
off <- sites |> filter(!intended) |>
  mutate(dangerous = mm <= 1 | (mm == 2 & mm3 == 0))
primer_off <- off |> filter(role != "probe") |>
  group_by(set_id) |>
  summarise(n_dangerous_primer_sites = sum(dangerous), primer_min_offtarget_mm = min(mm),
            n_3mm_sites = sum(mm == 3))
probe_off <- off |> filter(role == "probe") |> group_by(set_id) |> summarise(probe_min_offtarget_mm = min(mm))

report <- d |>
  left_join(pcr, by = "set_id") |>
  left_join(primer_off, by = "set_id") |>
  left_join(probe_off, by = "set_id") |>
  mutate(across(c(n_products, n_on_target, n_off_products, n_dangerous_primer_sites, n_3mm_sites), \(x) replace_na(x, 0L)),
         # no site with <= 3 mismatches found -> at least 4
         across(c(primer_min_offtarget_mm, probe_min_offtarget_mm), \(x) replace_na(x, 4L)),
         off_products = replace_na(off_products, ""),
         pass = n_products == 1 & n_on_target == 1 & n_dangerous_primer_sites == 0)
write_tsv(report, file.path(RES, "specificity_report.tsv"))

report |> group_by(strain) |>
  summarise(designed = n(), on_target_ok = sum(n_on_target == 1), no_off_products = sum(n_off_products == 0),
            no_dangerous_sites = sum(n_dangerous_primer_sites == 0), passed = sum(pass)) |> print()
