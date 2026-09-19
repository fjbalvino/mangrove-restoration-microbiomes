#!/usr/bin/env Rscript
# ============================================================
# 00_02_audit_read_ids.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Read-file and metadata matching; not read processing.
# Inputs (source expressions; complete list in docs/contracts/00_02_audit_read_ids.json):
#   meta_raw <- read_csv(args$meta_csv, show_col_types = FALSE, guess_max = 100000)
#   cohort_raw <- read_csv(args$cohort_csv, show_col_types = FALSE, guess_max = 100000)
# Outputs (source expressions; complete list in contract):
#   writeLines(out_dir, latest_file)
#   write_csv(
#   write_csv(prefix_audit, file.path(out_dir, "tables", "002_metadata_filename_prefix_audit.csv"))
#   write_csv(prefix_summary, file.path(out_dir, "tables", "002_metadata_filename_prefix_summary.csv"))
#   write_csv(r1_files, file.path(out_dir, "tables", "002_raw_R1_files.csv"))
#   write_csv(r2_files, file.path(out_dir, "tables", "002_raw_R2_files.csv"))
#   write_csv(raw_pairs, file.path(out_dir, "tables", "002_raw_read_pairs_audit.csv"))
#   write_csv(raw_pair_summary, file.path(out_dir, "tables", "002_raw_read_pairs_summary.csv"))
#   write_csv(raw_vs_meta, file.path(out_dir, "tables", "002_raw_samples_vs_metadata.csv"))
#   write_csv(raw_not_in_metadata, file.path(out_dir, "tables", "002_raw_samples_not_in_metadata.csv"))
# Algorithmic provenance:
# FASTQ filename/sample-ID inventory and matching to metadata; no read processing.
#   McMurdie & Holmes (2013), phyloseq, doi:10.1371/journal.pone.0061217.
# Source SHA-256: 88aa419d8e3bb8f796676885ce52d2f4deed12412b91728e242d7308c4f85689
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

timestamp_now <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
stop2 <- function(...) stop(sprintf(...), call. = FALSE)
msg <- function(...) message(sprintf(...))

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  defaults <- list(
    meta_csv = "/data/Ciencia-Frontera/Metadata/Metricas_mangrove_metadata_clean.csv",
    cohort_csv = "/home/fjbalvino/Tipping_points/Results/044_build_protein_filtered_matrix_HI_20260529_000843/tables/044_metadata_aligned_51samples.csv",
    expected_n = "51",
    r1_dir = "/data/Ciencia-Frontera/Raw/P240203/Carpeta_Reads/R1",
    r2_dir = "/data/Ciencia-Frontera/Raw/P240203/Carpeta_Reads/R2",
    out_root = "/home/fjbalvino/Tipping_points/resultados_finales"
  )

  if (length(args) == 0) return(defaults)

  if (length(args) %% 2 != 0) {
    stop2("Argumentos inválidos. Usa pares tipo --meta_csv VALUE --r1_dir VALUE.")
  }

  keys <- sub("^--", "", args[seq(1, length(args), by = 2)])
  vals <- args[seq(2, length(args), by = 2)]

  for (i in seq_along(keys)) {
    defaults[[keys[i]]] <- vals[i]
  }

  defaults
}

clean_names <- function(x) {
  x |>
    str_replace_all("\\s+", "_") |>
    str_replace_all("[^A-Za-z0-9_\\.]+", "_") |>
    str_replace_all("_+", "_") |>
    str_replace_all("^_|_$", "") |>
    tolower()
}

first_existing <- function(nms, candidates, label) {
  hit <- intersect(candidates, nms)

  if (length(hit) == 0) {
    stop2(
      "No encontré %s. Candidatas: %s. Columnas disponibles: %s",
      label,
      paste(candidates, collapse = ", "),
      paste(nms, collapse = ", ")
    )
  }

  hit[1]
}

strip_fastq_extensions <- function(x) {
  basename(x) |>
    str_replace("\\.gz$", "") |>
    str_replace("\\.fastq$", "") |>
    str_replace("\\.fq$", "")
}

fallback_extract_sample_id <- function(filename) {
  x <- strip_fastq_extensions(filename)

  # Para nombres tipo:
  # CC1_1_15493AAG_TCCGAATACG-CGGAGAACCT_R1.fastq.gz
  # CZ7_3_15555AAG_XXXX-YYYY_R2.fastq.gz
  # TX2_1_...
  #
  # Extrae solo CC1_1 / CZ7_3 / TX2_1.
  out <- str_match(x, "^([A-Za-z]+[0-9]+_[0-9]+)")[, 2]

  # Fallback más general si no cumple el patrón anterior.
  out2 <- x |>
    str_replace("([._-])L[0-9]{3}([._-])R[12]([._-])001$", "") |>
    str_replace("([._-])L[0-9]{3}([._-])R[12]$", "") |>
    str_replace("([._-])R[12]([._-])001$", "") |>
    str_replace("([._-])R[12]$", "") |>
    str_replace("([._-])[12]([._-])001$", "") |>
    str_replace("([._-])[12]$", "") |>
    str_replace("_154[0-9]+[A-Za-z]+_.*$", "") |>
    str_replace("_[0-9]+[A-Za-z]+_.*$", "")

  ifelse(is.na(out) | !nzchar(out), out2, out)
}

match_metadata_id_from_filename <- function(filename, metadata_ids) {
  x <- strip_fastq_extensions(filename)

  ids <- metadata_ids[!is.na(metadata_ids) & nzchar(metadata_ids)]
  ids <- ids[order(nchar(ids), decreasing = TRUE)]

  hits <- ids[vapply(ids, function(id) {
    if (!startsWith(x, id)) return(FALSE)

    # Asegura frontera: después del ID debe venir fin, _, -, .
    next_pos <- nchar(id) + 1L
    if (nchar(x) < next_pos) return(TRUE)

    next_char <- substr(x, next_pos, next_pos)
    next_char %in% c("_", "-", ".")
  }, logical(1))]

  if (length(hits) > 0) hits[1] else NA_character_
}

infer_sample_code <- function(filename, metadata_ids) {
  matched <- match_metadata_id_from_filename(filename, metadata_ids)

  if (!is.na(matched) && nzchar(matched)) {
    return(matched)
  }

  fallback_extract_sample_id(filename)
}

list_fastq_dir <- function(dir_path, read_label, metadata_ids) {
  if (!dir.exists(dir_path)) {
    stop2("No existe directorio %s: %s", read_label, dir_path)
  }

  files <- list.files(
    dir_path,
    pattern = "\\.(fastq|fq)(\\.gz)?$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )

  if (length(files) == 0) {
    return(tibble(
      read = character(),
      file_path = character(),
      filename = character(),
      sample_code = character(),
      sample_code_source = character(),
      file_size_bytes = numeric()
    ))
  }

  filenames <- basename(files)

  matched_ids <- vapply(
    filenames,
    match_metadata_id_from_filename,
    metadata_ids = metadata_ids,
    FUN.VALUE = character(1)
  )

  fallback_ids <- vapply(
    filenames,
    fallback_extract_sample_id,
    FUN.VALUE = character(1)
  )

  sample_codes <- ifelse(!is.na(matched_ids) & nzchar(matched_ids), matched_ids, fallback_ids)

  tibble(
    read = read_label,
    file_path = normalizePath(files, mustWork = TRUE),
    filename = filenames,
    sample_code = sample_codes,
    sample_code_source = ifelse(!is.na(matched_ids) & nzchar(matched_ids), "metadata_prefix_match", "fallback_regex"),
    matched_metadata_id = matched_ids,
    fallback_sample_code = fallback_ids,
    file_size_bytes = file.info(files)$size
  ) |>
    arrange(sample_code, filename)
}

starts_with_id <- function(filename, id) {
  filename <- as.character(filename)
  id <- as.character(id)

  !is.na(filename) &
    nzchar(filename) &
    !is.na(id) &
    nzchar(id) &
    str_starts(filename, fixed(id))
}

args <- parse_args()

if (!file.exists(args$meta_csv)) {
  stop2("No existe metadata: %s", args$meta_csv)
}
if (!file.exists(args$cohort_csv)) {
  stop2("No existe metadata canonica de cohorte: %s", args$cohort_csv)
}

run_id <- paste0("002_auditar_lecturas_vs_metadata_", timestamp_now())
out_dir <- file.path(args$out_root, run_id)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

latest_file <- file.path(args$out_root, "LATEST_002_auditar_lecturas_vs_metadata.txt")
writeLines(out_dir, latest_file)

log_file <- file.path(out_dir, "logs", "002_log.txt")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

msg("===== 002 auditar lecturas vs metadata =====")
msg("Start: %s", as.character(Sys.time()))
msg("Metadata: %s", args$meta_csv)
msg("Cohorte canonica: %s", args$cohort_csv)
msg("R1 dir:   %s", args$r1_dir)
msg("R2 dir:   %s", args$r2_dir)
msg("Output:   %s", out_dir)

# ------------------------------------------------------------------------------
# Metadata
# ------------------------------------------------------------------------------

meta_raw <- read_csv(args$meta_csv, show_col_types = FALSE, guess_max = 100000)
names(meta_raw) <- clean_names(names(meta_raw))

cohort_raw <- read_csv(args$cohort_csv, show_col_types = FALSE, guess_max = 100000)
names(cohort_raw) <- clean_names(names(cohort_raw))
cohort_id_col <- first_existing(
  names(cohort_raw),
  c("sample_id", "id_metagenomics", "run_accession"),
  "canonical cohort ID column"
)
cohort_ids <- unique(trimws(as.character(cohort_raw[[cohort_id_col]])))
cohort_ids <- cohort_ids[!is.na(cohort_ids) & nzchar(cohort_ids)]
expected_n <- suppressWarnings(as.integer(args$expected_n))
if (!is.finite(expected_n) || expected_n < 1L) stop2("expected_n debe ser entero positivo")
if (expected_n != 51L) stop2("El canon esta congelado en 51 muestras; no se permite redefinirlo")
if (length(cohort_ids) != expected_n || anyDuplicated(cohort_ids)) {
  stop2(
    "La cohorte canonica debe contener exactamente %s IDs unicos; encontro %s",
    expected_n,
    length(cohort_ids)
  )
}

id_col <- first_existing(
  names(meta_raw),
  c("id_metagenomics", "sample_id", "run_accession"),
  "metadata ID column"
)

r1_col <- first_existing(
  names(meta_raw),
  c("r1_filename", "r1", "r1_file", "read1", "forward_reads"),
  "r1_filename column"
)

r2_col <- first_existing(
  names(meta_raw),
  c("r2_filename", "r2", "r2_file", "read2", "reverse_reads"),
  "r2_filename column"
)

meta <- meta_raw |>
  mutate(
    .row_order_original = row_number(),
    id_metagenomics = as.character(.data[[id_col]]),
    r1_filename = as.character(.data[[r1_col]]),
    r2_filename = as.character(.data[[r2_col]])
  )

meta_excluded <- meta |>
  filter(!id_metagenomics %in% cohort_ids)
write_csv(
  meta_excluded,
  file.path(out_dir, "tables", "002_metadata_excluida_fuera_canon_51.csv")
)
meta <- meta |>
  filter(id_metagenomics %in% cohort_ids)

if (n_distinct(meta$id_metagenomics) != expected_n || nrow(meta) != expected_n) {
  stop2(
    "La metadata de lecturas no cubre exactamente el canon de %s muestras",
    expected_n
  )
}

metadata_ids <- unique(meta$id_metagenomics)
metadata_ids <- metadata_ids[!is.na(metadata_ids) & nzchar(metadata_ids)]

msg("Metadata rows retained in canon: %s", nrow(meta))
msg("Unique canonical metadata IDs: %s", length(metadata_ids))
msg("Metadata rows excluded outside canon: %s", nrow(meta_excluded))

# ------------------------------------------------------------------------------
# Check metadata r1/r2 filename prefix
# ------------------------------------------------------------------------------

prefix_audit <- meta |>
  transmute(
    row_order = .row_order_original,
    id_metagenomics,
    r1_filename,
    r2_filename,
    r1_starts_with_id = starts_with_id(r1_filename, id_metagenomics),
    r2_starts_with_id = starts_with_id(r2_filename, id_metagenomics),
    r1_fallback_sample_code = fallback_extract_sample_id(r1_filename),
    r2_fallback_sample_code = fallback_extract_sample_id(r2_filename),
    r1_fallback_equals_id = r1_fallback_sample_code == id_metagenomics,
    r2_fallback_equals_id = r2_fallback_sample_code == id_metagenomics,
    status = case_when(
      r1_starts_with_id & r2_starts_with_id ~ "OK_prefix",
      !r1_starts_with_id & r2_starts_with_id ~ "R1_prefix_mismatch",
      r1_starts_with_id & !r2_starts_with_id ~ "R2_prefix_mismatch",
      TRUE ~ "R1_R2_prefix_mismatch"
    )
  )

prefix_summary <- prefix_audit |>
  count(status, name = "n") |>
  arrange(desc(n))

write_csv(prefix_audit, file.path(out_dir, "tables", "002_metadata_filename_prefix_audit.csv"))
write_csv(prefix_summary, file.path(out_dir, "tables", "002_metadata_filename_prefix_summary.csv"))

msg("Prefix audit summary:")
print(prefix_summary)

# ------------------------------------------------------------------------------
# Raw FASTQ files
# ------------------------------------------------------------------------------

r1_files <- list_fastq_dir(args$r1_dir, "R1", metadata_ids)
r2_files <- list_fastq_dir(args$r2_dir, "R2", metadata_ids)

write_csv(r1_files, file.path(out_dir, "tables", "002_raw_R1_files.csv"))
write_csv(r2_files, file.path(out_dir, "tables", "002_raw_R2_files.csv"))

msg("R1 FASTQ files: %s", nrow(r1_files))
msg("R2 FASTQ files: %s", nrow(r2_files))
msg("R1 unique sample codes: %s", n_distinct(r1_files$sample_code))
msg("R2 unique sample codes: %s", n_distinct(r2_files$sample_code))

sample_code_source_summary <- bind_rows(r1_files, r2_files) |>
  count(read, sample_code_source, name = "n") |>
  arrange(read, sample_code_source)

write_csv(
  sample_code_source_summary,
  file.path(out_dir, "tables", "002_sample_code_source_summary.csv")
)

# ------------------------------------------------------------------------------
# Pair R1/R2 by corrected sample_code
# ------------------------------------------------------------------------------

r1_by_sample <- r1_files |>
  group_by(sample_code) |>
  summarise(
    n_R1_files = n(),
    R1_filenames = paste(filename, collapse = ";"),
    R1_paths = paste(file_path, collapse = ";"),
    R1_total_size_bytes = sum(file_size_bytes, na.rm = TRUE),
    R1_sample_code_sources = paste(unique(sample_code_source), collapse = ";"),
    .groups = "drop"
  )

r2_by_sample <- r2_files |>
  group_by(sample_code) |>
  summarise(
    n_R2_files = n(),
    R2_filenames = paste(filename, collapse = ";"),
    R2_paths = paste(file_path, collapse = ";"),
    R2_total_size_bytes = sum(file_size_bytes, na.rm = TRUE),
    R2_sample_code_sources = paste(unique(sample_code_source), collapse = ";"),
    .groups = "drop"
  )

raw_pairs <- full_join(r1_by_sample, r2_by_sample, by = "sample_code") |>
  mutate(
    n_R1_files = replace_na(n_R1_files, 0L),
    n_R2_files = replace_na(n_R2_files, 0L),
    has_R1 = n_R1_files > 0,
    has_R2 = n_R2_files > 0,
    pair_status = case_when(
      has_R1 & has_R2 & n_R1_files == 1 & n_R2_files == 1 ~ "OK_single_pair",
      has_R1 & has_R2 & (n_R1_files > 1 | n_R2_files > 1) ~ "multiple_files_for_sample",
      has_R1 & !has_R2 ~ "R1_only",
      !has_R1 & has_R2 ~ "R2_only",
      TRUE ~ "unexpected"
    )
  ) |>
  arrange(sample_code)

raw_pair_summary <- raw_pairs |>
  count(pair_status, name = "n") |>
  arrange(desc(n))

write_csv(raw_pairs, file.path(out_dir, "tables", "002_raw_read_pairs_audit.csv"))
write_csv(raw_pair_summary, file.path(out_dir, "tables", "002_raw_read_pairs_summary.csv"))

msg("Raw pair summary:")
print(raw_pair_summary)

# ------------------------------------------------------------------------------
# Compare raw sample codes against metadata IDs
# ------------------------------------------------------------------------------

metadata_ids_tbl <- meta |>
  transmute(
    id_metagenomics,
    in_metadata = TRUE,
    metadata_r1_filename = r1_filename,
    metadata_r2_filename = r2_filename
  ) |>
  distinct(id_metagenomics, .keep_all = TRUE)

raw_vs_meta <- raw_pairs |>
  rename(raw_sample_code = sample_code) |>
  left_join(metadata_ids_tbl, by = c("raw_sample_code" = "id_metagenomics")) |>
  mutate(
    in_metadata = replace_na(in_metadata, FALSE),
    status_metadata = if_else(in_metadata, "present_in_canon_51", "raw_sample_outside_canon_51")
  ) |>
  arrange(status_metadata, raw_sample_code)

raw_not_in_metadata <- raw_vs_meta |>
  filter(!in_metadata) |>
  select(
    raw_sample_code,
    pair_status,
    n_R1_files,
    n_R2_files,
    R1_filenames,
    R2_filenames,
    R1_paths,
    R2_paths,
    R1_total_size_bytes,
    R2_total_size_bytes,
    R1_sample_code_sources,
    R2_sample_code_sources
  ) |>
  arrange(raw_sample_code)

metadata_missing_raw_reads <- metadata_ids_tbl |>
  left_join(raw_pairs, by = c("id_metagenomics" = "sample_code")) |>
  mutate(
    n_R1_files = replace_na(n_R1_files, 0L),
    n_R2_files = replace_na(n_R2_files, 0L),
    has_R1 = n_R1_files > 0,
    has_R2 = n_R2_files > 0,
    raw_status = case_when(
      has_R1 & has_R2 ~ "has_R1_R2",
      has_R1 & !has_R2 ~ "metadata_sample_R1_only",
      !has_R1 & has_R2 ~ "metadata_sample_R2_only",
      TRUE ~ "metadata_sample_missing_raw_reads"
    )
  ) |>
  filter(raw_status != "has_R1_R2") |>
  select(
    id_metagenomics,
    raw_status,
    metadata_r1_filename,
    metadata_r2_filename,
    n_R1_files,
    n_R2_files,
    R1_filenames,
    R2_filenames,
    R1_paths,
    R2_paths
  ) |>
  arrange(id_metagenomics)

write_csv(raw_vs_meta, file.path(out_dir, "tables", "002_raw_samples_vs_metadata.csv"))
write_csv(raw_not_in_metadata, file.path(out_dir, "tables", "002_raw_samples_not_in_metadata.csv"))
write_csv(metadata_missing_raw_reads, file.path(out_dir, "tables", "002_metadata_samples_missing_raw_reads.csv"))

# ------------------------------------------------------------------------------
# Exact declared filename presence
# ------------------------------------------------------------------------------

metadata_declared_file_presence <- meta |>
  transmute(
    id_metagenomics,
    r1_filename,
    r2_filename,
    r1_filename_present_in_R1_dir = r1_filename %in% r1_files$filename,
    r2_filename_present_in_R2_dir = r2_filename %in% r2_files$filename,
    status = case_when(
      r1_filename_present_in_R1_dir & r2_filename_present_in_R2_dir ~ "OK_declared_files_present",
      !r1_filename_present_in_R1_dir & r2_filename_present_in_R2_dir ~ "declared_R1_missing",
      r1_filename_present_in_R1_dir & !r2_filename_present_in_R2_dir ~ "declared_R2_missing",
      TRUE ~ "declared_R1_R2_missing"
    )
  )

declared_file_presence_summary <- metadata_declared_file_presence |>
  count(status, name = "n") |>
  arrange(desc(n))

write_csv(
  metadata_declared_file_presence,
  file.path(out_dir, "tables", "002_metadata_declared_file_presence.csv")
)

write_csv(
  declared_file_presence_summary,
  file.path(out_dir, "tables", "002_metadata_declared_file_presence_summary.csv")
)

msg("Declared file presence summary:")
print(declared_file_presence_summary)

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

run_summary <- tibble(
  key = c(
    "script",
    "timestamp",
    "metadata_csv",
    "canonical_cohort_csv",
    "expected_canonical_samples",
    "r1_dir",
    "r2_dir",
    "out_dir",
    "n_metadata_rows",
    "n_metadata_unique_ids",
    "n_metadata_rows_excluded_outside_canon",
    "n_r1_files",
    "n_r2_files",
    "n_raw_unique_sample_codes",
    "n_raw_samples_not_in_metadata",
    "n_metadata_samples_missing_raw_reads",
    "n_prefix_mismatch_rows",
    "latest_file"
  ),
  value = c(
    "002_auditar_lecturas_vs_metadata.R",
    as.character(Sys.time()),
    args$meta_csv,
    args$cohort_csv,
    as.character(expected_n),
    args$r1_dir,
    args$r2_dir,
    out_dir,
    as.character(nrow(meta)),
    as.character(n_distinct(meta$id_metagenomics)),
    as.character(nrow(meta_excluded)),
    as.character(nrow(r1_files)),
    as.character(nrow(r2_files)),
    as.character(nrow(raw_pairs)),
    as.character(nrow(raw_not_in_metadata)),
    as.character(nrow(metadata_missing_raw_reads)),
    as.character(sum(prefix_audit$status != "OK_prefix")),
    latest_file
  )
)

write_csv(run_summary, file.path(out_dir, "tables", "002_run_summary.csv"))

report <- c(
  "002 auditar lecturas vs metadata compartida",
  "======================================",
  "",
  paste0("Metadata: ", args$meta_csv),
  paste0("R1 dir:   ", args$r1_dir),
  paste0("R2 dir:   ", args$r2_dir),
  paste0("Output:   ", out_dir),
  "",
  "Run summary:",
  paste0(run_summary$key, ": ", run_summary$value),
  "",
  "Prefix audit summary:",
  capture.output(print(prefix_summary)),
  "",
  "Raw pair summary:",
  capture.output(print(raw_pair_summary)),
  "",
  "Declared file presence summary:",
  capture.output(print(declared_file_presence_summary)),
  "",
  "Key output tables:",
  "tables/002_raw_samples_not_in_metadata.csv",
  "tables/002_metadata_samples_missing_raw_reads.csv",
  "tables/002_metadata_filename_prefix_audit.csv",
  "tables/002_raw_read_pairs_audit.csv",
  "tables/002_metadata_declared_file_presence.csv"
)

writeLines(report, file.path(out_dir, "002_audit_report.txt"))

msg("===== DONE 002 =====")
msg("Output: %s", out_dir)
msg("Latest: %s", latest_file)
msg("CSV principal muestras raw no metadata:")
msg("  %s", file.path(out_dir, "tables", "002_raw_samples_not_in_metadata.csv"))
msg("CSV metadata sin raw:")
msg("  %s", file.path(out_dir, "tables", "002_metadata_samples_missing_raw_reads.csv"))
msg("Reporte:")
msg("  %s", file.path(out_dir, "002_audit_report.txt"))
