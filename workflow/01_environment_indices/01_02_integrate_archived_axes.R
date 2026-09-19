#!/usr/bin/env Rscript
# ============================================================
# 01_02_integrate_archived_axes.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Historical integration step; its MHI output is superseded by 003b. Never substitute a newly fitted HI.
# Inputs (source expressions; complete list in docs/contracts/01_02_integrate_archived_axes.json):
#   trimws(readLines(f, warn = FALSE, n = 1L))
#   x <- read_csv(path, show_col_types = FALSE, guess_max = 100000)
# Outputs (source expressions; complete list in contract):
#   writeLines(
#   write_csv(
#   write_csv(profile_check, file.path(out_dir, "tables", "003_auditoria_perfiles_canon_51.csv"))
#   write_csv(updated_out, out_csv)
#   writeLines(out_dir, latest_file)
# Algorithmic provenance:
# Join archived environmental axes and index scores by sample ID; historical integration.
#   Study-specific environmental calibration and data integration; see docs/METHODS.md.
# Source SHA-256: e91fe93ce7113a7de81d8275f2836ac9c30ac003283fa3f5bab022a2193a8dee
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


# ==============================================================================
# 003_integrar_ejes_ECI_HI_en_metadata.R
#
# Objetivo:
#   Construir la metadata integrada del canon unico de 51 muestras.
#
#   Incorporando:
#     - ECI_PC1, ECI_PC2, ECI_PC3
#     - HI / RBMA1
#     - vegetation_landscape_PC1
#     - water_inundation_PC1
#     - moisture_stress_PC1
#     - physicochemical_PC1
#     - nutrients_redox_PC1
#
#   Manteniendo el orden estricto de las 51 filas canonicas.
#
#   Matching canónico:
#     base$sample_id = aux$sample_id (con aliases de ID auditados)
#
#   Además:
#     - Recodea stage: Preserved -> Conserved
#     - Elimina aliases ECI heredados env_global_PC1-3 del producto integrado
#     - Crea backup antes de sobrescribir
#     - Exporta auditorías de matching y columnas
#     - Exige evidencia de que ECI_PC1-3 fueron recalculados por el script 101
#       después de congelar la cohorte de 51 muestras.
#
# Default: escribe un archivo nuevo dentro de resultados_finales.
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

timestamp_now <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

message2 <- function(...) {
  message(sprintf(...))
}

stop2 <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  results_root <- "/home/fjbalvino/Tipping_points/resultados_finales"
  latest_output <- function(id) {
    f <- file.path(results_root, paste0("LATEST_", id, ".txt"))
    if (!file.exists(f)) return(NA_character_)
    trimws(readLines(f, warn = FALSE, n = 1L))
  }
  run101 <- latest_output("101_construir_ejes_ambientales_por_bloque")
  run203 <- latest_output("203_construir_HI_residualizado")

  defaults <- list(
    base_csv = "/home/fjbalvino/Tipping_points/Results/044_build_protein_filtered_matrix_HI_20260529_000843/tables/044_metadata_aligned_51samples.csv",
    eci_csv = if (!is.na(run101)) file.path(run101, "tables", "101_metadata_with_env_axes.csv") else NA_character_,
    env_axes_csv = if (!is.na(run101)) file.path(run101, "tables", "101_metadata_with_env_axes.csv") else NA_character_,
    rbma_csv = if (!is.na(run203)) file.path(run203, "tables", "203_HI_por_muestra.csv") else NA_character_,
    mhi_csv = "/home/fjbalvino/R/Ciencia-Frontera/resultados/47_mhi_gam_pco1_singletrend_20260216_223915/mhi_local_global_by_sample.csv",
    out_csv = NA_character_,
    out_root = results_root,
    overwrite = "0",
    stop_if_duplicate_ids = "1",
    expected_n = "51",
    expected_profiles = "17"
  )

  if (length(args) == 0) return(defaults)

  if (length(args) %% 2 != 0) {
    stop2("Argumentos inválidos. Usa pares tipo --base_csv VALUE --overwrite 1.")
  }

  parsed <- defaults
  keys <- args[seq(1, length(args), by = 2)]
  vals <- args[seq(2, length(args), by = 2)]
  keys <- sub("^--", "", keys)

  for (i in seq_along(keys)) {
    parsed[[keys[i]]] <- vals[i]
  }

  parsed
}

clean_names <- function(x) {
  x %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("[^A-Za-z0-9_\\.]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    tolower()
}

first_existing <- function(nms, candidates, required = TRUE, label = "column") {
  hit <- intersect(candidates, nms)

  if (length(hit) == 0) {
    if (required) {
      stop2(
        "No encontré %s. Candidatas: %s. Columnas disponibles: %s",
        label,
        paste(candidates, collapse = ", "),
        paste(nms, collapse = ", ")
      )
    } else {
      return(NA_character_)
    }
  }

  hit[1]
}

read_csv_clean <- function(path, label) {
  if (!file.exists(path)) {
    stop2("No existe %s: %s", label, path)
  }

  x <- read_csv(path, show_col_types = FALSE, guess_max = 100000)
  names(x) <- clean_names(names(x))
  x
}

audit_ids <- function(df, id_col, label) {
  x <- as.character(df[[id_col]])

  tibble(
    source = label,
    id_col = id_col,
    n_rows = nrow(df),
    n_non_missing = sum(!is.na(x) & nzchar(x)),
    n_missing_or_blank = sum(is.na(x) | !nzchar(x)),
    n_unique = n_distinct(x[!is.na(x) & nzchar(x)]),
    n_duplicated_values = sum(duplicated(x[!is.na(x) & nzchar(x)]))
  )
}

prepare_aux <- function(df, source_label, requested_cols) {
  nms <- names(df)

  sample_col <- first_existing(
    nms,
    c("sample_id", ".sample_id", "id_metagenomics", "run_accession"),
    label = paste0("sample_id in ", source_label)
  )

  keep_cols <- intersect(requested_cols, nms)

  missing_cols <- setdiff(requested_cols, nms)

  if (length(missing_cols) > 0) {
    message2(
      "Aviso: en %s faltan columnas solicitadas: %s",
      source_label,
      paste(missing_cols, collapse = ", ")
    )
  }

  out <- df %>%
    transmute(
      sample_id_join = as.character(.data[[sample_col]]),
      across(all_of(keep_cols))
    ) %>%
    filter(!is.na(sample_id_join), nzchar(sample_id_join)) %>%
    distinct(sample_id_join, .keep_all = TRUE)

  attr(out, "sample_col") <- sample_col
  attr(out, "missing_cols") <- missing_cols

  out
}

safe_overwrite_column <- function(base, aux, join_col = "sample_id_join", cols) {
  cols_present <- intersect(cols, names(aux))

  if (length(cols_present) == 0) {
    return(base)
  }

  for (cc in cols_present) {
    base[[cc]] <- aux[[cc]][match(base[[join_col]], aux[[join_col]])]
  }

  base
}

# ------------------------------------------------------------------------------
# Args and paths
# ------------------------------------------------------------------------------

args <- parse_args()

overwrite <- as.integer(args$overwrite) == 1
stop_if_duplicate_ids <- as.integer(args$stop_if_duplicate_ids) == 1
expected_n <- suppressWarnings(as.integer(args$expected_n))
expected_profiles <- suppressWarnings(as.integer(args$expected_profiles))
if (!is.finite(expected_n) || !is.finite(expected_profiles)) {
  stop2("expected_n y expected_profiles deben ser enteros")
}
if (expected_n != 51L || expected_profiles != 17L) {
  stop2("El canon esta congelado en 51 muestras y 17 perfiles; no se permite redefinirlo")
}

run_id <- paste0("003_integrar_ejes_ECI_HI_en_metadata_", timestamp_now())
out_dir <- file.path(args$out_root, run_id)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

latest_file <- file.path(args$out_root, "LATEST_003_integrar_ejes_ECI_HI_en_metadata.txt")
writeLines(
  out_dir,
  file.path(args$out_root, "LATEST_ATTEMPT_003_integrar_ejes_ECI_HI_en_metadata.txt")
)

log_file <- file.path(out_dir, "logs", "003_log.txt")
sink(log_file, split = TRUE)

on.exit({
  sink()
}, add = TRUE)

message2("===== 003 Update shared metadata with ECI / HI / environmental axes =====")
message2("Start: %s", as.character(Sys.time()))
message2("Base CSV:     %s", args$base_csv)
message2("ECI CSV:      %s", args$eci_csv)
message2("Env axes CSV: %s", args$env_axes_csv)
message2("RBMA CSV:     %s", args$rbma_csv)
message2("MHI CSV:      %s", args$mhi_csv)
message2("Output dir:   %s", out_dir)
message2("Overwrite:    %s", overwrite)

# ------------------------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------------------------

base <- read_csv_clean(args$base_csv, "base_csv")
eci <- read_csv_clean(args$eci_csv, "eci_csv")
env_axes <- read_csv_clean(args$env_axes_csv, "env_axes_csv")
rbma <- read_csv_clean(args$rbma_csv, "rbma_csv")
mhi <- read_csv_clean(args$mhi_csv, "mhi_csv")

# ECI y los ejes por bloque deben proceder del mismo archivo generado por 101.
# Así no puede combinarse accidentalmente un ECI heredado del universo de 57
# muestras con ejes ambientales recalculados sobre las 51 muestras canónicas.
eci_path_norm <- normalizePath(args$eci_csv, winslash = "/", mustWork = TRUE)
env_axes_path_norm <- normalizePath(
  args$env_axes_csv,
  winslash = "/",
  mustWork = TRUE
)
same_101_metadata_file <- identical(eci_path_norm, env_axes_path_norm)

eci_tables_dir <- dirname(eci_path_norm)
eci_gate_path <- file.path(eci_tables_dir, "101_eci_canon_gate.csv")
axis_provenance_path <- file.path(
  eci_tables_dir,
  "101_axis_provenance.csv"
)

if (!same_101_metadata_file) {
  stop2(
    "eci_csv y env_axes_csv deben ser exactamente el mismo 101_metadata_with_env_axes.csv"
  )
}
if (!file.exists(eci_gate_path) || !file.exists(axis_provenance_path)) {
  stop2(
    paste0(
      "Falta la procedencia ECI-51 junto al metadata 101. Esperaba: ",
      eci_gate_path, " y ", axis_provenance_path
    )
  )
}

eci_gate_101 <- read_csv_clean(eci_gate_path, "101_eci_canon_gate.csv")
axis_provenance_101 <- read_csv_clean(
  axis_provenance_path,
  "101_axis_provenance.csv"
)

missing_eci_gate_cols <- setdiff(c("gate", "pass"), names(eci_gate_101))
missing_provenance_cols <- setdiff(
  c("axis", "n_samples_fitted", "analytical_universe"),
  names(axis_provenance_101)
)
if (length(missing_eci_gate_cols) > 0L || length(missing_provenance_cols) > 0L) {
  stop2(
    "Archivos de procedencia 101 incompletos; gate faltantes: %s; provenance faltantes: %s",
    paste(missing_eci_gate_cols, collapse = ","),
    paste(missing_provenance_cols, collapse = ",")
  )
}

required_eci_axes_lower <- c("eci_pc1", "eci_pc2", "eci_pc3")
legacy_eci_aliases <- c(
  "env_global_pc1", "env_global_pc2", "env_global_pc3"
)
forbidden_eci_aliases_101 <- intersect(
  legacy_eci_aliases,
  names(eci)
)
provenance_eci <- axis_provenance_101 %>%
  filter(tolower(axis) %in% required_eci_axes_lower)

eci_contract_audit <- tibble(
  gate = c(
    "eci_and_block_axes_same_101_file",
    "all_101_eci_gates_pass",
    "canonical_ECI_columns_present",
    "no_env_global_aliases_in_101_source",
    "three_ECI_axes_in_provenance",
    "ECI_provenance_sample_count_51",
    "ECI_provenance_universe_canon51"
  ),
  pass = c(
    same_101_metadata_file,
    nrow(eci_gate_101) > 0L && all(eci_gate_101$pass %in% TRUE),
    all(required_eci_axes_lower %in% names(eci)),
    length(forbidden_eci_aliases_101) == 0L,
    nrow(provenance_eci) == 3L,
    nrow(provenance_eci) == 3L &&
      all(suppressWarnings(as.integer(provenance_eci$n_samples_fitted)) == expected_n),
    nrow(provenance_eci) == 3L &&
      all(
        provenance_eci$analytical_universe ==
          "canonical_51_samples_17_complete_profiles"
      )
  ),
  severity = "critical",
  detail = c(
    eci_path_norm,
    paste0("n_101_gates=", nrow(eci_gate_101)),
    paste(intersect(required_eci_axes_lower, names(eci)), collapse = ";"),
    ifelse(
      length(forbidden_eci_aliases_101) == 0L,
      "0 forbidden aliases",
      paste(forbidden_eci_aliases_101, collapse = ";")
    ),
    paste(provenance_eci$axis, collapse = ";"),
    paste(provenance_eci$n_samples_fitted, collapse = ";"),
    paste(unique(provenance_eci$analytical_universe), collapse = ";")
  )
)
write_csv(
  eci_contract_audit,
  file.path(out_dir, "tables", "003_compuerta_procedencia_ECI_canon_51.csv")
)
if (!all(eci_contract_audit$pass)) {
  stop2(
    paste0(
      "Fallo la compuerta de procedencia ECI-51; revise ",
      "003_compuerta_procedencia_ECI_canon_51.csv"
    )
  )
}

# ------------------------------------------------------------------------------
# Check base ID
# ------------------------------------------------------------------------------

base_id_col <- first_existing(
  names(base),
  c("sample_id", "id_metagenomics", "run_accession"),
  label = "sample ID in canonical base metadata"
)

base <- base %>%
  mutate(
    .row_order_original = row_number(),
    sample_id_join = as.character(.data[[base_id_col]])
  )

id_audit <- bind_rows(
  audit_ids(base, "sample_id_join", "base_metadata"),
  audit_ids(eci, first_existing(names(eci), c("sample_id", ".sample_id", "id_metagenomics", "run_accession"), label = "ECI ID"), "eci_csv"),
  audit_ids(env_axes, first_existing(names(env_axes), c("sample_id", ".sample_id", "id_metagenomics", "run_accession"), label = "env_axes ID"), "env_axes_csv"),
  audit_ids(rbma, first_existing(names(rbma), c("sample_id", ".sample_id", "id_metagenomics", "run_accession"), label = "RBMA ID"), "rbma_csv"),
  audit_ids(mhi, first_existing(names(mhi), c("sample_id", ".sample_id", "id_metagenomics", "run_accession"), label = "MHI ID"), "mhi_csv")
)

write_csv(
  id_audit,
  file.path(out_dir, "tables", "003_id_audit.csv")
)

print(id_audit)

if (stop_if_duplicate_ids && any(id_audit$n_duplicated_values > 0)) {
  stop2("Hay IDs duplicados en alguna fuente. Revisa 003_id_audit.csv antes de actualizar.")
}

if (nrow(base) != expected_n || n_distinct(base$sample_id_join) != expected_n) {
  stop2("La metadata base no representa exactamente el canon de %s muestras", expected_n)
}

profile_col_base <- first_existing(names(base), c("lat_block", "profile_id"), label = "profile ID in base")
depth_col_base <- first_existing(names(base), c("depth_cm"), label = "depth in base")
profile_check <- base %>%
  transmute(
    sample_id = sample_id_join,
    profile_id = as.character(.data[[profile_col_base]]),
    depth_cm = suppressWarnings(as.numeric(.data[[depth_col_base]]))
  ) %>%
  group_by(profile_id) %>%
  summarise(
    n = n(),
    depths = paste(sort(unique(depth_cm)), collapse = "|"),
    valid = n == 3L && setequal(depth_cm, c(5, 20, 40)),
    .groups = "drop"
  )
write_csv(profile_check, file.path(out_dir, "tables", "003_auditoria_perfiles_canon_51.csv"))
if (nrow(profile_check) != expected_profiles || !all(profile_check$valid)) {
  stop2("La metadata base no representa exactamente %s perfiles completos", expected_profiles)
}

# ------------------------------------------------------------------------------
# Recode Preserved -> Conserved
# ------------------------------------------------------------------------------

stage_cols <- intersect(c("stage", "condition", "restoration4"), names(base))

for (cc in stage_cols) {
  base[[cc]] <- ifelse(
    as.character(base[[cc]]) == "Preserved",
    "Conserved",
    as.character(base[[cc]])
  )
}

# If restoration4 does not exist, create it from stage if possible.
if (!"restoration4" %in% names(base) && "stage" %in% names(base)) {
  base$restoration4 <- base$stage
}

# ------------------------------------------------------------------------------
# Prepare auxiliary tables
# ------------------------------------------------------------------------------

eci_cols <- c("eci_pc1", "eci_pc2", "eci_pc3")

block_axis_cols <- c(
  "vegetation_landscape_pc1",
  "water_inundation_pc1",
  "moisture_stress_pc1",
  "physicochemical_pc1",
  "nutrients_redox_pc1"
)

rbma_candidates <- c("rbma1", "rbma1_raw", "hi")

eci_aux <- prepare_aux(
  eci,
  source_label = "ECI metadata",
  requested_cols = eci_cols
)
if (length(attr(eci_aux, "missing_cols")) > 0L) {
  stop2(
    "El metadata 101 no contiene los tres ECI canónicos recalculados: %s",
    paste(attr(eci_aux, "missing_cols"), collapse = ", ")
  )
}

env_axes_aux <- prepare_aux(
  env_axes,
  source_label = "environmental axes metadata",
  requested_cols = block_axis_cols
)
if (length(attr(env_axes_aux, "missing_cols")) > 0L) {
  stop2(
    "El metadata 101 no contiene todos los ejes ambientales por bloque: %s",
    paste(attr(env_axes_aux, "missing_cols"), collapse = ", ")
  )
}

rbma_sample_col <- first_existing(
  names(rbma),
  c("sample_id", ".sample_id", "id_metagenomics", "run_accession"),
  label = "sample_id in RBMA file"
)

rbma_col <- first_existing(
  names(rbma),
  rbma_candidates,
  label = "RBMA1 / HI column in RBMA file"
)

rbma_aux <- rbma %>%
  transmute(
    sample_id_join = as.character(.data[[rbma_sample_col]]),
    rbma1 = suppressWarnings(as.numeric(.data[[rbma_col]])),
    hi = suppressWarnings(as.numeric(.data[[rbma_col]]))
  ) %>%
  filter(!is.na(sample_id_join), nzchar(sample_id_join)) %>%
  distinct(sample_id_join, .keep_all = TRUE)

mhi_sample_col <- first_existing(
  names(mhi),
  c("sample_id", ".sample_id", "id_metagenomics", "run_accession"),
  label = "sample_id in MHI file"
)
mhi_local_col <- first_existing(
  names(mhi),
  c("mhi_local", "mhi"),
  label = "MHI_local in MHI file"
)
mhi_global_col <- first_existing(
  names(mhi),
  c("mhi_global"),
  required = FALSE,
  label = "MHI_global in MHI file"
)
mhi_aux <- mhi %>%
  transmute(
    sample_id_join = as.character(.data[[mhi_sample_col]]),
    mhi_local = suppressWarnings(as.numeric(.data[[mhi_local_col]])),
    mhi_global = if (!is.na(mhi_global_col)) {
      suppressWarnings(as.numeric(.data[[mhi_global_col]]))
    } else {
      NA_real_
    }
  ) %>%
  filter(!is.na(sample_id_join), nzchar(sample_id_join)) %>%
  distinct(sample_id_join, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# Merge while preserving strict base order
# ------------------------------------------------------------------------------

updated <- base

updated <- safe_overwrite_column(
  updated,
  eci_aux,
  join_col = "sample_id_join",
  cols = eci_cols
)

updated <- safe_overwrite_column(
  updated,
  env_axes_aux,
  join_col = "sample_id_join",
  cols = block_axis_cols
)

updated <- safe_overwrite_column(
  updated,
  rbma_aux,
  join_col = "sample_id_join",
  cols = c("rbma1", "hi")
)

updated <- safe_overwrite_column(
  updated,
  mhi_aux,
  join_col = "sample_id_join",
  cols = c("mhi_local", "mhi_global")
)

# El metadata base 044 puede contener ejes globales heredados calculados sobre
# otro universo analitico. Se conserva el archivo base como snapshot, pero esos
# aliases ambiguos no pueden coexistir con los ECI_PC1-3 canonicos de 101 en el
# producto integrado.
legacy_eci_aliases_in_base <- intersect(legacy_eci_aliases, names(base))
legacy_eci_aliases_before_cleanup <- intersect(
  legacy_eci_aliases,
  names(updated)
)
legacy_eci_aliases_removed <- legacy_eci_aliases_before_cleanup

updated <- updated %>%
  select(-any_of(legacy_eci_aliases))

legacy_eci_aliases_after_cleanup <- intersect(
  legacy_eci_aliases,
  names(updated)
)

legacy_alias_cleanup_audit <- tibble(
  alias = legacy_eci_aliases,
  present_in_base = alias %in% legacy_eci_aliases_in_base,
  present_before_cleanup = alias %in% legacy_eci_aliases_before_cleanup,
  removed_from_integrated_metadata = alias %in% legacy_eci_aliases_removed,
  present_after_cleanup = alias %in% legacy_eci_aliases_after_cleanup
)

write_csv(
  legacy_alias_cleanup_audit,
  file.path(out_dir, "tables", "003_legacy_ECI_alias_cleanup_audit.csv")
)

cleanup_gate <- tibble(
  gate = c(
    "legacy_env_global_alias_cleanup_consistent",
    "no_env_global_aliases_in_integrated_metadata"
  ),
  pass = c(
    setequal(
      legacy_eci_aliases_before_cleanup,
      legacy_eci_aliases_removed
    ) && length(legacy_eci_aliases_after_cleanup) == 0L,
    length(legacy_eci_aliases_after_cleanup) == 0L
  ),
  severity = "critical",
  detail = c(
    paste0(
      "detected=",
      ifelse(
        length(legacy_eci_aliases_before_cleanup) == 0L,
        "none",
        paste(legacy_eci_aliases_before_cleanup, collapse = ";")
      ),
      ";removed=",
      ifelse(
        length(legacy_eci_aliases_removed) == 0L,
        "none",
        paste(legacy_eci_aliases_removed, collapse = ";")
      )
    ),
    ifelse(
      length(legacy_eci_aliases_after_cleanup) == 0L,
      "0 forbidden aliases",
      paste(legacy_eci_aliases_after_cleanup, collapse = ";")
    )
  )
)

eci_contract_audit <- bind_rows(eci_contract_audit, cleanup_gate)
write_csv(
  eci_contract_audit,
  file.path(out_dir, "tables", "003_compuerta_procedencia_ECI_canon_51.csv")
)
if (!all(eci_contract_audit$pass)) {
  stop2(
    paste0(
      "Fallo la limpieza de aliases ECI heredados; revise ",
      "003_compuerta_procedencia_ECI_canon_51.csv"
    )
  )
}

message2(
  "Legacy ECI aliases removed: %s",
  ifelse(
    length(legacy_eci_aliases_removed) == 0L,
    "none",
    paste(legacy_eci_aliases_removed, collapse = ", ")
  )
)

# Rename columns to preferred manuscript/shared naming.
rename_if_present <- function(df, old, new) {
  if (old %in% names(df)) {
    names(df)[names(df) == old] <- new
  }
  df
}

updated <- updated %>%
  rename_if_present("eci_pc1", "ECI_PC1") %>%
  rename_if_present("eci_pc2", "ECI_PC2") %>%
  rename_if_present("eci_pc3", "ECI_PC3") %>%
  rename_if_present("vegetation_landscape_pc1", "vegetation_landscape_PC1") %>%
  rename_if_present("water_inundation_pc1", "water_inundation_PC1") %>%
  rename_if_present("moisture_stress_pc1", "moisture_stress_PC1") %>%
  rename_if_present("physicochemical_pc1", "physicochemical_PC1") %>%
  rename_if_present("nutrients_redox_pc1", "nutrients_redox_PC1") %>%
  rename_if_present("rbma1", "RBMA1") %>%
  rename_if_present("hi", "HI") %>%
  rename_if_present("mhi_local", "MHI_local") %>%
  rename_if_present("mhi_global", "MHI_global")

updated <- updated %>%
  arrange(.row_order_original)

# ------------------------------------------------------------------------------
# Matching audit
# ------------------------------------------------------------------------------

columns_added <- c(
  "ECI_PC1", "ECI_PC2", "ECI_PC3",
  "vegetation_landscape_PC1",
  "water_inundation_PC1",
  "moisture_stress_PC1",
  "physicochemical_PC1",
  "nutrients_redox_PC1",
  "RBMA1", "HI", "MHI_local", "MHI_global"
)

columns_added_present <- intersect(columns_added, names(updated))

matching_audit <- tibble(
  column = columns_added_present,
  n_non_missing = sapply(columns_added_present, function(cc) sum(!is.na(updated[[cc]]))),
  n_missing = sapply(columns_added_present, function(cc) sum(is.na(updated[[cc]]))),
  pct_missing = round(100 * sapply(columns_added_present, function(cc) mean(is.na(updated[[cc]]))), 3),
  min = sapply(columns_added_present, function(cc) suppressWarnings(min(updated[[cc]], na.rm = TRUE))),
  median = sapply(columns_added_present, function(cc) suppressWarnings(median(updated[[cc]], na.rm = TRUE))),
  max = sapply(columns_added_present, function(cc) suppressWarnings(max(updated[[cc]], na.rm = TRUE)))
)

write_csv(
  matching_audit,
  file.path(out_dir, "tables", "003_added_columns_missingness_summary.csv")
)

print(matching_audit)

# IDs missing by source
source_match_audit <- tibble(
  source = c("ECI metadata", "environmental axes metadata", "RBMA / HI metadata", "MHI metadata"),
  n_base_rows = nrow(base),
  n_matched = c(
    sum(base$sample_id_join %in% eci_aux$sample_id_join),
    sum(base$sample_id_join %in% env_axes_aux$sample_id_join),
    sum(base$sample_id_join %in% rbma_aux$sample_id_join),
    sum(base$sample_id_join %in% mhi_aux$sample_id_join)
  ),
  n_unmatched = n_base_rows - n_matched
)

write_csv(
  source_match_audit,
  file.path(out_dir, "tables", "003_source_matching_summary.csv")
)

print(source_match_audit)

if (nrow(base) != expected_n || n_distinct(base$sample_id_join) != expected_n) {
  stop2("La base integrada debe conservar exactamente %s muestras", expected_n)
}
if (any(source_match_audit$n_matched != expected_n)) {
  stop2("Alguna fuente no cubre las %s muestras canonicas; revise 003_source_matching_summary.csv", expected_n)
}

# La coincidencia de IDs no basta: cada variable analitica critica debe existir
# y contener un valor numerico finito para las 51 muestras. MHI_global se deja
# fuera de esta compuerta porque es una salida opcional del archivo MHI.
required_integrated_cols <- c(
  "ECI_PC1", "ECI_PC2", "ECI_PC3",
  "vegetation_landscape_PC1",
  "water_inundation_PC1",
  "moisture_stress_PC1",
  "physicochemical_PC1",
  "nutrients_redox_PC1",
  "RBMA1", "HI", "MHI_local"
)
missing_integrated_cols <- setdiff(required_integrated_cols, names(updated))
if (length(missing_integrated_cols) > 0L) {
  stop2(
    "Faltan columnas integradas criticas: %s",
    paste(missing_integrated_cols, collapse = ", ")
  )
}

forbidden_eci_aliases_final <- intersect(
  legacy_eci_aliases,
  clean_names(names(updated))
)
if (length(forbidden_eci_aliases_final) > 0L) {
  stop2(
    "Persisten aliases ECI heredados en el metadata integrado: %s",
    paste(forbidden_eci_aliases_final, collapse = ", ")
  )
}

integrated_completeness <- tibble(
  column = required_integrated_cols,
  n_rows = nrow(updated),
  n_numeric_finite = vapply(
    required_integrated_cols,
    function(cc) sum(is.finite(suppressWarnings(as.numeric(updated[[cc]])))),
    integer(1)
  ),
  n_missing_or_nonfinite = n_rows - n_numeric_finite,
  complete_canon_51 = n_numeric_finite == expected_n
)
write_csv(
  integrated_completeness,
  file.path(out_dir, "tables", "003_compuerta_completitud_variables_canon_51.csv")
)
if (!all(integrated_completeness$complete_canon_51)) {
  stop2(
    "Hay variables analiticas incompletas; revise 003_compuerta_completitud_variables_canon_51.csv"
  )
}

missing_by_source <- bind_rows(
  tibble(
    source = "ECI metadata",
    id_metagenomics = setdiff(base$sample_id_join, eci_aux$sample_id_join)
  ),
  tibble(
    source = "environmental axes metadata",
    id_metagenomics = setdiff(base$sample_id_join, env_axes_aux$sample_id_join)
  ),
  tibble(
    source = "RBMA / HI metadata",
    id_metagenomics = setdiff(base$sample_id_join, rbma_aux$sample_id_join)
  ),
  tibble(
    source = "MHI metadata",
    id_metagenomics = setdiff(base$sample_id_join, mhi_aux$sample_id_join)
  )
)

write_csv(
  missing_by_source,
  file.path(out_dir, "tables", "003_unmatched_ids_by_source.csv")
)

# ------------------------------------------------------------------------------
# Column order and output
# ------------------------------------------------------------------------------

updated_out <- updated %>%
  select(-sample_id_join, -.row_order_original)

# Force sample_id first. Keep id_metagenomics only as an optional alias.
if ("sample_id" %in% names(updated_out)) {
  updated_out <- updated_out %>%
    select(sample_id, everything())
} else {
  updated_out <- updated_out %>%
    mutate(sample_id = as.character(.data[[base_id_col]]), .before = 1)
}

# Decide output path
if (is.na(args$out_csv) || !nzchar(args$out_csv)) {
  out_csv <- file.path(out_dir, "tables", "003_metadata_integrada_canon_51.csv")
} else {
  out_csv <- args$out_csv
}

# Snapshot reproducible de la base dentro del mismo output; no se escribe en
# directorios historicos de inputs.
backup_csv <- file.path(out_dir, "tables", "003_base_canon_51_snapshot.csv")
backup_ok <- file.copy(args$base_csv, backup_csv, overwrite = TRUE)
if (!isTRUE(backup_ok)) {
  stop2("No se pudo crear el snapshot reproducible: %s", backup_csv)
}

write_csv(updated_out, out_csv)

message2("Updated metadata written: %s", out_csv)

# Save a copy in results audit dir too
write_csv(
  updated_out,
  file.path(out_dir, "tables", "003_updated_metadata_preview_full.csv")
)

# Export column names and first rows
write_csv(
  tibble(
    col_index = seq_along(names(updated_out)),
    column = names(updated_out)
  ),
  file.path(out_dir, "tables", "003_updated_column_names.csv")
)

write_csv(
  updated_out %>% slice_head(n = 20),
  file.path(out_dir, "tables", "003_updated_first_20_rows.csv")
)

# ------------------------------------------------------------------------------
# Metadata note
# ------------------------------------------------------------------------------

note <- c(
  "003 integrar ejes ECI HI en metadata",
  "==========================",
  "",
  paste0("Base CSV: ", args$base_csv),
  paste0("Output CSV: ", out_csv),
  paste0("Backup CSV: ", backup_csv),
  paste0("ECI source: ", args$eci_csv),
  paste0("Environmental axes source: ", args$env_axes_csv),
  paste0("RBMA / HI source: ", args$rbma_csv),
  paste0("MHI source: ", args$mhi_csv),
  "",
  "Matching rule:",
  "  base$sample_id = auxiliary$sample_id (aliases accepted and audited)",
  "  analytical universe = exactly 51 samples in 17 complete depth profiles",
  "",
  "Stage recoding:",
  "  Preserved -> Conserved",
  "",
  "Legacy ECI cleanup:",
  paste0(
    "  Removed from integrated metadata: ",
    ifelse(
      length(legacy_eci_aliases_removed) == 0L,
      "none",
      paste(legacy_eci_aliases_removed, collapse = ", ")
    )
  ),
  "  The original base is preserved unchanged in the snapshot file.",
  "",
  "Added / updated columns:",
  paste0("  - ", columns_added_present),
  "",
  "Environmental summary axes:",
  "  vegetation_landscape_PC1: shdi, ndvi_p50, ndvi_sd",
  "  water_inundation_PC1: ndwi_p50, ndwi_sd, mndwi_p50, mndwi_sd, water_prop_mndwi",
  "  moisture_stress_PC1: ndmi_p50, ndmi_sd",
  "  physicochemical_PC1: temperature_c, salinity_ups, ph, redox_shallow_mv",
  "  nutrients_redox_PC1: n_no2_umol_l, n_no3_umol_l, n_nh3_mg_l_campo, n_nh4plus_umol_l, p_po4_3_umol_l, s_2_umol",
  "",
  "Index definitions:",
  "  ECI_PC1-3: global environmental PCA refitted exclusively on the canonical 51 samples; no microbiome used.",
  "  PCA axis signs are arbitrary and must be interpreted with the 101 loadings.",
  "  MHI_local: joined from the canonical 20260216 MHI run and restricted to the 51-sample universe; it is not refitted here.",
  "  HI / RBMA1: microbiome-derived residual bimodal axis; higher values indicate more Preserved-like microbiome."
)

writeLines(
  note,
  file.path(out_dir, "003_update_note.txt")
)

# ------------------------------------------------------------------------------
# Final summary
# ------------------------------------------------------------------------------

run_summary <- tibble(
  key = c(
    "script",
    "timestamp",
    "base_csv",
    "out_csv",
    "backup_csv",
    "out_dir",
    "eci_gate_path",
    "axis_provenance_path",
    "n_rows_original",
    "n_rows_updated",
    "n_cols_original",
    "n_cols_updated",
    "n_profiles",
    "legacy_eci_aliases_detected_in_base",
    "legacy_eci_aliases_removed",
    "legacy_eci_aliases_remaining",
    "overwrite",
    "latest_file"
  ),
  value = c(
    "003_integrar_ejes_ECI_HI_en_metadata.R",
    as.character(Sys.time()),
    args$base_csv,
    out_csv,
    backup_csv,
    out_dir,
    eci_gate_path,
    axis_provenance_path,
    as.character(nrow(base)),
    as.character(nrow(updated_out)),
    as.character(ncol(base) - 2),
    as.character(ncol(updated_out)),
    as.character(nrow(profile_check)),
    ifelse(
      length(legacy_eci_aliases_in_base) == 0L,
      "none",
      paste(legacy_eci_aliases_in_base, collapse = ";")
    ),
    ifelse(
      length(legacy_eci_aliases_removed) == 0L,
      "none",
      paste(legacy_eci_aliases_removed, collapse = ";")
    ),
    ifelse(
      length(legacy_eci_aliases_after_cleanup) == 0L,
      "none",
      paste(legacy_eci_aliases_after_cleanup, collapse = ";")
    ),
    as.character(overwrite),
    latest_file
  )
)

write_csv(
  run_summary,
  file.path(out_dir, "tables", "003_run_summary.csv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "logs", "003_sessionInfo.txt")
)

# LATEST se actualiza solo después de completar todas las compuertas y outputs.
writeLines(out_dir, latest_file)

message2("===== DONE 003 =====")
message2("Output CSV: %s", out_csv)
message2("Backup CSV: %s", backup_csv)
message2("Audit dir:   %s", out_dir)
message2("Latest:      %s", latest_file)

# ==============================================================================
# End
# ==============================================================================
