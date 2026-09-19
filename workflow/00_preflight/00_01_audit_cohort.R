#!/usr/bin/env Rscript
# ============================================================
# 00_01_audit_cohort.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Cohort identity, complete profiles and reference coverage; phyloseq object required.
# Inputs (source expressions; complete list in docs/contracts/00_01_audit_cohort.json):
#   first <- readLines(path, n = 1L, warn = FALSE)
#   read_delim(path, delim = delim, show_col_types = FALSE, guess_max = 100000)
#   ps_source <- readRDS(PS_RDS)
# Outputs (source expressions; complete list in contract):
#   writeLines(OUT_DIR, file.path(OUT_ROOT, paste0("LATEST_ATTEMPT_", SCRIPT_ID, ".txt")))
#   write_csv(tibble(missing_column = missing_required),
#   write_csv(column_audit, file.path(TABLE_DIR, "001_columnas_y_faltantes.csv"))
#   write_csv(invalid_rows, file.path(TABLE_DIR, "001_filas_con_diseno_invalido.csv"))
#   write_csv(duplicate_samples, file.path(TABLE_DIR, "001_ids_muestra_duplicados.csv"))
#   write_csv(profile_audit, file.path(TABLE_DIR, "001_auditoria_perfiles.csv"))
#   write_csv(design_counts, file.path(TABLE_DIR, "001_conteos_diseno_localidad_estadio.csv"))
#   write_csv(reference_audit, file.path(TABLE_DIR, "001_disponibilidad_referencias_localidad_profundidad.csv"))
#   write_csv(ps_matching, file.path(TABLE_DIR, "001_matching_metadata_phyloseq.csv"))
#   write_csv(md_canon, file.path(TABLE_DIR, "001_metadata_canon_51.csv"))
# Algorithmic provenance:
# Cohort membership, identifier integrity and complete three-depth profile checks.
#   McMurdie & Holmes (2013), phyloseq, doi:10.1371/journal.pone.0061217.
# Source SHA-256: a2d30465616b3540fb025a455cb1cb9f3b1ade2a66860dc14524e2f7a846ee00
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


# 001_auditar_estructura_metadata.R
# Auditoria del diseno antes de cualquier contraste estadistico.
# La unidad experimental primaria es el perfil (tres profundidades anidadas).

suppressPackageStartupMessages({
  library(phyloseq)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (!length(i)) return(default)
  if (i == length(args)) stop("Argumento sin valor: ", flag, call. = FALSE)
  args[[i + 1L]]
}
as_int <- function(x, default) {
  z <- suppressWarnings(as.integer(x)); if (is.na(z)) default else z
}
stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S", tz = "UTC")

IN_CSV <- get_arg(
  "--in_csv",
  "/home/fjbalvino/Tipping_points/Results/044_build_protein_filtered_matrix_HI_20260529_000843/tables/044_metadata_aligned_51samples.csv"
)
PS_RDS <- get_arg(
  "--ps_rds",
  "/home/fjbalvino/R/Ciencia-Frontera/resultados/phyloseq_withConservedDegraded_20260130_1730/rank_S/phyloseq_FILTRADO_rank_S.rds"
)
OUT_ROOT <- get_arg("--out_root", "/home/fjbalvino/Tipping_points/resultados_finales")
SAMPLE_COL <- get_arg("--sample_col", "sample_id")
PROFILE_COL <- get_arg("--profile_col", "lat_block")
LOCALITY_COL <- get_arg("--locality_col", "locality")
STAGE_COL <- get_arg("--stage_col", "restoration4")
DEPTH_COL <- get_arg("--depth_col", "depth_cm")
EXPECTED_PROFILES <- as_int(get_arg("--expected_profiles", "17"), 17L)
EXPECTED_SAMPLES <- as_int(get_arg("--expected_samples", "51"), 51L)
if (EXPECTED_SAMPLES != 51L || EXPECTED_PROFILES != 17L) {
  stop("El canon esta congelado en 51 muestras y 17 perfiles; no se permite redefinirlo.", call. = FALSE)
}

SCRIPT_ID <- "001_auditar_estructura_metadata"
OUT_DIR <- file.path(OUT_ROOT, paste0(SCRIPT_ID, "_", stamp()))
TABLE_DIR <- file.path(OUT_DIR, "tables")
LOG_DIR <- file.path(OUT_DIR, "logs")
RDS_DIR <- file.path(OUT_DIR, "rds")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
LATEST_FILE <- file.path(OUT_ROOT, paste0("LATEST_", SCRIPT_ID, ".txt"))
writeLines(OUT_DIR, file.path(OUT_ROOT, paste0("LATEST_ATTEMPT_", SCRIPT_ID, ".txt")))

log_con <- file(file.path(LOG_DIR, "001_log.txt"), "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0L) sink(type = "message")
  while (sink.number() > 0L) sink()
  close(log_con)
}, add = TRUE)

if (!file.exists(IN_CSV)) stop("No existe metadata canonica: ", IN_CSV, call. = FALSE)
if (!file.exists(PS_RDS)) stop("No existe phyloseq fuente: ", PS_RDS, call. = FALSE)

STAGE_LEVELS <- c(
  "Degraded", "Early restoration", "Intermediate restoration",
  "Advanced restoration", "Conserved"
)
DEPTH_LEVELS <- c(5, 20, 40)

normalize_depth <- function(x) {
  z <- str_extract(as.character(x), "[0-9]+(?:\\.[0-9]+)?")
  suppressWarnings(as.numeric(z))
}
normalize_stage <- function(x) {
  z <- str_squish(as.character(x))
  key <- str_to_lower(z)
  case_when(
    key %in% c("degraded", "degradado") ~ "Degraded",
    str_detect(key, "early|tempran") ~ "Early restoration",
    str_detect(key, "intermediate|intermedi") ~ "Intermediate restoration",
    str_detect(key, "advanced|avanzad") ~ "Advanced restoration",
    key %in% c("conserved", "preserved", "conservado", "preservado") ~ "Conserved",
    TRUE ~ z
  )
}
safe_read <- function(path) {
  first <- readLines(path, n = 1L, warn = FALSE)
  counts <- c(
    comma = str_count(first, fixed(",")),
    tab = str_count(first, fixed("\t")),
    semicolon = str_count(first, fixed(";")),
    pipe = str_count(first, fixed("|"))
  )
  delim <- c(comma = ",", tab = "\t", semicolon = ";", pipe = "|")[[names(which.max(counts))]]
  read_delim(path, delim = delim, show_col_types = FALSE, guess_max = 100000)
}

cat("001 metadata/design audit\nInput:", IN_CSV, "\nOutput:", OUT_DIR, "\n")
md_raw <- safe_read(IN_CSV)

required <- c(SAMPLE_COL, PROFILE_COL, LOCALITY_COL, STAGE_COL, DEPTH_COL)
missing_required <- setdiff(required, names(md_raw))
if (length(missing_required)) {
  write_csv(tibble(missing_column = missing_required),
            file.path(TABLE_DIR, "001_missing_required_columns.csv"))
  stop("Faltan columnas de diseno: ", paste(missing_required, collapse = ", "), call. = FALSE)
}

column_audit <- tibble(
  column = names(md_raw),
  class = vapply(md_raw, function(x) paste(class(x), collapse = ";"), character(1)),
  n_missing = vapply(md_raw, function(x) sum(is.na(x) | !nzchar(trimws(as.character(x)))), integer(1)),
  n_unique = vapply(md_raw, dplyr::n_distinct, integer(1), na.rm = TRUE)
) %>% mutate(pct_missing = 100 * n_missing / nrow(md_raw))
write_csv(column_audit, file.path(TABLE_DIR, "001_columnas_y_faltantes.csv"))

md <- md_raw %>%
  transmute(
    sample_id = str_squish(as.character(.data[[SAMPLE_COL]])),
    profile_raw = str_squish(as.character(.data[[PROFILE_COL]])),
    locality = str_squish(as.character(.data[[LOCALITY_COL]])),
    restoration4 = normalize_stage(.data[[STAGE_COL]]),
    depth_cm = normalize_depth(.data[[DEPTH_COL]]),
    source_row = row_number()
  ) %>%
  mutate(
    profile_id = if_else(
      str_detect(profile_raw, fixed(locality)),
      profile_raw,
      paste(locality, profile_raw, sep = "::")
    ),
    restoration4 = factor(restoration4, levels = STAGE_LEVELS, ordered = FALSE),
    depth_cm = as.numeric(depth_cm)
  )

invalid_rows <- md %>%
  filter(
    is.na(sample_id) | !nzchar(sample_id) |
      is.na(profile_id) | !nzchar(profile_id) |
      is.na(locality) | !nzchar(locality) |
      is.na(restoration4) | is.na(depth_cm)
  )
write_csv(invalid_rows, file.path(TABLE_DIR, "001_filas_con_diseno_invalido.csv"))

duplicate_samples <- md %>% count(sample_id, name = "n") %>% filter(n != 1L)
write_csv(duplicate_samples, file.path(TABLE_DIR, "001_ids_muestra_duplicados.csv"))

profile_audit <- md %>%
  group_by(profile_id) %>%
  summarise(
    locality = paste(sort(unique(locality)), collapse = "|"),
    restoration4 = paste(sort(unique(as.character(restoration4))), collapse = "|"),
    n_samples = n(),
    n_sample_ids = n_distinct(sample_id),
    n_localities = n_distinct(locality),
    n_stages = n_distinct(restoration4),
    depths = paste(sort(unique(depth_cm)), collapse = "|"),
    n_depths = n_distinct(depth_cm),
    has_5 = 5 %in% depth_cm,
    has_20 = 20 %in% depth_cm,
    has_40 = 40 %in% depth_cm,
    valid_complete_profile =
      n_samples == 3L && n_sample_ids == 3L && n_localities == 1L &&
      n_stages == 1L && setequal(depth_cm, DEPTH_LEVELS),
    .groups = "drop"
  ) %>% arrange(locality, restoration4, profile_id)
write_csv(profile_audit, file.path(TABLE_DIR, "001_auditoria_perfiles.csv"))

design_counts <- md %>%
  filter(!is.na(restoration4), !is.na(locality)) %>%
  distinct(locality, restoration4, profile_id, sample_id, depth_cm) %>%
  group_by(locality, restoration4) %>%
  summarise(
    n_profiles = n_distinct(profile_id),
    n_samples = n_distinct(sample_id),
    n_depths_observed = n_distinct(depth_cm),
    profile_ids = paste(sort(unique(profile_id)), collapse = ";"),
    .groups = "drop"
  ) %>%
  complete(
    locality = sort(unique(md$locality[!is.na(md$locality)])),
    restoration4 = factor(STAGE_LEVELS, levels = STAGE_LEVELS),
    fill = list(n_profiles = 0L, n_samples = 0L, n_depths_observed = 0L, profile_ids = "")
  ) %>%
  mutate(
    empty_cell = n_profiles == 0L,
    single_profile_cell = n_profiles == 1L,
    inferential_warning = case_when(
      empty_cell ~ "empty_locality_stage_cell",
      single_profile_cell ~ "one_independent_profile",
      TRUE ~ "none"
    )
  )
write_csv(design_counts, file.path(TABLE_DIR, "001_conteos_diseno_localidad_estadio.csv"))

reference_audit <- md %>%
  filter(restoration4 %in% c("Conserved", "Degraded")) %>%
  group_by(locality, depth_cm, restoration4) %>%
  summarise(
    n_profiles = n_distinct(profile_id),
    n_samples = n_distinct(sample_id),
    profile_ids = paste(sort(unique(profile_id)), collapse = ";"),
    sample_ids = paste(sort(unique(sample_id)), collapse = ";"),
    .groups = "drop"
  ) %>%
  complete(
    locality = sort(unique(md$locality[!is.na(md$locality)])),
    depth_cm = DEPTH_LEVELS,
    restoration4 = factor(c("Degraded", "Conserved"), levels = STAGE_LEVELS),
    fill = list(n_profiles = 0L, n_samples = 0L, profile_ids = "", sample_ids = "")
  ) %>%
  mutate(
    reference_available = n_profiles >= 1L,
    loo_reference_possible = n_profiles >= 2L,
    primary_rule = case_when(
      restoration4 == "Conserved" & n_profiles == 1L ~
        "usable_as_anchor_for_nonconserved_only; exclude_self_from_inference",
      restoration4 == "Conserved" & n_profiles >= 2L ~
        "leave_one_profile_out_available",
      n_profiles == 0L ~ "missing_reference",
      TRUE ~ "available"
    )
  )
write_csv(reference_audit, file.path(TABLE_DIR, "001_disponibilidad_referencias_localidad_profundidad.csv"))

model_frame <- md %>%
  filter(!is.na(restoration4), !is.na(locality), !is.na(depth_cm)) %>%
  mutate(
    locality = factor(locality), restoration4 = droplevels(restoration4),
    depth_cm = factor(depth_cm)
  )
X <- model.matrix(~ locality + depth_cm + restoration4, data = model_frame)
rank_x <- qr(X)$rank

# El objeto de 57 muestras se usa unicamente como fuente taxonomica. La lista
# de IDs de IN_CSV define el universo analitico y se aplica antes de cualquier
# analisis de taxonomia.
ps_source <- readRDS(PS_RDS)
ps_source_ids <- as.character(phyloseq::sample_names(ps_source))
canonical_ids <- as.character(md$sample_id)
matching_ids <- union(canonical_ids, ps_source_ids)
ps_matching <- tibble(
  sample_id = matching_ids,
  in_canonical_metadata = matching_ids %in% canonical_ids,
  in_source_phyloseq = matching_ids %in% ps_source_ids
) %>%
  mutate(
    status = case_when(
      in_canonical_metadata & in_source_phyloseq ~ "retained_canon_51",
      in_canonical_metadata & !in_source_phyloseq ~ "canonical_missing_in_phyloseq",
      !in_canonical_metadata & in_source_phyloseq ~ "excluded_outside_canon_51",
      TRUE ~ "unexpected"
    )
  ) %>%
  arrange(status, sample_id)
write_csv(ps_matching, file.path(TABLE_DIR, "001_matching_metadata_phyloseq.csv"))

ps_keep <- intersect(ps_source_ids, canonical_ids)
ps_canon <- phyloseq::prune_samples(ps_keep, ps_source)
ps_canon <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_canon) > 0, ps_canon)

md_canon <- md_raw
md_canon$sample_id <- md$sample_id
md_canon$profile_id <- md$profile_id
md_canon$locality <- md$locality
md_canon$restoration4 <- as.character(md$restoration4)
md_canon$depth_cm <- md$depth_cm
write_csv(md_canon, file.path(TABLE_DIR, "001_metadata_canon_51.csv"))

gates <- tibble(
  gate = c(
    "required_columns", "unique_sample_ids", "no_invalid_design_rows",
    "all_profiles_complete_5_20_40", "expected_sample_count",
    "expected_profile_count", "model_matrix_full_rank",
    "conserved_reference_each_locality_depth", "loo_conserved_reference_each_locality_depth",
    "canonical_samples_present_in_phyloseq", "phyloseq_pruned_sample_count"
  ),
  pass = c(
    !length(missing_required), nrow(duplicate_samples) == 0L, nrow(invalid_rows) == 0L,
    nrow(profile_audit) > 0L && all(profile_audit$valid_complete_profile),
    n_distinct(md$sample_id) == EXPECTED_SAMPLES,
    n_distinct(md$profile_id) == EXPECTED_PROFILES,
    rank_x == ncol(X),
    reference_audit %>% filter(restoration4 == "Conserved") %>% pull(reference_available) %>% all(),
    reference_audit %>% filter(restoration4 == "Conserved") %>% pull(loo_reference_possible) %>% all(),
    all(canonical_ids %in% ps_source_ids),
    phyloseq::nsamples(ps_canon) == EXPECTED_SAMPLES
  ),
  severity = c("critical", "critical", "critical", "critical", "critical", "critical",
               "critical", "critical", "informative", "critical", "critical"),
  detail = c(
    paste(required, collapse = ";"),
    paste(nrow(duplicate_samples), "duplicated IDs"),
    paste(nrow(invalid_rows), "invalid rows"),
    paste(sum(!profile_audit$valid_complete_profile), "invalid profiles"),
    paste(n_distinct(md$sample_id), "observed;", EXPECTED_SAMPLES, "expected"),
    paste(n_distinct(md$profile_id), "observed;", EXPECTED_PROFILES, "expected"),
    paste("rank", rank_x, "of", ncol(X)),
    "Required for local-by-depth reference analyses",
    "Not required: if FALSE, exclude conserved anchor from inference",
    paste(sum(!canonical_ids %in% ps_source_ids), "canonical IDs absent from source phyloseq"),
    paste(phyloseq::nsamples(ps_canon), "samples retained in pruned phyloseq")
  )
)
write_csv(gates, file.path(TABLE_DIR, "001_compuertas_diseno.csv"))
write_csv(md, file.path(TABLE_DIR, "001_metadata_diseno_normalizada.csv"))

summary_tbl <- tibble(
  metric = c(
    "n_rows", "n_unique_samples", "n_profiles", "n_localities", "n_stages",
    "n_depths", "n_empty_locality_stage_cells", "n_single_profile_cells",
    "n_invalid_profiles", "model_matrix_rank", "model_matrix_columns"
  ),
  value = c(
    nrow(md), n_distinct(md$sample_id), n_distinct(md$profile_id),
    n_distinct(md$locality), n_distinct(md$restoration4, na.rm = TRUE),
    n_distinct(md$depth_cm, na.rm = TRUE), sum(design_counts$empty_cell),
    sum(design_counts$single_profile_cell), sum(!profile_audit$valid_complete_profile),
    rank_x, ncol(X)
  )
)
write_csv(summary_tbl, file.path(TABLE_DIR, "001_resumen_diseno.csv"))

critical_fail <- gates %>% filter(severity == "critical", !pass)
cat("\nDesign summary\n")
print(summary_tbl)
cat("\nGates\n")
print(gates)

if (nrow(critical_fail)) {
  stop(
    "AUDITORIA FALLIDA. No ejecutar inferencia. Compuertas criticas: ",
    paste(critical_fail$gate, collapse = ", "), call. = FALSE
  )
}

sample_data_canon <- as.data.frame(md_canon, stringsAsFactors = FALSE)
rownames(sample_data_canon) <- sample_data_canon$sample_id
sample_data_canon <- sample_data_canon[
  phyloseq::sample_names(ps_canon), , drop = FALSE
]
phyloseq::sample_data(ps_canon) <- phyloseq::sample_data(sample_data_canon)
saveRDS(ps_canon, file.path(RDS_DIR, "001_phyloseq_canon_51.rds"))
writeLines(OUT_DIR, LATEST_FILE)

writeLines(
  c(
    "AUDITORIA APROBADA PARA CONTINUAR.",
    "Universo analitico unico: 51 muestras en 17 perfiles.",
    "El phyloseq fuente fue podado por sample_id antes de cualquier analisis.",
    "Unidad experimental primaria: profile_id.",
    "Las profundidades 5/20/40 cm son submuestras repetidas del mismo perfil.",
    "Celdas con un solo perfil no sostienen contrastes locales independientes.",
    "Una unica referencia Conserved puede anclar muestras no conservadas, pero debe excluirse de su propia inferencia."
  ),
  file.path(OUT_DIR, "001_LEEME_resultado_auditoria.txt")
)

cat("\nAUDITORIA APROBADA\n")
cat("Finished UTC:", format(Sys.time(), tz = "UTC"), "\n")
