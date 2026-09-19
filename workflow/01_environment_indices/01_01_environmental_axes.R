#!/usr/bin/env Rscript
# ============================================================
# 01_01_environmental_axes.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# PCA of environmental blocks and the 19-variable ECI on the analytical cohort.
# Inputs (source expressions; complete list in docs/contracts/01_01_environmental_axes.json):
#   meta_raw <- readr::read_csv(META_CSV, show_col_types = FALSE)
#   cohort_raw <- readr::read_csv(COHORT_CSV, show_col_types = FALSE)
#   cap <- readr::read_csv(CAP_SCORES, show_col_types = FALSE)
#   readr::read_csv(global_var, show_col_types = FALSE)
# Outputs (source expressions; complete list in contract):
#   readr::write_csv(scores, file.path(out_tables, paste0(prefix, "_scores.csv")))
#   readr::write_csv(loadings, file.path(out_tables, paste0(prefix, "_loadings.csv")))
#   readr::write_csv(variance, file.path(out_tables, paste0(prefix, "_variance_explained.csv")))
#   ggsave(
#   writeLines(OUT_DIR, LATEST_ATTEMPT_FILE)
#   readr::write_csv(profile_audit, file.path(DIR_TABLES, "101_auditoria_perfiles_canon_51.csv"))
#   readr::write_csv(
#   readr::write_csv(audit, file.path(DIR_TABLES, "101_env_variable_audit.csv"))
#   readr::write_csv(cor_df, file.path(DIR_TABLES, "101_env_correlation_matrix_long.csv"))
#   readr::write_csv(cor_wide, file.path(DIR_TABLES, "101_env_correlation_matrix.csv"))
# Algorithmic provenance:
# Median imputation, mean/SD scaling, block PCA and 19-variable global ECI PCA.
#   Study-specific environmental calibration and data integration; see docs/METHODS.md.
# Source SHA-256: 24b8b1fc3d5fc2f85837e201fd1422b760bdf5b8ce69f71d39534ffedc868a82
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


# ==============================================================================
# 101_construir_ejes_ambientales_por_bloque.R
#
# Objetivo:
#   Auditar, ordenar y sintetizar variables ambientales/espaciales asociadas
#   exclusivamente al canon de 51 muestras y 17 perfiles completos.
#
#   Este script NO analiza la comunidad completa.
#   Prepara los gradientes ambientales para análisis posteriores de:
#
#     bimodal_CLR ~ restoration4 + depth_cm + env_axes + Condition(locality)
#
# Inputs principales:
#   1) Metadata con sample_id, restoration4, depth_cm, locality y variables env.
#   2) Cohorte canónica de 51 muestras; es la única fuente autorizada para
#      definir el universo analítico.
#   3) Opcional: CAP scores, usados solo como comprobación de IDs. Nunca
#      redefinen la cohorte.
#
# Contrato ECI-51:
#   - descarta cualquier ECI o eje ambiental heredado de la metadata de entrada;
#   - recalcula el PCA global después de restringir a las 51 muestras;
#   - conserva la definición histórica del ECI: 19 variables, sin shdi ni
#     depth_cm, con exclusión de variables que superen 20% de faltantes;
#   - exporta esos scores nuevos como ECI_PC1, ECI_PC2 y ECI_PC3;
#   - conserva loadings y procedencia para interpretar el signo arbitrario de
#     cada componente.
#
# Outputs:
#   101_env_variable_audit.csv
#   101_env_correlation_matrix.csv
#   101_env_correlation_heatmap.png
#   101_env_vif_global.csv
#   101_env_pca_global_scores.csv
#   101_env_pca_global_loadings.csv
#   101_env_pca_global_PC1_correlations.csv
#   101_env_within_profile_variation.csv
#   101_eci_variable_contract.csv
#   101_env_pca_by_block_scores.csv
#   101_env_pca_by_block_loadings.csv
#   101_axis_provenance.csv
#   101_eci_canon_gate.csv
#   101_metadata_with_env_axes.csv
#
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(purrr)
})

# ==============================================================================
# 0. Helpers
# ==============================================================================

timestamp_now <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
  cat("\n")
  flush.console()
}

safe_mkdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[[hit + 1]]
}

stop_if_missing_file <- function(path, label = "file") {
  if (is.null(path) || is.na(path) || !file.exists(path)) {
    stop(sprintf("Missing %s: %s", label, path), call. = FALSE)
  }
}

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

normalize_depth <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\.0$", "", x)
  x
}

clean_colnames <- function(x) {
  x <- make.names(x, unique = TRUE)
  x
}

median_impute <- function(df) {
  out <- df
  for (nm in names(out)) {
    x <- safe_num(out[[nm]])
    med <- median(x, na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    x[!is.finite(x)] <- med
    out[[nm]] <- x
  }
  out
}

scale_numeric_df <- function(df) {
  out <- as.data.frame(scale(df))
  out[] <- lapply(out, function(x) {
    x[!is.finite(x)] <- 0
    x
  })
  out
}

safe_cor <- function(df) {
  mat <- suppressWarnings(cor(df, use = "pairwise.complete.obs", method = "spearman"))
  mat[!is.finite(mat)] <- NA_real_
  mat
}

calc_vif <- function(df) {
  df <- as.data.frame(df)
  df <- df[, sapply(df, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]

  if (ncol(df) < 2) {
    return(tibble(variable = names(df), vif = NA_real_, r2 = NA_real_, note = "too_few_variables"))
  }

  map_dfr(names(df), function(v) {
    others <- setdiff(names(df), v)

    dat <- df[, c(v, others), drop = FALSE]
    dat <- dat[stats::complete.cases(dat), , drop = FALSE]

    if (nrow(dat) < length(others) + 3) {
      return(tibble(variable = v, vif = NA_real_, r2 = NA_real_, note = "too_few_complete_cases"))
    }

    form <- as.formula(paste(v, "~", paste(others, collapse = " + ")))

    fit <- tryCatch(lm(form, data = dat), error = function(e) NULL)

    if (is.null(fit)) {
      return(tibble(variable = v, vif = NA_real_, r2 = NA_real_, note = "lm_failed"))
    }

    r2 <- summary(fit)$r.squared
    vif <- ifelse(is.finite(r2) && r2 < 1, 1 / (1 - r2), Inf)

    tibble(variable = v, vif = vif, r2 = r2, note = "computed_not_acceptability")
  }) %>%
    arrange(desc(vif))
}

run_pca <- function(df_scaled, sample_ids, prefix, out_tables, out_plots, meta_plot = NULL) {
  if (ncol(df_scaled) < 2 || nrow(df_scaled) < 3) {
    msg("Skipping PCA for ", prefix, ": insufficient dimensions")
    return(NULL)
  }

  pca <- prcomp(df_scaled, center = FALSE, scale. = FALSE)

  eig <- pca$sdev^2
  var_exp <- eig / sum(eig)

  scores <- as_tibble(pca$x, rownames = NA) %>%
    mutate(sample_id = sample_ids, .before = 1)

  loadings <- as_tibble(pca$rotation, rownames = "variable") %>%
    mutate(variable = as.character(variable))

  variance <- tibble(
    axis = paste0("PC", seq_along(var_exp)),
    eigenvalue = eig,
    variance_explained = var_exp,
    cumulative_variance = cumsum(var_exp)
  )

  readr::write_csv(scores, file.path(out_tables, paste0(prefix, "_scores.csv")))
  readr::write_csv(loadings, file.path(out_tables, paste0(prefix, "_loadings.csv")))
  readr::write_csv(variance, file.path(out_tables, paste0(prefix, "_variance_explained.csv")))

  if (!is.null(meta_plot) && all(c("sample_id", "restoration4", "depth_cm") %in% names(meta_plot))) {
    plot_df <- scores %>%
      left_join(meta_plot, by = "sample_id")

    if (all(c("PC1", "PC2") %in% names(plot_df))) {
      p <- ggplot(plot_df, aes(x = PC1, y = PC2, shape = depth_cm)) +
        geom_point(size = 3, alpha = 0.85) +
        labs(
          title = paste0(prefix, " PCA"),
          subtitle = paste0(
            "PC1 = ", round(100 * var_exp[1], 1), "%; ",
            "PC2 = ", round(100 * var_exp[2], 1), "%"
          ),
          x = paste0("PC1 (", round(100 * var_exp[1], 1), "%)"),
          y = paste0("PC2 (", round(100 * var_exp[2], 1), "%)")
        ) +
        theme_bw(base_size = 12) +
        theme(panel.grid.minor = element_blank())

      ggsave(
        file.path(out_plots, paste0(prefix, "_PCA_PC1_PC2.png")),
        p,
        width = 7.5,
        height = 5.5,
        dpi = 300
      )
    }
  }

  list(pca = pca, scores = scores, loadings = loadings, variance = variance)
}

# ==============================================================================
# 1. Parámetros
# ==============================================================================

RESULTS_DIR <- "/home/fjbalvino/Tipping_points/resultados_finales"

DEFAULT_OUT12 <- NA_character_

DEFAULT_META <- "/data/Ciencia-Frontera/Metadata/Metricas_mangrove_metadata_clean.csv"
DEFAULT_COHORT <- "/home/fjbalvino/Tipping_points/Results/044_build_protein_filtered_matrix_HI_20260529_000843/tables/044_metadata_aligned_51samples.csv"

OUT12 <- get_arg("--out12", DEFAULT_OUT12)
META_CSV <- get_arg("--meta_csv", DEFAULT_META)
COHORT_CSV <- get_arg("--cohort_csv", DEFAULT_COHORT)
OUT_ROOT <- get_arg("--out_root", RESULTS_DIR)
EXPECTED_N <- suppressWarnings(as.integer(get_arg("--expected_n", "51")))
EXPECTED_PROFILES <- suppressWarnings(as.integer(get_arg("--expected_profiles", "17")))
DROP_IF_NA_GT <- suppressWarnings(as.numeric(get_arg("--drop_if_na_gt", "0.20")))
if (EXPECTED_N != 51L || EXPECTED_PROFILES != 17L) {
  stop("El canon esta congelado en 51 muestras y 17 perfiles; no se permite redefinirlo.", call. = FALSE)
}
if (!is.finite(DROP_IF_NA_GT) || DROP_IF_NA_GT < 0 || DROP_IF_NA_GT > 1) {
  stop("--drop_if_na_gt debe estar entre 0 y 1.", call. = FALSE)
}

SCRIPT_BASENAME <- "101_construir_ejes_ambientales_por_bloque"

OUT_DIR <- file.path(
  OUT_ROOT,
  paste0(SCRIPT_BASENAME, "_", timestamp_now())
)

DIR_TABLES <- file.path(OUT_DIR, "tables")
DIR_PLOTS <- file.path(OUT_DIR, "plots")
DIR_LOGS <- file.path(OUT_DIR, "logs")

safe_mkdir(OUT_DIR)
safe_mkdir(DIR_TABLES)
safe_mkdir(DIR_PLOTS)
safe_mkdir(DIR_LOGS)

LOG_FILE <- file.path(OUT_DIR, "log.txt")
sink(LOG_FILE, split = TRUE)
on.exit(sink(), add = TRUE)

LATEST_FILE <- file.path(OUT_ROOT, paste0("LATEST_", SCRIPT_BASENAME, ".txt"))
LATEST_ATTEMPT_FILE <- file.path(
  OUT_ROOT,
  paste0("LATEST_ATTEMPT_", SCRIPT_BASENAME, ".txt")
)
writeLines(OUT_DIR, LATEST_ATTEMPT_FILE)

msg("Starting script 101")
msg("META_CSV: ", META_CSV)
msg("COHORT_CSV: ", COHORT_CSV)
msg("OUT12: ", OUT12)
msg("DROP_IF_NA_GT: ", DROP_IF_NA_GT)
msg("OUT_DIR: ", OUT_DIR)

stop_if_missing_file(META_CSV, "metadata csv")
stop_if_missing_file(COHORT_CSV, "canonical cohort csv")

STAGE_LEVELS <- c(
  "Degraded",
  "Early restoration",
  "Intermediate restoration",
  "Advanced restoration",
  "Conserved"
)

DEPTH_LEVELS <- c("5", "20", "40")

# ==============================================================================
# 2. Leer metadata ambiental y restringir por el canon externo de 51 muestras
# ==============================================================================

meta_raw <- readr::read_csv(META_CSV, show_col_types = FALSE)
names(meta_raw) <- clean_colnames(names(meta_raw))

cohort_raw <- readr::read_csv(COHORT_CSV, show_col_types = FALSE)
names(cohort_raw) <- clean_colnames(names(cohort_raw))

sample_col <- intersect(
  c("sample_id", "SampleID", "sample", "id_metagenomics", "run_accession"),
  names(meta_raw)
)[1]

if (is.na(sample_col)) {
  stop("No pude identificar columna de muestra en metadata. Esperaba sample_id/id_metagenomics/run_accession.", call. = FALSE)
}

cohort_sample_col <- intersect(
  c("sample_id", "SampleID", "sample", "id_metagenomics", "run_accession"),
  names(cohort_raw)
)[1]
if (is.na(cohort_sample_col)) {
  stop("No pude identificar sample_id en la cohorte canonica.", call. = FALSE)
}
required_cohort <- c("restoration4", "depth_cm", "locality", "lat_block")
missing_cohort <- setdiff(required_cohort, names(cohort_raw))
if (length(missing_cohort)) {
  stop("Faltan columnas en cohorte canonica: ", paste(missing_cohort, collapse = ", "), call. = FALSE)
}

cohort <- cohort_raw %>%
  mutate(
    sample_id = as.character(.data[[cohort_sample_col]]),
    restoration4 = as.character(restoration4),
    depth_cm = normalize_depth(depth_cm),
    locality = as.character(locality),
    lat_block = as.character(lat_block)
  ) %>%
  mutate(
    restoration4 = factor(restoration4, levels = STAGE_LEVELS, ordered = FALSE),
    depth_cm = factor(depth_cm, levels = DEPTH_LEVELS, ordered = TRUE)
  ) %>%
  select(sample_id, restoration4, depth_cm, locality, lat_block)

if (
  nrow(cohort) != EXPECTED_N || n_distinct(cohort$sample_id) != EXPECTED_N ||
  anyNA(cohort)
) {
  stop("La cohorte no contiene exactamente 51 muestras completas y unicas.", call. = FALSE)
}

profile_audit <- cohort %>%
  group_by(lat_block) %>%
  summarise(
    n = n(),
    n_stage = n_distinct(restoration4),
    n_locality = n_distinct(locality),
    depths = paste(sort(as.character(depth_cm)), collapse = "|"),
    valid = n == 3L && n_stage == 1L && n_locality == 1L &&
      setequal(as.character(depth_cm), DEPTH_LEVELS),
    .groups = "drop"
  )
readr::write_csv(profile_audit, file.path(DIR_TABLES, "101_auditoria_perfiles_canon_51.csv"))
if (nrow(profile_audit) != EXPECTED_PROFILES || !all(profile_audit$valid)) {
  stop("La cohorte no contiene exactamente 17 perfiles completos 5/20/40.", call. = FALSE)
}

# Los ejes derivados presentes en META_CSV pueden provenir del universo previo
# de 57 muestras. Se eliminan antes del join para impedir su propagación.
derived_axis_aliases <- c(
  "ECI_PC1", "ECI_PC2", "ECI_PC3",
  "ECI.PC1", "ECI.PC2", "ECI.PC3",
  "env_global_PC1", "env_global_PC2", "env_global_PC3",
  "vegetation_landscape_PC1", "water_inundation_PC1",
  "inundation_PC1", "moisture_stress_PC1",
  "physicochemical_PC1", "nutrients_redox_PC1"
)
derived_axis_cols_in_input <- names(meta_raw)[
  tolower(names(meta_raw)) %in% tolower(derived_axis_aliases)
]

readr::write_csv(
  tibble(
    column = derived_axis_cols_in_input,
    action = rep(
      "dropped_before_refitting_on_canonical_51",
      length(derived_axis_cols_in_input)
    )
  ),
  file.path(DIR_TABLES, "101_inherited_axis_columns_dropped.csv")
)

meta_environment <- meta_raw %>%
  mutate(sample_id = as.character(.data[[sample_col]])) %>%
  select(-any_of(c(
    "restoration4", "stage", "condition", "depth_cm", "locality",
    "lat_block", "profile_id", derived_axis_cols_in_input
  )))
if (anyDuplicated(meta_environment$sample_id)) {
  stop("La metadata ambiental contiene sample_id duplicados.", call. = FALSE)
}

missing_environment_ids <- setdiff(cohort$sample_id, meta_environment$sample_id)
readr::write_csv(
  tibble(sample_id = missing_environment_ids),
  file.path(DIR_TABLES, "101_ids_canon_sin_metadata_ambiental.csv")
)
if (length(missing_environment_ids)) {
  stop("Faltan muestras canonicas en la metadata ambiental.", call. = FALSE)
}

meta <- cohort %>%
  left_join(meta_environment, by = "sample_id")

CAP_SCORES <- first_existing_file(c(
  file.path(OUT12, "tables", "302_CAP_scores_source_data.csv"),
  file.path(OUT12, "tables", "012_CAP_scores_restoration4_Condition_locality.csv"),
  file.path(OUT12, "012_CAP_scores_restoration4_Condition_locality.csv")
))

if (!is.na(CAP_SCORES)) {
  cap <- readr::read_csv(CAP_SCORES, show_col_types = FALSE)
  names(cap) <- clean_colnames(names(cap))

  cap_sample_col <- intersect(
    c("sample_id", "SampleID", "sample", "id_metagenomics", "run_accession"),
    names(cap)
  )[1]

  if (!is.na(cap_sample_col)) {
    cap_samples <- unique(as.character(cap[[cap_sample_col]]))
    if (!setequal(cap_samples, cohort$sample_id)) {
      stop("OUT12 no coincide exactamente con el canon de 51 muestras; no se permite redefinir la cohorte por un outcome.", call. = FALSE)
    }
    msg("CAP sample IDs agree exactly with the frozen 51-sample cohort.")
  } else {
    msg("CAP score file found, but sample column not recognized. Not restricting.")
  }
} else {
  msg("CAP score file not found. Using the frozen 51-sample cohort.")
}

# ==============================================================================
# 3. Definir variables ambientales y bloques
# ==============================================================================

env_blocks <- list(
  vegetation_landscape = c(
    "shdi",
    "ndvi_p50",
    "ndvi_sd"
  ),
  water_inundation = c(
    "ndwi_p50",
    "ndwi_sd",
    "mndwi_p50",
    "mndwi_sd",
    "water_prop_mndwi"
  ),
  moisture_stress = c(
    "ndmi_p50",
    "ndmi_sd"
  ),
  physicochemical = c(
    "temperature_c",
    "salinity_ups",
    "ph",
    "redox_shallow_mv"
  ),
  nutrients_redox = c(
    "n_no2_umol_l",
    "n_no3_umol_l",
    "n_nh3_mg_l_campo",
    "n_nh4plus_umol_l",
    "p_po4_3_umol_l",
    "s_2_umol"
  )
)

all_env_vars_requested <- unique(unlist(env_blocks))

# El ECI global conserva el contrato del script histórico 29: shdi puede
# participar en el bloque vegetation_landscape, pero no en el PCA global.
# depth_cm tampoco entra en ningún PCA ambiental.
eci_global_vars_requested <- c(
  "ndvi_p50", "ndvi_sd",
  "ndwi_p50", "ndwi_sd", "mndwi_p50", "mndwi_sd", "water_prop_mndwi",
  "ndmi_p50", "ndmi_sd",
  "temperature_c", "salinity_ups", "ph", "redox_shallow_mv",
  "n_no2_umol_l", "n_no3_umol_l", "n_nh3_mg_l_campo",
  "n_nh4plus_umol_l", "p_po4_3_umol_l", "s_2_umol"
)

eci_variable_contract <- tibble(
  variable = unique(c(eci_global_vars_requested, "shdi", "depth_cm")),
  role = c(
    rep("included_in_global_ECI", length(eci_global_vars_requested)),
    "excluded_from_global_ECI_block_only",
    "excluded_from_all_environmental_PCA"
  )
)
readr::write_csv(
  eci_variable_contract,
  file.path(DIR_TABLES, "101_eci_variable_contract.csv")
)

env_vars_present <- intersect(all_env_vars_requested, names(meta))
env_vars_missing <- setdiff(all_env_vars_requested, names(meta))

msg("Environmental variables requested: ", length(all_env_vars_requested))
msg("Environmental variables present: ", length(env_vars_present))
msg("Environmental variables missing: ", length(env_vars_missing))

if (length(env_vars_present) < 2) {
  stop(
    "Muy pocas variables ambientales presentes. Revisa nombres de columnas en metadata.\nFaltantes: ",
    paste(env_vars_missing, collapse = ", "),
    call. = FALSE
  )
}

# ==============================================================================
# 4. Auditoría de variables ambientales
# ==============================================================================

audit <- map_dfr(all_env_vars_requested, function(v) {
  if (!v %in% names(meta)) {
    return(tibble(
      variable = v,
      present = FALSE,
      n = nrow(meta),
      n_missing = NA_integer_,
      prop_missing = NA_real_,
      n_unique = NA_integer_,
      mean = NA_real_,
      sd = NA_real_,
      min = NA_real_,
      median = NA_real_,
      max = NA_real_
    ))
  }

  x <- safe_num(meta[[v]])

  tibble(
    variable = v,
    present = TRUE,
    n = length(x),
    n_missing = sum(!is.finite(x)),
    prop_missing = mean(!is.finite(x)),
    n_unique = length(unique(x[is.finite(x)])),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}) %>%
  mutate(
    block = case_when(
      variable %in% env_blocks$vegetation_landscape ~ "vegetation_landscape",
      variable %in% env_blocks$water_inundation ~ "water_inundation",
      variable %in% env_blocks$moisture_stress ~ "moisture_stress",
      variable %in% env_blocks$physicochemical ~ "physicochemical",
      variable %in% env_blocks$nutrients_redox ~ "nutrients_redox",
      TRUE ~ "unassigned"
    ),
    included_in_global_ECI = variable %in% eci_global_vars_requested,
    passes_missingness_threshold = present & is.finite(prop_missing) &
      prop_missing <= DROP_IF_NA_GT
  ) %>%
  select(block, everything())

readr::write_csv(audit, file.path(DIR_TABLES, "101_env_variable_audit.csv"))

msg("Environmental variable audit:")
print(audit)

# ==============================================================================
# 5. Matriz ambiental limpia
# ==============================================================================

env_df_raw <- meta %>%
  select(sample_id, restoration4, depth_cm, locality, all_of(env_vars_present))

env_num <- env_df_raw %>%
  select(all_of(env_vars_present)) %>%
  mutate(across(everything(), safe_num))

# Remover variables con más de 20% de faltantes, completamente vacías o
# constantes. El umbral puede cambiarse solo de forma explícita por argumento.
keep_vars <- names(env_num)[sapply(env_num, function(x) {
  prop_missing <- mean(!is.finite(x))
  prop_missing <= DROP_IF_NA_GT &&
    sum(is.finite(x)) >= 5 &&
    is.finite(sd(x, na.rm = TRUE)) &&
    sd(x, na.rm = TRUE) > 0
})]

env_num <- env_num %>% select(all_of(keep_vars))

msg("Environmental variables retained after missing/variance filter: ", ncol(env_num))
msg("Retained variables: ", paste(names(env_num), collapse = ", "))

readr::write_csv(
  tibble(variable = names(env_num)),
  file.path(DIR_TABLES, "101_env_variables_retained.csv")
)

eci_global_vars_retained <- intersect(eci_global_vars_requested, names(env_num))
eci_global_vars_not_retained <- setdiff(
  eci_global_vars_requested,
  eci_global_vars_retained
)

if (length(eci_global_vars_not_retained)) {
  stop(
    "No se puede reconstruir el ECI histórico de 19 variables en el canon 51. ",
    "Variables ausentes o filtradas: ",
    paste(eci_global_vars_not_retained, collapse = ", "),
    call. = FALSE
  )
}

env_global_scaled <- env_num %>%
  select(all_of(eci_global_vars_requested)) %>%
  median_impute() %>%
  scale_numeric_df()

# Auditoría de variación intraperfil con el lat_block congelado en META_51.
# No se reconstruyen perfiles redondeando latitude.
within_profile_source <- bind_cols(
  meta %>% select(sample_id, lat_block),
  env_num
)

within_profile_variation <- map_dfr(names(env_num), function(v) {
  by_profile <- within_profile_source %>%
    group_by(lat_block) %>%
    summarise(
      n_finite = sum(is.finite(.data[[v]])),
      within_sd = ifelse(
        n_finite >= 2L,
        sd(.data[[v]][is.finite(.data[[v]])]),
        NA_real_
      ),
      .groups = "drop"
    )

  n_estimable <- sum(is.finite(by_profile$within_sd))
  n_varying <- sum(by_profile$within_sd > 0, na.rm = TRUE)

  tibble(
    variable = v,
    n_profiles_total = nrow(by_profile),
    n_profiles_estimable = n_estimable,
    n_profiles_sd_gt0 = n_varying,
    frac_profiles_sd_gt0 = ifelse(
      n_estimable > 0,
      n_varying / n_estimable,
      NA_real_
    )
  )
}) %>%
  arrange(desc(frac_profiles_sd_gt0), variable)

readr::write_csv(
  within_profile_variation,
  file.path(DIR_TABLES, "101_env_within_profile_variation.csv")
)

# ==============================================================================
# 6. Correlación ambiental
# ==============================================================================

cor_mat <- safe_cor(env_num)

cor_df <- as.data.frame(cor_mat) %>%
  rownames_to_column("variable_1") %>%
  pivot_longer(-variable_1, names_to = "variable_2", values_to = "spearman_rho")

readr::write_csv(cor_df, file.path(DIR_TABLES, "101_env_correlation_matrix_long.csv"))

cor_wide <- as.data.frame(cor_mat) %>%
  rownames_to_column("variable")

readr::write_csv(cor_wide, file.path(DIR_TABLES, "101_env_correlation_matrix.csv"))

p_cor <- cor_df %>%
  mutate(
    variable_1 = factor(variable_1, levels = names(env_num)),
    variable_2 = factor(variable_2, levels = names(env_num))
  ) %>%
  ggplot(aes(x = variable_1, y = variable_2, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    limits = c(-1, 1),
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "grey90"
  ) +
  labs(
    title = "Environmental correlation matrix",
    subtitle = "Spearman correlation among environmental variables",
    x = NULL,
    y = NULL,
    fill = "rho"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  )

ggsave(
  file.path(DIR_PLOTS, "101_env_correlation_heatmap.png"),
  p_cor,
  width = 9,
  height = 8,
  dpi = 300
)

# ==============================================================================
# 7. VIF global
# ==============================================================================

vif_global <- calc_vif(env_global_scaled)

readr::write_csv(vif_global, file.path(DIR_TABLES, "101_env_vif_global.csv"))

msg("Global VIF:")
print(vif_global)

# ==============================================================================
# 8. PCA global ambiental
# ==============================================================================

meta_plot <- env_df_raw %>%
  select(sample_id, restoration4, depth_cm, locality)

pca_global <- run_pca(
  df_scaled = env_global_scaled,
  sample_ids = env_df_raw$sample_id,
  prefix = "101_env_pca_global",
  out_tables = DIR_TABLES,
  out_plots = DIR_PLOTS,
  meta_plot = meta_plot
)

if (!is.null(pca_global)) {
  pc1_correlations <- map_dfr(names(env_global_scaled), function(v) {
    tibble(
      variable = v,
      pearson_r = suppressWarnings(cor(
        env_global_scaled[[v]], pca_global$scores$PC1,
        use = "complete.obs", method = "pearson"
      )),
      spearman_rho = suppressWarnings(cor(
        env_global_scaled[[v]], pca_global$scores$PC1,
        use = "complete.obs", method = "spearman"
      ))
    )
  }) %>%
    mutate(abs_pearson_r = abs(pearson_r)) %>%
    arrange(desc(abs_pearson_r))

  readr::write_csv(
    pc1_correlations,
    file.path(DIR_TABLES, "101_env_pca_global_PC1_correlations.csv")
  )
}

# ==============================================================================
# 9. PCA por bloques ambientales
# ==============================================================================

block_scores_list <- list()
block_loadings_list <- list()
block_variance_list <- list()

for (block_name in names(env_blocks)) {
  vars_block <- intersect(env_blocks[[block_name]], names(env_num))

  if (length(vars_block) < 2) {
    msg("Skipping block ", block_name, ": fewer than 2 retained variables")
    next
  }

  block_raw <- env_num %>% select(all_of(vars_block))
  block_scaled <- block_raw %>%
    median_impute() %>%
    scale_numeric_df()

  pca_block <- run_pca(
    df_scaled = block_scaled,
    sample_ids = env_df_raw$sample_id,
    prefix = paste0("101_env_pca_block_", block_name),
    out_tables = DIR_TABLES,
    out_plots = DIR_PLOTS,
    meta_plot = meta_plot
  )

  if (!is.null(pca_block)) {
    block_scores <- pca_block$scores %>%
      select(sample_id, PC1, PC2) %>%
      rename(
        !!paste0(block_name, "_PC1") := PC1,
        !!paste0(block_name, "_PC2") := PC2
      )

    block_loadings <- pca_block$loadings %>%
      mutate(block = block_name, .before = 1)

    block_variance <- pca_block$variance %>%
      mutate(block = block_name, .before = 1)

    block_scores_list[[block_name]] <- block_scores
    block_loadings_list[[block_name]] <- block_loadings
    block_variance_list[[block_name]] <- block_variance
  }
}

block_scores_joined <- reduce(block_scores_list, full_join, by = "sample_id")
block_loadings_all <- bind_rows(block_loadings_list)
block_variance_all <- bind_rows(block_variance_list)

if (!is.null(block_scores_joined) && nrow(block_scores_joined) > 0) {
  readr::write_csv(
    block_scores_joined,
    file.path(DIR_TABLES, "101_env_pca_by_block_scores.csv")
  )
}

if (nrow(block_loadings_all) > 0) {
  readr::write_csv(
    block_loadings_all,
    file.path(DIR_TABLES, "101_env_pca_by_block_loadings.csv")
  )
}

if (nrow(block_variance_all) > 0) {
  readr::write_csv(
    block_variance_all,
    file.path(DIR_TABLES, "101_env_pca_by_block_variance_explained.csv")
  )
}

# ==============================================================================
# 10. Metadata final con ejes ambientales
# ==============================================================================

meta_env_axes <- meta %>%
  select(sample_id, restoration4, depth_cm, locality, everything())

if (!is.null(pca_global)) {
  global_scores_small <- pca_global$scores %>%
    select(sample_id, PC1, PC2, PC3) %>%
    rename(
      ECI_PC1 = PC1,
      ECI_PC2 = PC2,
      ECI_PC3 = PC3
    )

  meta_env_axes <- meta_env_axes %>%
    left_join(global_scores_small, by = "sample_id")
}

if (!is.null(block_scores_joined) && nrow(block_scores_joined) > 0) {
  meta_env_axes <- meta_env_axes %>%
    left_join(block_scores_joined, by = "sample_id")
}

canonical_eci_cols <- c("ECI_PC1", "ECI_PC2", "ECI_PC3")
missing_eci_cols <- setdiff(canonical_eci_cols, names(meta_env_axes))
forbidden_eci_aliases <- names(meta_env_axes)[
  tolower(names(meta_env_axes)) %in% c(
    "env_global_pc1", "env_global_pc2", "env_global_pc3"
  )
]

eci_finite_counts <- if (length(missing_eci_cols) == 0L) {
  vapply(
    canonical_eci_cols,
    function(cc) sum(is.finite(suppressWarnings(as.numeric(meta_env_axes[[cc]])))),
    integer(1)
  )
} else {
  setNames(rep(0L, length(canonical_eci_cols)), canonical_eci_cols)
}

eci_gate <- tibble(
  gate = c(
    "canonical_sample_count",
    "canonical_sample_order",
    "global_ECI_exact_19_variable_contract",
    "shdi_excluded_from_global_ECI",
    "depth_excluded_from_environmental_PCA",
    "canonical_ECI_columns_present",
    "canonical_ECI_columns_complete",
    "no_env_global_aliases",
    "PCA_fitted_after_canonical_subset"
  ),
  pass = c(
    nrow(meta_env_axes) == EXPECTED_N && n_distinct(meta_env_axes$sample_id) == EXPECTED_N,
    identical(as.character(meta_env_axes$sample_id), as.character(cohort$sample_id)),
    identical(names(env_global_scaled), eci_global_vars_requested),
    !"shdi" %in% names(env_global_scaled),
    !"depth_cm" %in% all_env_vars_requested && !"depth_cm" %in% names(env_global_scaled),
    length(missing_eci_cols) == 0L,
    length(missing_eci_cols) == 0L && all(eci_finite_counts == EXPECTED_N),
    length(forbidden_eci_aliases) == 0L,
    !is.null(pca_global) && nrow(pca_global$scores) == EXPECTED_N
  ),
  severity = "critical",
  detail = c(
    paste0(nrow(meta_env_axes), " rows; expected ", EXPECTED_N),
    "Output order must exactly match the frozen cohort order",
    paste(names(env_global_scaled), collapse = ";"),
    "shdi is reserved for the vegetation_landscape block",
    "depth_cm is a design factor, never a PCA input",
    ifelse(
      length(missing_eci_cols) == 0L,
      paste(canonical_eci_cols, collapse = ";"),
      paste("missing", paste(missing_eci_cols, collapse = ";"))
    ),
    paste(names(eci_finite_counts), eci_finite_counts, sep = "=", collapse = ";"),
    ifelse(
      length(forbidden_eci_aliases) == 0L,
      "0 ambiguous env_global aliases",
      paste(forbidden_eci_aliases, collapse = ";")
    ),
    "Global environmental PCA refitted after exact 51-sample cohort join"
  )
)

readr::write_csv(eci_gate, file.path(DIR_TABLES, "101_eci_canon_gate.csv"))
if (!all(eci_gate$pass)) {
  stop(
    "Fallo la compuerta de procedencia ECI-51; revise 101_eci_canon_gate.csv.",
    call. = FALSE
  )
}

axis_columns <- c(
  canonical_eci_cols,
  setdiff(names(block_scores_joined), "sample_id")
)
axis_provenance <- tibble(
  axis = axis_columns,
  source = if_else(
    axis %in% canonical_eci_cols,
    "global_PCA_refit_on_historical_19_variables_without_shdi_or_depth",
    "block_PCA_refit_on_raw_environmental_variables"
  ),
  analytical_universe = "canonical_51_samples_17_complete_profiles",
  n_samples_fitted = EXPECTED_N,
  n_profiles = EXPECTED_PROFILES,
  missingness_exclusion_threshold = DROP_IF_NA_GT,
  metadata_source = META_CSV,
  cohort_source = COHORT_CSV,
  inherited_axis_columns_dropped = paste(derived_axis_cols_in_input, collapse = ";"),
  sign_interpretation = "PCA sign is arbitrary; interpret scores jointly with loadings"
)
readr::write_csv(
  axis_provenance,
  file.path(DIR_TABLES, "101_axis_provenance.csv")
)

readr::write_csv(
  meta_env_axes,
  file.path(DIR_TABLES, "101_metadata_with_env_axes.csv")
)

# ==============================================================================
# 11. Figuras resumen de ejes ambientales por estadio/profundidad
# ==============================================================================

if ("ECI_PC1" %in% names(meta_env_axes)) {
  p_pc1_stage <- ggplot(meta_env_axes, aes(x = restoration4, y = ECI_PC1)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(aes(shape = factor(depth_cm)), width = 0.12, size = 2.4, alpha = 0.8) +
    labs(
      title = "Global environmental PC1 across restoration stages",
      subtitle = "Environmental axis derived from measured variables",
      x = "Restoration stage",
      y = "Environmental PC1",
      shape = "Depth"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1)
    )

  ggsave(
    file.path(DIR_PLOTS, "101_ECI_PC1_by_restoration4.png"),
    p_pc1_stage,
    width = 8,
    height = 5,
    dpi = 300
  )

  p_pc1_depth <- ggplot(meta_env_axes, aes(x = depth_cm, y = ECI_PC1)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, size = 2.4, alpha = 0.8) +
    labs(
      title = "Global environmental PC1 by sediment depth",
      subtitle = "Environmental axis derived from measured variables",
      x = "Depth",
      y = "Environmental PC1"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  ggsave(
    file.path(DIR_PLOTS, "101_ECI_PC1_by_depth_cm.png"),
    p_pc1_depth,
    width = 6.5,
    height = 5,
    dpi = 300
  )
}

# ==============================================================================
# 12. Resumen interpretativo automático
# ==============================================================================

global_var <- file.path(DIR_TABLES, "101_env_pca_global_variance_explained.csv")

global_variance_tbl <- if (file.exists(global_var)) {
  readr::read_csv(global_var, show_col_types = FALSE)
} else {
  tibble(axis = character(), variance_explained = numeric())
}

top_global_loadings <- if (!is.null(pca_global)) {
  pca_global$loadings %>%
    select(variable, PC1, PC2) %>%
    mutate(abs_PC1 = abs(PC1), abs_PC2 = abs(PC2)) %>%
    arrange(desc(abs_PC1)) %>%
    slice_head(n = 10)
} else {
  tibble()
}

readr::write_csv(
  top_global_loadings,
  file.path(DIR_TABLES, "101_env_pca_global_top_loadings_PC1.csv")
)

interpretation <- tibble(
  item = c(
    "n_samples",
    "n_env_variables_requested",
    "n_env_variables_present",
    "n_env_variables_retained",
    "n_global_ECI_variables",
    "global_ECI_excludes_shdi",
    "missingness_exclusion_threshold",
    "global_PC1_variance",
    "global_PC2_variance",
    "recommended_next_script"
  ),
  value = c(
    as.character(nrow(meta_env_axes)),
    as.character(length(all_env_vars_requested)),
    as.character(length(env_vars_present)),
    as.character(ncol(env_num)),
    as.character(ncol(env_global_scaled)),
    as.character(!"shdi" %in% names(env_global_scaled)),
    as.character(DROP_IF_NA_GT),
    ifelse(nrow(global_variance_tbl) >= 1, as.character(round(global_variance_tbl$variance_explained[1], 4)), NA_character_),
    ifelse(nrow(global_variance_tbl) >= 2, as.character(round(global_variance_tbl$variance_explained[2], 4)), NA_character_),
    "102_modelar_direccion_gradientes_ambientales.R"
  )
)

readr::write_csv(
  interpretation,
  file.path(DIR_TABLES, "101_interpretation_summary.csv")
)

msg("Interpretation summary:")
print(interpretation)

# ==============================================================================
# 13. Run info
# ==============================================================================

run_info <- tibble(
  key = c(
    "script",
    "timestamp",
    "meta_csv",
    "cohort_csv",
    "out12",
    "cap_scores",
    "out_dir",
    "n_samples",
    "n_profiles",
    "n_env_vars_present",
    "n_env_vars_retained",
    "n_global_eci_vars",
    "drop_if_na_gt",
    "eci_definition",
    "inherited_axis_columns_dropped",
    "latest_attempt_file",
    "latest_success_file"
  ),
  value = c(
    SCRIPT_BASENAME,
    timestamp_now(),
    META_CSV,
    COHORT_CSV,
    OUT12,
    CAP_SCORES,
    OUT_DIR,
    as.character(nrow(meta_env_axes)),
    as.character(nrow(profile_audit)),
    as.character(length(env_vars_present)),
    as.character(ncol(env_num)),
    as.character(ncol(env_global_scaled)),
    as.character(DROP_IF_NA_GT),
    "ECI_PC1-3 from historical 19-variable global PCA refitted on canonical 51 samples; shdi and depth_cm excluded",
    paste(derived_axis_cols_in_input, collapse = ";"),
    LATEST_ATTEMPT_FILE,
    LATEST_FILE
  )
)

readr::write_csv(run_info, file.path(DIR_TABLES, "101_run_info.csv"))
capture.output(sessionInfo(), file = file.path(DIR_LOGS, "101_sessionInfo.txt"))

# LATEST apunta únicamente a una corrida que alcanzó todas las compuertas.
writeLines(OUT_DIR, LATEST_FILE)

msg("Finished script 101 successfully")
msg("Output directory: ", OUT_DIR)
