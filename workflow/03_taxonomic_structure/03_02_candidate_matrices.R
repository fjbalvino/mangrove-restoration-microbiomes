#!/usr/bin/env Rscript
# ============================================================
# 03_02_candidate_matrices.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Exports current raw1031 counts and CLR; also retains the earlier residual84 branch.
# Inputs (source expressions; complete list in docs/contracts/03_02_candidate_matrices.json):
#   out <- trimws(readLines(path, n = 1L, warn = FALSE))
#   ps <- readRDS(PS_RDS)
#   master <- read_csv(MASTER_CSV, show_col_types = FALSE)
#   metrics <- read_csv(METRICS_CSV, show_col_types = FALSE) %>%
#   conservative <- read_csv(CONSERVATIVE_CSV, show_col_types = FALSE) %>%
#   taxmap <- read_csv(TAXMAP_CSV, show_col_types = FALSE) %>%
#   residualization_models <- read_csv(
#   run201_summary <- read_csv(RUN201_SUMMARY_CSV, show_col_types = FALSE)
# Outputs (source expressions; complete list in contract):
#   writeLines(OUT_DIR, LATEST_ATTEMPT_FILE)
#   write_csv(parameter_gate, file.path(TABLE_DIR, "202_compuerta_parametros.csv"))
#   write_csv(input_gate, file.path(TABLE_DIR, "202_compuerta_inputs.csv"))
#   write_csv(profile_gate, file.path(TABLE_DIR, "202_auditoria_canon_51.csv"))
#   write_csv(design_gate, file.path(TABLE_DIR, "202_compuerta_diseno.csv"))
#   write_csv(universe_gate, file.path(TABLE_DIR, "202_compuerta_universos_taxonomicos.csv"))
#   write_csv(membership, file.path(TABLE_DIR, "202_membresia_universos_84_12_3.csv"))
#   write_csv(membership, file.path(TABLE_DIR, "202_auditoria_mapeo_taxa.csv"))
#   write_csv(
#   write_csv(meta, file.path(TABLE_DIR, "202_metadata_alineada.csv"))
# Algorithmic provenance:
# Build candidate taxonomic matrices and CLR subcompositions from explicit screen membership.
#   McArdle & Anderson (2001), doi:10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2; Peres-Neto et al. (2006), doi:10.1890/0012-9658(2006)87[2614:VPOESD]2.0.CO;2.
# Source SHA-256: 6100bbbdcb38b8afe2a4fbcbd9650b2984b1f6eebffb4a5d7a29b76280196e29
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


# 202_construir_matriz_CLR_y_eje_bimodal.R
#
# Prepara dos universos taxonomicos congelados por 201:
#   1) 84 candidatos residuales no unimodales: universo focal principal.
#   2) 12 taxa conservadores con soporte bootstrap: sensibilidad e
#      interpretacion ecologica.
#
# Este script construye subcomposiciones CLR y ejes PCA descriptivos. No prueba
# restauracion ni reemplaza el orden narrativo del manuscrito: la comunidad
# completa se presenta primero en 301 y los candidatos se comparan despues.

suppressPackageStartupMessages({
  library(phyloseq)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Funciones
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (!length(i)) return(default)
  if (i[[1L]] == length(args)) {
    stop("Argumento sin valor: ", flag, call. = FALSE)
  }
  args[[i[[1L]] + 1L]]
}

as_int <- function(x, name) {
  out <- suppressWarnings(as.integer(x))
  if (length(out) != 1L || is.na(out)) {
    stop(name, " debe ser un entero", call. = FALSE)
  }
  out
}

as_num <- function(x, name) {
  out <- suppressWarnings(as.numeric(x))
  if (length(out) != 1L || !is.finite(out)) {
    stop(name, " debe ser numerico y finito", call. = FALSE)
  }
  out
}

as_logical_strict <- function(x) {
  z <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(z))
  out[z %in% c("TRUE", "T", "1")] <- TRUE
  out[z %in% c("FALSE", "F", "0")] <- FALSE
  out
}

stamp <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S", tz = "UTC")
}

read_latest <- function(root, id) {
  path <- file.path(root, paste0("LATEST_", id, ".txt"))
  if (!file.exists(path)) stop("No existe: ", path, call. = FALSE)
  out <- trimws(readLines(path, n = 1L, warn = FALSE))
  if (!length(out) || !nzchar(out) || !dir.exists(out)) {
    stop("LATEST invalido: ", path, call. = FALSE)
  }
  out
}

require_file <- function(path) {
  if (!file.exists(path)) stop("Falta input: ", path, call. = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

coalesce_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) return(rep(NA_character_, nrow(df)))
  out <- as.character(df[[hit[[1L]]]])
  if (length(hit) > 1L) {
    for (nm in hit[-1L]) {
      out <- dplyr::coalesce(out, as.character(df[[nm]]))
    }
  }
  out
}

normalize_depth <- function(x) {
  suppressWarnings(
    as.numeric(str_extract(as.character(x), "[0-9]+(?:\\.[0-9]+)?"))
  )
}

clr_rows <- function(m, pseudocount) {
  z <- log(m + pseudocount)
  sweep(z, 1L, rowMeans(z), "-")
}

make_oriented_pca <- function(m, meta, prefix) {
  keep <- vapply(
    seq_len(ncol(m)),
    function(j) all(is.finite(m[, j])) && stats::sd(m[, j]) > 0,
    logical(1)
  )
  names(keep) <- colnames(m)
  m_use <- m[, keep, drop = FALSE]
  if (ncol(m_use) < 2L) {
    stop("Muy pocos taxa variables para PCA ", prefix, call. = FALSE)
  }

  fit <- stats::prcomp(m_use, center = TRUE, scale. = FALSE)
  variance <- fit$sdev^2 / sum(fit$sdev^2)
  n_axes <- min(3L, ncol(fit$x))
  score_matrix <- fit$x[, seq_len(n_axes), drop = FALSE]

  med_degraded <- stats::median(
    score_matrix[meta$restoration4 == "Degraded", 1L],
    na.rm = TRUE
  )
  med_conserved <- stats::median(
    score_matrix[meta$restoration4 == "Conserved", 1L],
    na.rm = TRUE
  )
  flip <- if (
    is.finite(med_degraded) &&
      is.finite(med_conserved) &&
      med_conserved < med_degraded
  ) -1 else 1

  scores <- as_tibble(score_matrix)
  names(scores) <- paste0(prefix, "_PC", seq_len(n_axes), "_raw")
  scores[[paste0(prefix, "_PC1")]] <- flip * score_matrix[, 1L]

  loadings <- as_tibble(
    fit$rotation,
    rownames = "Original_Taxon_ID"
  )
  loadings$PC1 <- flip * loadings$PC1
  loadings <- loadings %>% arrange(desc(abs(PC1)))

  variance_table <- tibble(
    axis = paste0("PC", seq_along(variance)),
    variance_explained = variance,
    cumulative_variance = cumsum(variance)
  )

  list(
    fit = fit,
    matrix = m_use,
    scores = scores,
    loadings = loadings,
    variance = variance_table,
    flip = flip,
    n_input_taxa = ncol(m),
    n_variable_taxa = ncol(m_use)
  )
}

# -----------------------------------------------------------------------------
# Parametros y rutas
# -----------------------------------------------------------------------------

OUT_ROOT <- get_arg(
  "--out_root",
  "/home/fjbalvino/Tipping_points/resultados_finales"
)
RUN001 <- get_arg(
  "--run001",
  read_latest(OUT_ROOT, "001_auditar_estructura_metadata")
)
RUN201 <- get_arg(
  "--run201",
  read_latest(OUT_ROOT, "201_detectar_especies_bimodales_bootstrap")
)
PS_RDS <- get_arg(
  "--ps_rds",
  file.path(RUN001, "rds", "001_phyloseq_canon_51.rds")
)
MASTER_CSV <- get_arg(
  "--master_csv",
  file.path(RUN001, "tables", "001_metadata_canon_51.csv")
)
METRICS_CSV <- get_arg(
  "--metrics_csv",
  file.path(RUN201, "tables", "201_metricas_bimodalidad_todas_especies.csv")
)
CONSERVATIVE_CSV <- get_arg(
  "--conservative_csv",
  file.path(RUN201, "tables", "201_especies_bimodales_conservadoras.csv")
)
TAXMAP_CSV <- get_arg(
  "--taxmap_csv",
  file.path(RUN201, "tables", "201_mapa_taxonomia_especies.csv")
)
RESIDUALIZATION_CSV <- get_arg(
  "--residualization_csv",
  file.path(RUN201, "tables", "201_modelos_residualizacion.csv")
)
RUN201_SUMMARY_CSV <- get_arg(
  "--run201_summary_csv",
  file.path(RUN201, "tables", "201_parametros_y_resumen.csv")
)

PSEUDOCOUNT <- as_num(get_arg("--pseudocount", "1"), "pseudocount")
EXPECTED_CANDIDATES <- as_int(
  get_arg("--expected_candidates", "84"),
  "expected_candidates"
)
EXPECTED_CONSERVATIVE <- as_int(
  get_arg("--expected_conservative", "12"),
  "expected_conservative"
)
EXPECTED_ROBUST <- as_int(
  get_arg("--expected_robust", "3"),
  "expected_robust"
)

SCRIPT_ID <- "202_construir_matriz_CLR_y_eje_bimodal"
OUT_DIR <- file.path(OUT_ROOT, paste0(SCRIPT_ID, "_", stamp()))
TABLE_DIR <- file.path(OUT_DIR, "tables")
PLOT_DIR <- file.path(OUT_DIR, "plots")
LOG_DIR <- file.path(OUT_DIR, "logs")
LATEST_ATTEMPT_FILE <- file.path(
  OUT_ROOT,
  paste0("LATEST_ATTEMPT_", SCRIPT_ID, ".txt")
)
LATEST_FILE <- file.path(OUT_ROOT, paste0("LATEST_", SCRIPT_ID, ".txt"))

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
writeLines(OUT_DIR, LATEST_ATTEMPT_FILE)

log_con <- file(file.path(LOG_DIR, "202_log.txt"), "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0L) sink(type = "message")
  while (sink.number() > 0L) sink()
  close(log_con)
}, add = TRUE)

cat("===== 202 Prepare candidate and conservative CLR universes =====\n")
cat("Start UTC:", format(Sys.time(), tz = "UTC"), "\n")
cat("RUN001:", RUN001, "\n")
cat("RUN201:", RUN201, "\n")
cat("Output:", OUT_DIR, "\n")

parameter_gate <- tibble(
  gate = c(
    "pseudocount_positive",
    "expected_candidates_positive",
    "expected_conservative_positive",
    "expected_robust_nonnegative",
    "nested_expected_counts"
  ),
  pass = c(
    PSEUDOCOUNT > 0,
    EXPECTED_CANDIDATES > 0L,
    EXPECTED_CONSERVATIVE > 0L,
    EXPECTED_ROBUST >= 0L,
    EXPECTED_CANDIDATES >= EXPECTED_CONSERVATIVE &&
      EXPECTED_CONSERVATIVE >= EXPECTED_ROBUST
  ),
  detail = as.character(c(
    PSEUDOCOUNT,
    EXPECTED_CANDIDATES,
    EXPECTED_CONSERVATIVE,
    EXPECTED_ROBUST,
    paste(EXPECTED_CANDIDATES, EXPECTED_CONSERVATIVE, EXPECTED_ROBUST, sep = ">=")
  ))
)
write_csv(parameter_gate, file.path(TABLE_DIR, "202_compuerta_parametros.csv"))
if (!all(parameter_gate$pass)) {
  stop("Fallo la compuerta de parametros de 202", call. = FALSE)
}

PS_RDS <- require_file(PS_RDS)
MASTER_CSV <- require_file(MASTER_CSV)
METRICS_CSV <- require_file(METRICS_CSV)
CONSERVATIVE_CSV <- require_file(CONSERVATIVE_CSV)
TAXMAP_CSV <- require_file(TAXMAP_CSV)
RESIDUALIZATION_CSV <- require_file(RESIDUALIZATION_CSV)
RUN201_SUMMARY_CSV <- require_file(RUN201_SUMMARY_CSV)

# -----------------------------------------------------------------------------
# Phyloseq, conteos y metadata canonica
# -----------------------------------------------------------------------------

ps <- readRDS(PS_RDS)
otu <- as(phyloseq::otu_table(ps), "matrix")
if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)
storage.mode(otu) <- "numeric"

ps_meta_df <- as.data.frame(
  phyloseq::sample_data(ps),
  stringsAsFactors = FALSE
)
class(ps_meta_df) <- "data.frame"
ps_meta_row_ids <- as.character(rownames(ps_meta_df))
ps_meta_has_sample_id <- "sample_id" %in% names(ps_meta_df)
ps_meta_column_ids <- if (ps_meta_has_sample_id) {
  as.character(ps_meta_df[["sample_id"]])
} else {
  character(0)
}
ps_meta_id_agreement <- if (ps_meta_has_sample_id) {
  identical(ps_meta_column_ids, ps_meta_row_ids)
} else {
  TRUE
}
if (ps_meta_has_sample_id) ps_meta_df[["sample_id"]] <- NULL
ps_meta <- data.frame(
  sample_id = ps_meta_row_ids,
  ps_meta_df,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = NULL
)

master <- read_csv(MASTER_CSV, show_col_types = FALSE)
if (!"sample_id" %in% names(master)) {
  stop("master_csv requiere sample_id", call. = FALSE)
}
master <- master %>% mutate(sample_id = as.character(sample_id))

input_gate <- tibble(
  gate = c(
    "phyloseq_has_51_samples",
    "sample_data_existing_id_agrees_with_rownames",
    "otu_and_sample_data_ids_match",
    "master_has_51_unique_ids",
    "phyloseq_and_master_ids_match",
    "taxon_ids_unique",
    "otu_values_finite",
    "otu_values_nonnegative",
    "otu_values_integer_like",
    "all_samples_have_positive_library_size"
  ),
  pass = c(
    ncol(otu) == 51L && nrow(ps_meta) == 51L,
    isTRUE(ps_meta_id_agreement),
    !anyDuplicated(colnames(otu)) &&
      !anyDuplicated(ps_meta$sample_id) &&
      setequal(colnames(otu), ps_meta$sample_id),
    nrow(master) == 51L &&
      !anyNA(master$sample_id) &&
      !anyDuplicated(master$sample_id),
    setequal(ps_meta$sample_id, master$sample_id),
    !is.null(rownames(otu)) && !anyDuplicated(rownames(otu)),
    all(is.finite(otu)),
    all(is.finite(otu)) && all(otu >= 0),
    all(is.finite(otu)) && all(abs(otu - round(otu)) < 1e-8),
    all(is.finite(otu)) && all(colSums(otu) > 0)
  ),
  detail = c(
    paste0("otu=", ncol(otu), "; sample_data=", nrow(ps_meta)),
    if (ps_meta_has_sample_id) {
      paste0("column_present; exact_ordered_match=", ps_meta_id_agreement)
    } else {
      "column_absent; phyloseq rownames used"
    },
    paste0("otu_ids=", ncol(otu), "; metadata_ids=", nrow(ps_meta)),
    paste0("n=", nrow(master), "; duplicated=", anyDuplicated(master$sample_id)),
    paste0("shared=", length(intersect(ps_meta$sample_id, master$sample_id))),
    paste0("n_taxa=", nrow(otu)),
    paste0("nonfinite=", sum(!is.finite(otu))),
    paste0("negative=", sum(otu < 0, na.rm = TRUE)),
    "raw count matrix required",
    paste0("empty_samples=", sum(colSums(otu) <= 0))
  )
)
write_csv(input_gate, file.path(TABLE_DIR, "202_compuerta_inputs.csv"))
if (!all(input_gate$pass)) {
  stop("Fallo la compuerta de inputs de 202", call. = FALSE)
}

meta0 <- inner_join(
  ps_meta,
  master,
  by = "sample_id",
  suffix = c(".ps", ".master")
)

stages <- c(
  "Degraded",
  "Early restoration",
  "Intermediate restoration",
  "Advanced restoration",
  "Conserved"
)

meta <- meta0 %>%
  transmute(
    sample_id = as.character(sample_id),
    locality = str_squish(
      coalesce_col(meta0, c("locality.master", "locality", "locality.ps"))
    ),
    restoration4 = str_squish(
      coalesce_col(
        meta0,
        c(
          "restoration4.master",
          "restoration4",
          "restoration4.ps",
          "collapsed_stage.master",
          "collapsed_stage.ps"
        )
      )
    ),
    depth_cm = normalize_depth(
      coalesce_col(meta0, c("depth_cm.master", "depth_cm", "depth_cm.ps"))
    ),
    lat_block = str_squish(
      coalesce_col(meta0, c("lat_block.master", "lat_block", "lat_block.ps"))
    )
  ) %>%
  mutate(
    restoration4 = recode(restoration4, Preserved = "Conserved"),
    restoration4 = factor(restoration4, levels = stages, ordered = FALSE),
    profile_id = if_else(
      str_detect(lat_block, fixed(locality)),
      lat_block,
      paste(locality, lat_block, sep = "::")
    )
  ) %>%
  filter(sample_id %in% colnames(otu)) %>%
  arrange(locality, profile_id, depth_cm)

profile_gate <- meta %>%
  group_by(profile_id) %>%
  summarise(
    n = n(),
    n_localities = n_distinct(locality),
    n_stages = n_distinct(restoration4),
    depths = paste(sort(unique(depth_cm)), collapse = "|"),
    valid = n == 3L &&
      setequal(depth_cm, c(5, 20, 40)) &&
      n_localities == 1L &&
      n_stages == 1L,
    .groups = "drop"
  )
write_csv(profile_gate, file.path(TABLE_DIR, "202_auditoria_canon_51.csv"))

design_gate <- tibble(
  gate = c(
    "canonical_51_samples",
    "canonical_17_complete_profiles",
    "canonical_three_depths",
    "canonical_three_localities",
    "canonical_five_restoration_stages"
  ),
  pass = c(
    nrow(meta) == 51L && !anyDuplicated(meta$sample_id),
    nrow(profile_gate) == 17L && all(profile_gate$valid),
    setequal(meta$depth_cm, c(5, 20, 40)),
    n_distinct(meta$locality) == 3L,
    setequal(as.character(meta$restoration4), stages)
  ),
  detail = c(
    paste0("n=", nrow(meta)),
    paste0("profiles=", nrow(profile_gate), "; valid=", sum(profile_gate$valid)),
    paste(sort(unique(meta$depth_cm)), collapse = ","),
    paste(sort(unique(meta$locality)), collapse = ","),
    paste(levels(meta$restoration4), collapse = " -> ")
  )
)
write_csv(design_gate, file.path(TABLE_DIR, "202_compuerta_diseno.csv"))
if (!all(design_gate$pass)) {
  stop("Fallo la compuerta de diseno de 202", call. = FALSE)
}

otu <- otu[, meta$sample_id, drop = FALSE]

# -----------------------------------------------------------------------------
# Congelamiento y auditoria de los universos 84 / 12 / 3
# -----------------------------------------------------------------------------

metrics <- read_csv(METRICS_CSV, show_col_types = FALSE) %>%
  mutate(Original_Taxon_ID = as.character(Original_Taxon_ID))
conservative <- read_csv(CONSERVATIVE_CSV, show_col_types = FALSE) %>%
  mutate(Original_Taxon_ID = as.character(Original_Taxon_ID))
taxmap <- read_csv(TAXMAP_CSV, show_col_types = FALSE) %>%
  mutate(Original_Taxon_ID = as.character(Original_Taxon_ID))
residualization_models <- read_csv(
  RESIDUALIZATION_CSV,
  show_col_types = FALSE
)
run201_summary <- read_csv(RUN201_SUMMARY_CSV, show_col_types = FALSE)

required_metrics <- c(
  "Original_Taxon_ID",
  "resid_Sarle_BC",
  "resid_Dip_Q_BH",
  "observed_joint_residual",
  "observed_joint_residual_interaction"
)
required_conservative <- c(
  "Original_Taxon_ID",
  "profile_bootstrap_joint_support",
  "interaction_residual_sensitivity_pass"
)
if (length(setdiff(required_metrics, names(metrics)))) {
  stop(
    "Faltan columnas en metricas 201: ",
    paste(setdiff(required_metrics, names(metrics)), collapse = ", "),
    call. = FALSE
  )
}
if (length(setdiff(required_conservative, names(conservative)))) {
  stop(
    "Faltan columnas en tabla conservadora 201: ",
    paste(setdiff(required_conservative, names(conservative)), collapse = ", "),
    call. = FALSE
  )
}
if (!"Original_Taxon_ID" %in% names(taxmap)) {
  stop("taxmap 201 requiere Original_Taxon_ID", call. = FALSE)
}
if (!all(c("role", "restoration_used") %in% names(residualization_models))) {
  stop("201_modelos_residualizacion.csv incompleto", call. = FALSE)
}

metrics <- metrics %>%
  mutate(
    observed_joint_residual = as_logical_strict(observed_joint_residual),
    observed_joint_residual_interaction =
      as_logical_strict(observed_joint_residual_interaction)
  )
conservative <- conservative %>%
  mutate(
    interaction_residual_sensitivity_pass =
      as_logical_strict(interaction_residual_sensitivity_pass)
  )

candidate_ids <- metrics %>%
  filter(observed_joint_residual %in% TRUE) %>%
  pull(Original_Taxon_ID)
conservative_ids <- unique(conservative$Original_Taxon_ID)
robust_ids <- conservative %>%
  filter(interaction_residual_sensitivity_pass %in% TRUE) %>%
  pull(Original_Taxon_ID) %>%
  unique()

restoration_flags <- as_logical_strict(
  residualization_models$restoration_used
)

universe_gate <- tibble(
  gate = c(
    "metrics_ids_unique",
    "candidate84_exact_count",
    "candidate84_matches_prespecified_observed_rule",
    "candidate84_all_map_to_phyloseq",
    "candidate84_all_map_to_taxonomy",
    "conservative12_ids_unique",
    "conservative12_exact_count",
    "conservative12_subset_of_candidate84",
    "conservative12_support_at_least_0_80",
    "robust3_exact_count",
    "robust3_subset_of_conservative12",
    "restoration_excluded_from_201_residualization_models"
  ),
  pass = c(
    !anyDuplicated(metrics$Original_Taxon_ID),
    length(candidate_ids) == EXPECTED_CANDIDATES &&
      !anyDuplicated(candidate_ids),
    all(
      metrics$resid_Sarle_BC[metrics$observed_joint_residual %in% TRUE] > 0.555 &
        metrics$resid_Dip_Q_BH[metrics$observed_joint_residual %in% TRUE] <= 0.05
    ),
    all(candidate_ids %in% rownames(otu)),
    all(candidate_ids %in% taxmap$Original_Taxon_ID),
    !anyDuplicated(conservative$Original_Taxon_ID),
    length(conservative_ids) == EXPECTED_CONSERVATIVE,
    all(conservative_ids %in% candidate_ids),
    all(
      is.finite(conservative$profile_bootstrap_joint_support) &
        conservative$profile_bootstrap_joint_support >= 0.80
    ),
    length(robust_ids) == EXPECTED_ROBUST,
    all(robust_ids %in% conservative_ids),
    length(restoration_flags) > 0L &&
      !anyNA(restoration_flags) &&
      all(!restoration_flags)
  ),
  detail = c(
    paste0("n_metrics=", nrow(metrics)),
    paste0("observed_additive_residual=", length(candidate_ids)),
    "resid_Sarle_BC>0.555 and resid_Dip_Q_BH<=0.05",
    paste0("mapped=", sum(candidate_ids %in% rownames(otu)), "/", length(candidate_ids)),
    paste0("mapped=", sum(candidate_ids %in% taxmap$Original_Taxon_ID), "/", length(candidate_ids)),
    paste0("duplicated=", anyDuplicated(conservative$Original_Taxon_ID)),
    paste0("bootstrap_conservative=", length(conservative_ids)),
    paste0("shared=", length(intersect(conservative_ids, candidate_ids))),
    "profile bootstrap joint support >= 0.80",
    paste0("conservative_plus_interaction=", length(robust_ids)),
    paste0("shared=", length(intersect(robust_ids, conservative_ids))),
    paste0("models=", nrow(residualization_models), "; all_false=", all(!restoration_flags))
  )
)
write_csv(universe_gate, file.path(TABLE_DIR, "202_compuerta_universos_taxonomicos.csv"))
if (!all(universe_gate$pass)) {
  stop("Fallo la compuerta de universos 84/12/3", call. = FALSE)
}

metric_columns <- intersect(
  c(
    "Original_Taxon_ID",
    "prevalence_n",
    "total_count",
    "raw_Sarle_BC",
    "raw_Dip_Q_BH",
    "resid_Sarle_BC",
    "resid_Dip_Q_BH",
    "resid_interaction_Sarle_BC",
    "resid_interaction_Dip_Q_BH",
    "observed_joint_raw",
    "observed_joint_residual",
    "observed_joint_residual_interaction"
  ),
  names(metrics)
)
support_columns <- intersect(
  c(
    "Original_Taxon_ID",
    "profile_bootstrap_joint_support",
    "Stability_Score",
    "bootstrap_support_MC_SE",
    "bootstrap_support_Wilson95_lower",
    "bootstrap_support_Wilson95_upper",
    "interaction_residual_sensitivity_pass",
    "Bimodality_Class"
  ),
  names(conservative)
)

membership <- tibble(Original_Taxon_ID = candidate_ids) %>%
  left_join(taxmap, by = "Original_Taxon_ID") %>%
  left_join(metrics[, metric_columns, drop = FALSE], by = "Original_Taxon_ID") %>%
  left_join(
    conservative[, support_columns, drop = FALSE],
    by = "Original_Taxon_ID"
  ) %>%
  mutate(
    candidate84 = TRUE,
    conservative12 = Original_Taxon_ID %in% conservative_ids,
    robust3 = Original_Taxon_ID %in% robust_ids,
    in_phyloseq = Original_Taxon_ID %in% rownames(otu),
    downstream_role = case_when(
      robust3 ~ "candidate84_plus_bootstrap12_plus_interaction3",
      conservative12 ~ "candidate84_plus_bootstrap12",
      observed_joint_residual_interaction %in% TRUE ~
        "candidate84_plus_interaction_sensitivity",
      TRUE ~ "candidate84_only"
    )
  ) %>%
  arrange(
    desc(conservative12),
    desc(robust3),
    desc(profile_bootstrap_joint_support),
    resid_Dip_Q_BH
  )

write_csv(membership, file.path(TABLE_DIR, "202_membresia_universos_84_12_3.csv"))
# Nombre conservado para compatibilidad con 301; ahora contiene exactamente 84.
write_csv(membership, file.path(TABLE_DIR, "202_auditoria_mapeo_taxa.csv"))

write_csv(
  tibble(
    universe = c(
      "full_species_universe",
      "residual_nonunimodal_candidates_84",
      "bootstrap_conservative_subset_12",
      "interaction_robust_subset_3"
    ),
    n_taxa = c(
      nrow(otu),
      length(candidate_ids),
      length(conservative_ids),
      length(robust_ids)
    ),
    manuscript_role = c(
      "first_primary_taxonomic_context",
      "second_focal_sentinel_candidate_analysis",
      "ecological_interpretation_and_sensitivity",
      "strongest_cross_residualization_evidence"
    ),
    selection_description = c(
      "all species in canonical phyloseq",
      "observed joint criteria after CLR residualization by locality plus depth",
      "candidate84 plus profile-bootstrap BH joint support at least 0.80",
      "bootstrap12 plus observed locality-by-depth residual sensitivity"
    )
  ),
  file.path(TABLE_DIR, "202_definicion_y_orden_universos.csv")
)

# -----------------------------------------------------------------------------
# Matrices CLR y PCA descriptivo para 84 y 12
# -----------------------------------------------------------------------------

counts_candidates <- t(otu[candidate_ids, meta$sample_id, drop = FALSE])
counts_conservative <- t(otu[conservative_ids, meta$sample_id, drop = FALSE])
clr_candidates <- clr_rows(counts_candidates, PSEUDOCOUNT)
clr_conservative <- clr_rows(counts_conservative, PSEUDOCOUNT)

clr_full <- clr_rows(t(otu[, meta$sample_id, drop = FALSE]), PSEUDOCOUNT)
clr_full_reference_candidates <- clr_full[, candidate_ids, drop = FALSE]
clr_full_reference_conservative <- clr_full[, conservative_ids, drop = FALSE]

pca84 <- make_oriented_pca(clr_candidates, meta, "BMA84")
pca12 <- make_oriented_pca(clr_conservative, meta, "BMA12")

scores84 <- bind_cols(meta, pca84$scores)
scores12 <- bind_cols(meta, pca12$scores)
loadings84 <- pca84$loadings %>%
  left_join(membership, by = "Original_Taxon_ID")
loadings12 <- pca12$loadings %>%
  left_join(membership, by = "Original_Taxon_ID")

write_csv(meta, file.path(TABLE_DIR, "202_metadata_alineada.csv"))
write_csv(
  as_tibble(clr_candidates, rownames = "sample_id"),
  file.path(TABLE_DIR, "202_CLR_candidates84_muestras_por_taxa.csv")
)
write_csv(
  as_tibble(t(clr_candidates), rownames = "Original_Taxon_ID"),
  file.path(TABLE_DIR, "202_CLR_candidates84_taxa_por_muestras.csv")
)
write_csv(
  as_tibble(clr_conservative, rownames = "sample_id"),
  file.path(TABLE_DIR, "202_CLR_conservative12_muestras_por_taxa.csv")
)
write_csv(
  as_tibble(t(clr_conservative), rownames = "Original_Taxon_ID"),
  file.path(TABLE_DIR, "202_CLR_conservative12_taxa_por_muestras.csv")
)
write_csv(
  as_tibble(clr_full_reference_candidates, rownames = "sample_id"),
  file.path(TABLE_DIR, "202_CLR_candidates84_referencia_universo_completo.csv")
)
write_csv(
  as_tibble(clr_full_reference_conservative, rownames = "sample_id"),
  file.path(TABLE_DIR, "202_CLR_conservative12_referencia_universo_completo.csv")
)
write_csv(scores84, file.path(TABLE_DIR, "202_eje_candidates84_por_muestra.csv"))
write_csv(scores12, file.path(TABLE_DIR, "202_eje_conservative12_por_muestra.csv"))
write_csv(loadings84, file.path(TABLE_DIR, "202_loadings_candidates84.csv"))
write_csv(loadings12, file.path(TABLE_DIR, "202_loadings_conservative12.csv"))
write_csv(pca84$variance, file.path(TABLE_DIR, "202_varianza_PCA_candidates84.csv"))
write_csv(pca12$variance, file.path(TABLE_DIR, "202_varianza_PCA_conservative12.csv"))

# -----------------------------------------------------------------------------
# Figuras diagnosticas; la inferencia y figura narrativa corresponden a 301
# -----------------------------------------------------------------------------

p84 <- ggplot(
  scores84,
  aes(restoration4, BMA84_PC1, color = locality, shape = factor(depth_cm))
) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_line(aes(group = profile_id), color = "grey75", alpha = 0.45) +
  geom_point(size = 2.4, alpha = 0.9) +
  labs(
    x = "Restoration stage",
    y = "Candidate-84 axis 1 (descriptive)",
    shape = "Depth (cm)",
    color = "Locality"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

p12 <- ggplot(
  scores12,
  aes(restoration4, BMA12_PC1, color = locality, shape = factor(depth_cm))
) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_line(aes(group = profile_id), color = "grey75", alpha = 0.45) +
  geom_point(size = 2.4, alpha = 0.9) +
  labs(
    x = "Restoration stage",
    y = "Bootstrap-12 axis 1 (descriptive sensitivity)",
    shape = "Depth (cm)",
    color = "Locality"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(
  file.path(PLOT_DIR, "202_eje_candidates84_descriptivo.png"),
  p84,
  width = 8.2,
  height = 5,
  dpi = 300
)
ggsave(
  file.path(PLOT_DIR, "202_eje_candidates84_descriptivo.pdf"),
  p84,
  width = 8.2,
  height = 5
)
ggsave(
  file.path(PLOT_DIR, "202_eje_conservative12_descriptivo.png"),
  p12,
  width = 8.2,
  height = 5,
  dpi = 300
)
ggsave(
  file.path(PLOT_DIR, "202_eje_conservative12_descriptivo.pdf"),
  p12,
  width = 8.2,
  height = 5
)

# -----------------------------------------------------------------------------

# PRIMARY_RAW1031_EXTENSION
# Universo principal actualizado; candidates84 se conserva como comparacion.
cat("\n===== Universo principal raw1031 =====\n")

stopifnot(all(c(
  "observed_joint_raw", "raw_Sarle_BC", "raw_Dip_Q_BH",
  "Original_Taxon_ID"
) %in% names(metrics)))

raw_flag <- as_logical_strict(metrics$observed_joint_raw)
raw_metrics <- metrics[which(raw_flag %in% TRUE), , drop = FALSE]
raw_ids <- as.character(raw_metrics$Original_Taxon_ID)

if (length(raw_ids) != 1031L || anyNA(raw_ids) ||
    anyDuplicated(raw_ids) ||
    !all(raw_ids %in% rownames(otu)) ||
    !all(candidate_ids %in% raw_ids)) {
  stop("Universo raw1031 o inclusion de los 84 incompatible.", call. = FALSE)
}
if (!all(is.finite(raw_metrics$raw_Sarle_BC)) ||
    !all(is.finite(raw_metrics$raw_Dip_Q_BH)) ||
    !all(raw_metrics$raw_Sarle_BC > 0.555 &
         raw_metrics$raw_Dip_Q_BH <= 0.05)) {
  stop("Los 1031 no cumplen el criterio crudo esperado.", call. = FALSE)
}

raw_taxmap <- taxmap[taxmap$Original_Taxon_ID %in% raw_ids, , drop = FALSE]
if (nrow(raw_taxmap) != 1031L ||
    anyDuplicated(raw_taxmap$Original_Taxon_ID)) {
  stop("Mapa taxonomico incompleto o ambiguo para raw1031.", call. = FALSE)
}

raw_membership <- raw_metrics %>%
  mutate(
    primary_raw1031 = TRUE,
    residual84 = Original_Taxon_ID %in% candidate_ids
  )
write_csv(raw_membership, file.path(
  TABLE_DIR, "202_membresia_raw1031.csv"
))

raw_counts_all <- t(otu[raw_ids, meta$sample_id, drop = FALSE])
raw_audit <- meta %>%
  mutate(
    raw1031_total = rowSums(raw_counts_all),
    raw1031_zero = raw1031_total == 0
  )

excluded_raw_profiles <- unique(
  raw_audit$profile_id[raw_audit$raw1031_zero]
)
raw_audit <- raw_audit %>%
  mutate(
    excluded_complete_profile = profile_id %in% excluded_raw_profiles
  )
write_csv(raw_audit, file.path(
  TABLE_DIR, "202_auditoria_composicion_raw1031.csv"
))

meta_raw1031 <- meta[
  !meta$profile_id %in% excluded_raw_profiles, , drop = FALSE
]
if (nrow(meta_raw1031) < 6L ||
    !all(table(meta_raw1031$profile_id) == 3L)) {
  stop("Poblacion raw1031 insuficiente o perfiles incompletos.", call. = FALSE)
}

counts_raw1031 <- raw_counts_all[
  meta_raw1031$sample_id, , drop = FALSE
]
stopifnot(all(rowSums(counts_raw1031) > 0))
clr_raw1031 <- clr_rows(counts_raw1031, PSEUDOCOUNT)
clr_raw1031_fullref <- clr_full[
  meta_raw1031$sample_id, raw_ids, drop = FALSE
]

pca_raw1031 <- make_oriented_pca(
  clr_raw1031, meta_raw1031, "BMA1031"
)
stopifnot(
  ncol(clr_raw1031) == 1031L,
  all(is.finite(clr_raw1031)),
  max(abs(rowMeans(clr_raw1031))) < 1e-10,
  identical(rownames(pca_raw1031$fit$x), meta_raw1031$sample_id)
)

write_csv(meta_raw1031, file.path(
  TABLE_DIR, "202_metadata_raw1031_alineada.csv"
))
write_csv(
  as_tibble(counts_raw1031, rownames = "sample_id"),
  file.path(TABLE_DIR, "202_counts_raw1031_muestras_por_taxa.csv")
)
write_csv(
  as_tibble(clr_raw1031, rownames = "sample_id"),
  file.path(TABLE_DIR, "202_CLR_raw1031_muestras_por_taxa.csv")
)
write_csv(
  as_tibble(t(clr_raw1031), rownames = "Original_Taxon_ID"),
  file.path(TABLE_DIR, "202_CLR_raw1031_taxa_por_muestras.csv")
)
write_csv(
  as_tibble(clr_raw1031_fullref, rownames = "sample_id"),
  file.path(TABLE_DIR, "202_CLR_raw1031_referencia_universo_completo.csv")
)
write_csv(
  bind_cols(meta_raw1031, pca_raw1031$scores),
  file.path(TABLE_DIR, "202_eje_raw1031_por_muestra.csv")
)
write_csv(
  pca_raw1031$loadings %>%
    left_join(raw_membership, by = "Original_Taxon_ID"),
  file.path(TABLE_DIR, "202_loadings_raw1031.csv")
)
write_csv(
  pca_raw1031$variance,
  file.path(TABLE_DIR, "202_varianza_PCA_raw1031.csv")
)

write_csv(
  tibble(
    universe = c("raw1031", "residual84"),
    role = c("primary_updated", "comparison_previous"),
    n_taxa = c(length(raw_ids), length(candidate_ids)),
    selection = c("observed_joint_raw", "observed_joint_residual")
  ),
  file.path(TABLE_DIR, "202_contrato_universo_principal_actualizado.csv")
)

writeLines(c(
  "Universo principal actualizado: 1031 candidatos por criterio CLR crudo.",
  "Los 84 residuales se conservan como comparacion.",
  "Este archivo actualiza la jerarquia de los reportes anteriores de 202.",
  paste("Muestras validas:", nrow(meta_raw1031)),
  paste("Perfiles validos:", n_distinct(meta_raw1031$profile_id)),
  paste("Muestras con suma cero:", sum(raw_audit$raw1031_zero)),
  paste("Perfiles excluidos:", length(excluded_raw_profiles)),
  paste("Taxa variables usados en PCA:", pca_raw1031$n_variable_taxa),
  "No se ha recalculado bootstrap de bimodalidad para el universo ampliado."
), file.path(OUT_DIR, "202_report_raw1031_primary.txt"))

cat("RAW1031 PASS\n")
cat("Taxa:", length(raw_ids), "\n")
cat("Muestras validas:", nrow(meta_raw1031), "\n")
cat("Perfiles validos:", n_distinct(meta_raw1031$profile_id), "\n")
cat("Muestras con suma cero:", sum(raw_audit$raw1031_zero), "\n")
cat("Perfiles excluidos:", length(excluded_raw_profiles), "\n")
cat("Taxa variables PCA:", pca_raw1031$n_variable_taxa, "\n")


# Compuerta de outputs, procedencia y cierre transaccional
# -----------------------------------------------------------------------------

output_gate <- tibble(
  gate = c(
    "candidate84_CLR_dimensions",
    "conservative12_CLR_dimensions",
    "candidate84_CLR_finite",
    "conservative12_CLR_finite",
    "candidate84_CLR_rows_centered",
    "conservative12_CLR_rows_centered",
    "candidate84_PCA_sample_alignment",
    "conservative12_PCA_sample_alignment",
    "loadings_membership_complete"
  ),
  pass = c(
    identical(dim(clr_candidates), c(51L, EXPECTED_CANDIDATES)),
    identical(dim(clr_conservative), c(51L, EXPECTED_CONSERVATIVE)),
    all(is.finite(clr_candidates)),
    all(is.finite(clr_conservative)),
    max(abs(rowMeans(clr_candidates))) < 1e-10,
    max(abs(rowMeans(clr_conservative))) < 1e-10,
    identical(rownames(pca84$fit$x), meta$sample_id),
    identical(rownames(pca12$fit$x), meta$sample_id),
    nrow(loadings84) == EXPECTED_CANDIDATES &&
      sum(loadings84$conservative12) == EXPECTED_CONSERVATIVE &&
      sum(loadings84$robust3) == EXPECTED_ROBUST
  ),
  detail = c(
    paste(dim(clr_candidates), collapse = "x"),
    paste(dim(clr_conservative), collapse = "x"),
    paste0("nonfinite=", sum(!is.finite(clr_candidates))),
    paste0("nonfinite=", sum(!is.finite(clr_conservative))),
    paste0("max_abs_row_mean=", signif(max(abs(rowMeans(clr_candidates))), 6)),
    paste0("max_abs_row_mean=", signif(max(abs(rowMeans(clr_conservative))), 6)),
    paste0("n=", nrow(pca84$fit$x)),
    paste0("n=", nrow(pca12$fit$x)),
    paste0(
      "loadings84=", nrow(loadings84),
      "; conservative=", sum(loadings84$conservative12),
      "; robust=", sum(loadings84$robust3)
    )
  )
)
write_csv(output_gate, file.path(TABLE_DIR, "202_compuerta_outputs.csv"))
if (!all(output_gate$pass)) {
  stop("Fallo la compuerta de outputs de 202", call. = FALSE)
}

input_paths <- c(
  PS_RDS,
  MASTER_CSV,
  METRICS_CSV,
  CONSERVATIVE_CSV,
  TAXMAP_CSV,
  RESIDUALIZATION_CSV,
  RUN201_SUMMARY_CSV
)
input_md5 <- tibble(
  input = c(
    "ps_rds",
    "master_csv",
    "metrics_csv",
    "conservative_csv",
    "taxmap_csv",
    "residualization_csv",
    "run201_summary_csv"
  ),
  path = input_paths,
  md5 = unname(as.character(tools::md5sum(input_paths)))
)
write_csv(input_md5, file.path(TABLE_DIR, "202_input_md5.csv"))

run_info <- tibble(
  parameter = c(
    "script",
    "run001",
    "run201",
    "ps_rds",
    "master_csv",
    "pseudocount",
    "n_samples",
    "n_profiles",
    "n_full_taxa",
    "n_candidate84",
    "n_conservative12",
    "n_robust3",
    "candidate84_definition",
    "conservative12_definition",
    "manuscript_order",
    "candidate84_PC1_variance",
    "conservative12_PC1_variance",
    "candidate84_orientation_flip",
    "conservative12_orientation_flip",
    "stage_inference",
    "decision_timing",
    "latest_attempt_file",
    "latest_file"
  ),
  value = as.character(c(
    SCRIPT_ID,
    RUN001,
    RUN201,
    PS_RDS,
    MASTER_CSV,
    PSEUDOCOUNT,
    nrow(meta),
    n_distinct(meta$profile_id),
    nrow(otu),
    length(candidate_ids),
    length(conservative_ids),
    length(robust_ids),
    "observed joint criteria after CLR residualization by locality plus depth",
    "candidate84 with profile-bootstrap BH joint support at least 0.80",
    "full community first; candidate84 second; conservative12 interpretation and sensitivity",
    pca84$variance$variance_explained[[1L]],
    pca12$variance$variance_explained[[1L]],
    pca84$flip,
    pca12$flip,
    "none_descriptive_only",
    "84/12 hierarchy frozen before downstream ecological tests",
    LATEST_ATTEMPT_FILE,
    LATEST_FILE
  ))
)
write_csv(run_info, file.path(TABLE_DIR, "202_parametros_y_resumen.csv"))

report <- c(
  "202 preparation of taxonomic universes",
  "=======================================",
  "",
  paste0("RUN001: ", RUN001),
  paste0("RUN201: ", RUN201),
  paste0("Output: ", OUT_DIR),
  "",
  "Frozen hierarchy:",
  paste0("  - full phyloseq universe: ", nrow(otu), " taxa"),
  paste0("  - residual non-unimodal candidates: ", length(candidate_ids), " taxa"),
  paste0("  - bootstrap-conservative subset: ", length(conservative_ids), " taxa"),
  paste0("  - interaction-robust subset: ", length(robust_ids), " taxa"),
  "",
  "Manuscript order:",
  "  1. full-community Hill diversity and taxonomic ordination",
  "  2. test whether the candidate-84 subset retains ecological signal",
  "  3. use the conservative-12 subset for ecological interpretation and sensitivity",
  "",
  "Interpretation limit:",
  "  - 202 performs no restoration-stage inference",
  "  - PCA sign orientation does not change distances, loadings magnitude or fit",
  "  - candidate84 are not called bootstrap-conservative taxa"
)
writeLines(report, file.path(OUT_DIR, "202_report.txt"))

writeLines(
  capture.output(sessionInfo()),
  file.path(LOG_DIR, "202_sessionInfo.txt")
)

writeLines(OUT_DIR, LATEST_FILE)

cat("===== DONE 202 =====\n")
cat("Full species universe:", nrow(otu), "\n")
cat("Residual non-unimodal candidates:", length(candidate_ids), "\n")
cat("Bootstrap-conservative subset:", length(conservative_ids), "\n")
cat("Interaction-robust subset:", length(robust_ids), "\n")
cat("Candidate84 PC1 variance:", pca84$variance$variance_explained[[1L]], "\n")
cat("Conservative12 PC1 variance:", pca12$variance$variance_explained[[1L]], "\n")
cat("Output:", OUT_DIR, "\n")
cat("Latest:", LATEST_FILE, "\n")
cat("Finished UTC:", format(Sys.time(), tz = "UTC"), "\n")
