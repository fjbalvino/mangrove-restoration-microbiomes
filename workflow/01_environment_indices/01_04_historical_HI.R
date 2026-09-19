#!/usr/bin/env Rscript
# ============================================================
# 01_04_historical_HI.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Archived Model A defines current HI: 1,467 taxa, categorical locality/depth, unscaled residual PCA. Other models and tests are historical.
# Inputs (source expressions; complete list in docs/contracts/01_04_historical_HI.json):
#   DIR_002 <- readLines(LATEST_002_FILE, warn = FALSE)[1]
#   clr_taxa <- read_csv(CLR_TAXA_BY_SAMPLES, show_col_types = FALSE)
#   meta_raw <- read_csv(META_JOINED, show_col_types = FALSE)
# Outputs (source expressions; complete list in contract):
#   writeLines(OUT_DIR, file.path(BASE_RESULTS, paste0("LATEST_", SCRIPT_ID, ".txt")))
#   write_csv(axis_tbl, file.path(out_dir, paste0("003_", prefix, "_axis_by_sample.csv")))
#   write_csv(loadings_tbl, file.path(out_dir, paste0("003_", prefix, "_loadings.csv")))
#   write_csv(variance_tbl, file.path(out_dir, paste0("003_", prefix, "_variance_explained.csv")))
#   write_csv(stats_tbl, file.path(out_dir, paste0("003_", prefix, "_nonparametric_tests.csv")))
#   write_csv(lm_anova, file.path(out_dir, paste0("003_", prefix, "_lm_anova.csv")))
#   write_csv(pairwise_tbl, file.path(out_dir, paste0("003_", prefix, "_pairwise_wilcox_stage.csv")))
#   write_csv(adonis_tbl, file.path(out_dir, paste0("003_", prefix, "_adonis2_margin.csv")))
#   ggsave(
#   write_csv(meta_use, file.path(OUT_DIR, "003_metadata_used.csv"))
# Algorithmic provenance:
# Historical Model A: categorical locality/depth residualisation and centred unscaled PCA.
#   Aitchison (1982), doi:10.1111/j.2517-6161.1982.tb01195.x. Study-specific index calibration is not an externally validated recovery index.
# Source SHA-256: 8a6962c5c6d914e4bbe53ad43bc48fe8f6d318f6da78153314bbc510e29d1eb5
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


# ============================================================
# Script 003: residualized_bimodal_axis_locality_depth_ECI
# Proyecto: Tipping_points
#
# Objetivo:
#   Reconstruir ejes comunitarios residuales de especies
#   bimodales conservadoras después de remover:
#
#   Modelo A: locality + depth_cm
#   Modelo B: locality + depth_cm + ECI, si existe
#
# Pregunta:
#   ¿queda una señal sucesional Degraded-Intermediate-Preserved
#   después de controlar estructura espacial, vertical y ambiental?
#
# Inputs:
#   - Outputs del Script 002
#
# Outputs:
#   - Ejes residuales por muestra
#   - Loadings residuales
#   - Estadística por estadio
#   - PERMANOVA marginal
#   - Figuras
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# -----------------------------
# 1. Parámetros
# -----------------------------

BASE_RESULTS <- "/home/fjbalvino/Tipping_points/Results"

SCRIPT_ID <- "003_residualized_bimodal_axis_locality_depth_ECI"
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

OUT_DIR <- file.path(BASE_RESULTS, paste0(SCRIPT_ID, "_", TIMESTAMP))
PLOT_DIR <- file.path(OUT_DIR, "plots")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

writeLines(OUT_DIR, file.path(BASE_RESULTS, paste0("LATEST_", SCRIPT_ID, ".txt")))

LOG_FILE <- file.path(OUT_DIR, "log.txt")
sink(LOG_FILE, split = TRUE)

STAGE_LEVELS <- c("Degraded", "Intermediate", "Preserved")
SEED <- 123
set.seed(SEED)

# Usar último output del Script 002
LATEST_002_FILE <- file.path(BASE_RESULTS, "LATEST_002_build_bimodal_species_axis.txt")

if (!file.exists(LATEST_002_FILE)) {
  stop("No existe LATEST_002_build_bimodal_species_axis.txt")
}

DIR_002 <- readLines(LATEST_002_FILE, warn = FALSE)[1]

CLR_TAXA_BY_SAMPLES <- file.path(DIR_002, "002_clr_bimodal_species_taxa_by_samples.csv")
META_JOINED <- file.path(DIR_002, "002_meta_axis_after_join.csv")
LOADINGS_002 <- file.path(DIR_002, "002_bimodal_species_axis_loadings.csv")

cat("============================================================\n")
cat("Script:", SCRIPT_ID, "\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("Output:", OUT_DIR, "\n")
cat("Using Script 002 output:", DIR_002, "\n")
cat("============================================================\n\n")

# -----------------------------
# 2. Funciones
# -----------------------------

coalesce_existing <- function(df, candidates) {
  candidates <- candidates[candidates %in% colnames(df)]

  if (length(candidates) == 0) {
    return(rep(NA_character_, nrow(df)))
  }

  out <- as.character(df[[candidates[1]]])

  if (length(candidates) > 1) {
    for (nm in candidates[-1]) {
      out <- dplyr::coalesce(out, as.character(df[[nm]]))
    }
  }

  out
}

find_numeric_covariate <- function(df, candidates) {
  candidates <- candidates[candidates %in% colnames(df)]

  if (length(candidates) == 0) {
    return(NULL)
  }

  for (nm in candidates) {
    x <- suppressWarnings(as.numeric(df[[nm]]))
    if (sum(is.finite(x)) >= 10 && stats::sd(x, na.rm = TRUE) > 0) {
      return(nm)
    }
  }

  NULL
}

residualize_matrix <- function(mat_samples_by_taxa, meta, formula_rhs) {
  stopifnot(nrow(mat_samples_by_taxa) == nrow(meta))

  residual_mat <- matrix(
    NA_real_,
    nrow = nrow(mat_samples_by_taxa),
    ncol = ncol(mat_samples_by_taxa),
    dimnames = dimnames(mat_samples_by_taxa)
  )

  for (j in seq_len(ncol(mat_samples_by_taxa))) {
    y <- mat_samples_by_taxa[, j]

    df <- meta
    df$y <- as.numeric(y)

    fit <- tryCatch(
      lm(as.formula(paste("y ~", formula_rhs)), data = df),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      residual_mat[, j] <- NA_real_
    } else {
      residual_mat[, j] <- residuals(fit)
    }
  }

  keep <- apply(residual_mat, 2, function(z) all(is.finite(z)) && sd(z) > 0)
  residual_mat[, keep, drop = FALSE]
}

run_residual_pca <- function(resid_mat, meta, model_label, out_dir, plot_dir) {
  cat("\n[MODEL]", model_label, "\n")
  cat("Residual matrix samples x taxa:", nrow(resid_mat), "x", ncol(resid_mat), "\n")

  pca <- prcomp(resid_mat, center = TRUE, scale. = FALSE)

  var_explained <- (pca$sdev^2) / sum(pca$sdev^2)

  axis_tbl <- meta %>%
    mutate(
      model = model_label,
      RBMA1_raw = pca$x[, 1],
      RBMA2_raw = pca$x[, 2],
      RBMA3_raw = pca$x[, 3]
    )

  med_degraded <- median(axis_tbl$RBMA1_raw[axis_tbl$collapsed_stage == "Degraded"], na.rm = TRUE)
  med_preserved <- median(axis_tbl$RBMA1_raw[axis_tbl$collapsed_stage == "Preserved"], na.rm = TRUE)

  flip <- ifelse(
    is.finite(med_degraded) &&
      is.finite(med_preserved) &&
      med_preserved < med_degraded,
    -1,
    1
  )

  axis_tbl <- axis_tbl %>%
    mutate(
      RBMA1 = flip * RBMA1_raw,
      RBMA2 = RBMA2_raw,
      RBMA3 = RBMA3_raw
    )

  loadings_tbl <- as.data.frame(pca$rotation[, 1:min(5, ncol(pca$rotation)), drop = FALSE]) %>%
    rownames_to_column("Original_Taxon_ID") %>%
    mutate(
      model = model_label,
      Original_Taxon_ID = as.character(Original_Taxon_ID),
      PC1 = flip * PC1
    ) %>%
    arrange(desc(abs(PC1)))

  variance_tbl <- tibble(
    model = model_label,
    axis = paste0("PC", seq_along(var_explained)),
    variance_explained = var_explained,
    cumulative_variance = cumsum(var_explained)
  )

  kw_stage <- kruskal.test(RBMA1 ~ collapsed_stage, data = axis_tbl)
  kw_depth <- kruskal.test(RBMA1 ~ depth_cm, data = axis_tbl)
  kw_locality <- kruskal.test(RBMA1 ~ locality, data = axis_tbl)

  lm_fit <- lm(RBMA1 ~ collapsed_stage + depth_cm + locality, data = axis_tbl)
  lm_anova <- as.data.frame(anova(lm_fit)) %>%
    rownames_to_column("term") %>%
    mutate(model = model_label)

  pairwise_stage <- pairwise.wilcox.test(
    x = axis_tbl$RBMA1,
    g = axis_tbl$collapsed_stage,
    p.adjust.method = "BH",
    exact = FALSE
  )

  pairwise_tbl <- as.data.frame(as.table(pairwise_stage$p.value)) %>%
    filter(!is.na(Freq)) %>%
    rename(group1 = Var1, group2 = Var2, p_adj_BH = Freq) %>%
    mutate(model = model_label)

  stats_tbl <- tibble(
    model = model_label,
    test = c(
      "Kruskal_RBMA1_by_collapsed_stage",
      "Kruskal_RBMA1_by_depth_cm",
      "Kruskal_RBMA1_by_locality"
    ),
    statistic = c(
      unname(kw_stage$statistic),
      unname(kw_depth$statistic),
      unname(kw_locality$statistic)
    ),
    p_value = c(
      kw_stage$p.value,
      kw_depth$p.value,
      kw_locality$p.value
    ),
    RBMA1_variance_explained = var_explained[1],
    RBMA2_variance_explained = var_explained[2],
    orientation_flip = flip
  )

  cat("RBMA1 variance explained:", round(var_explained[1], 4), "\n")
  cat("RBMA2 variance explained:", round(var_explained[2], 4), "\n")
  cat("Kruskal RBMA1 ~ collapsed_stage p:", kw_stage$p.value, "\n")
  cat("Kruskal RBMA1 ~ depth_cm p:", kw_depth$p.value, "\n")
  cat("Kruskal RBMA1 ~ locality p:", kw_locality$p.value, "\n")

  # PERMANOVA marginal sobre matriz residual
  adonis_tbl <- NULL

  if (requireNamespace("vegan", quietly = TRUE)) {
    dist_resid <- dist(resid_mat, method = "euclidean")

    adonis_fit <- vegan::adonis2(
      dist_resid ~ collapsed_stage + depth_cm + locality,
      data = axis_tbl,
      permutations = 999,
      by = "margin"
    )

    adonis_tbl <- as.data.frame(adonis_fit) %>%
      rownames_to_column("term") %>%
      mutate(model = model_label)

    cat("PERMANOVA marginal:\n")
    print(adonis_tbl)
  }

  prefix <- gsub("[^A-Za-z0-9]+", "_", model_label)

  write_csv(axis_tbl, file.path(out_dir, paste0("003_", prefix, "_axis_by_sample.csv")))
  write_csv(loadings_tbl, file.path(out_dir, paste0("003_", prefix, "_loadings.csv")))
  write_csv(variance_tbl, file.path(out_dir, paste0("003_", prefix, "_variance_explained.csv")))
  write_csv(stats_tbl, file.path(out_dir, paste0("003_", prefix, "_nonparametric_tests.csv")))
  write_csv(lm_anova, file.path(out_dir, paste0("003_", prefix, "_lm_anova.csv")))
  write_csv(pairwise_tbl, file.path(out_dir, paste0("003_", prefix, "_pairwise_wilcox_stage.csv")))

  if (!is.null(adonis_tbl)) {
    write_csv(adonis_tbl, file.path(out_dir, paste0("003_", prefix, "_adonis2_margin.csv")))
  }

  p_stage <- ggplot(axis_tbl, aes(x = collapsed_stage, y = RBMA1)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65) +
    geom_jitter(aes(shape = depth_cm), width = 0.14, size = 2.8, alpha = 0.85) +
    facet_wrap(~ depth_cm, scales = "free_y") +
    theme_bw(base_size = 13) +
    labs(
      title = paste0("Residual Bimodal Microbiome Axis 1: ", model_label),
      subtitle = "Testing restoration signal after residualization",
      x = "Collapsed restoration stage",
      y = "RBMA1"
    )

  ggsave(
    file.path(plot_dir, paste0("003_", prefix, "_RBMA1_by_stage_depth.png")),
    p_stage,
    width = 10,
    height = 6,
    dpi = 300
  )

  ggsave(
    file.path(plot_dir, paste0("003_", prefix, "_RBMA1_by_stage_depth.pdf")),
    p_stage,
    width = 10,
    height = 6
  )

  p_density <- ggplot(axis_tbl, aes(x = RBMA1, fill = collapsed_stage)) +
    geom_density(alpha = 0.4, na.rm = TRUE) +
    facet_wrap(~ depth_cm, scales = "free_y") +
    theme_bw(base_size = 13) +
    labs(
      title = paste0("Density of residual BMA1: ", model_label),
      subtitle = "Intermediate broadening after residualization",
      x = "RBMA1",
      y = "Density",
      fill = "Stage"
    )

  ggsave(
    file.path(plot_dir, paste0("003_", prefix, "_RBMA1_density_by_depth.png")),
    p_density,
    width = 10,
    height = 6,
    dpi = 300
  )

  ggsave(
    file.path(plot_dir, paste0("003_", prefix, "_RBMA1_density_by_depth.pdf")),
    p_density,
    width = 10,
    height = 6
  )

  p_pca <- ggplot(axis_tbl, aes(x = RBMA1, y = RBMA2)) +
    geom_point(aes(color = collapsed_stage, shape = depth_cm), size = 3.5, alpha = 0.9) +
    theme_bw(base_size = 13) +
    labs(
      title = paste0("Residual bimodal species PCA: ", model_label),
      subtitle = paste0(
        "RBMA1 = ", round(100 * var_explained[1], 1),
        "%; RBMA2 = ", round(100 * var_explained[2], 1), "%"
      ),
      x = "RBMA1",
      y = "RBMA2",
      color = "Stage",
      shape = "Depth"
    )

  if (all(table(axis_tbl$collapsed_stage) >= 3)) {
    p_pca <- p_pca +
      stat_ellipse(aes(color = collapsed_stage), type = "norm", linewidth = 0.7, na.rm = TRUE)
  }

  ggsave(
    file.path(plot_dir, paste0("003_", prefix, "_RBMA1_RBMA2_pca.png")),
    p_pca,
    width = 8,
    height = 6,
    dpi = 300
  )

  ggsave(
    file.path(plot_dir, paste0("003_", prefix, "_RBMA1_RBMA2_pca.pdf")),
    p_pca,
    width = 8,
    height = 6
  )

  list(
    axis = axis_tbl,
    loadings = loadings_tbl,
    variance = variance_tbl,
    stats = stats_tbl,
    lm_anova = lm_anova,
    pairwise = pairwise_tbl,
    adonis = adonis_tbl
  )
}

# -----------------------------
# 3. Validar inputs
# -----------------------------

input_files <- c(CLR_TAXA_BY_SAMPLES, META_JOINED)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop("Faltan archivos requeridos del Script 002:\n", paste(missing_files, collapse = "\n"))
}

# -----------------------------
# 4. Cargar datos
# -----------------------------

cat("[1] Loading Script 002 outputs\n")

clr_taxa <- read_csv(CLR_TAXA_BY_SAMPLES, show_col_types = FALSE)
meta_raw <- read_csv(META_JOINED, show_col_types = FALSE)

cat("CLR taxa rows:", nrow(clr_taxa), "\n")
cat("Meta rows:", nrow(meta_raw), "\n\n")

if (!"Original_Taxon_ID" %in% colnames(clr_taxa)) {
  stop("La matriz CLR no contiene Original_Taxon_ID.")
}

taxa_ids <- as.character(clr_taxa$Original_Taxon_ID)

clr_mat_taxa_by_samples <- clr_taxa %>%
  select(-Original_Taxon_ID) %>%
  as.data.frame()

rownames(clr_mat_taxa_by_samples) <- taxa_ids

clr_mat_taxa_by_samples <- as.matrix(clr_mat_taxa_by_samples)
storage.mode(clr_mat_taxa_by_samples) <- "numeric"

clr_samples_by_taxa <- t(clr_mat_taxa_by_samples)

# -----------------------------
# 5. Preparar metadata
# -----------------------------

cat("[2] Preparing metadata\n")

meta <- meta_raw %>%
  mutate(
    sample_id = as.character(sample_id),
    collapsed_stage = coalesce_existing(
      .,
      c("collapsed_stage", "collapsed_stage_clean", "collapsed_stage.master", "collapsed_stage.ps")
    ),
    depth_cm = coalesce_existing(
      .,
      c("depth_cm", "depth_cm_clean", "depth_cm.master", "depth_cm.ps", "Depth_cm", "depth")
    ),
    locality = coalesce_existing(
      .,
      c("locality", "locality_clean", "locality.master", "locality.ps", "Locality")
    )
  )

eci_col <- find_numeric_covariate(
  meta,
  c(
    "ECI_PC1",
    "ECI",
    "ECI1",
    "env_PC1",
    "Env_PC1",
    "environmental_PC1",
    "PC1_env",
    "CAP1_env"
  )
)

if (is.null(eci_col)) {
  cat("No ECI column detected. Model B will be skipped.\n")
} else {
  cat("Detected ECI column:", eci_col, "\n")
  meta$ECI_used <- suppressWarnings(as.numeric(meta[[eci_col]]))
}

common_samples <- intersect(rownames(clr_samples_by_taxa), meta$sample_id)

meta_use <- meta %>%
  filter(sample_id %in% common_samples) %>%
  mutate(
    collapsed_stage = factor(collapsed_stage, levels = STAGE_LEVELS),
    depth_cm = factor(depth_cm),
    locality = factor(locality)
  ) %>%
  filter(!is.na(collapsed_stage), !is.na(depth_cm), !is.na(locality)) %>%
  arrange(match(sample_id, rownames(clr_samples_by_taxa)))

clr_use <- clr_samples_by_taxa[meta_use$sample_id, , drop = FALSE]

if (!all(rownames(clr_use) == meta_use$sample_id)) {
  stop("Error de alineación entre matriz CLR y metadata.")
}

cat("Samples retained:", nrow(meta_use), "\n")
cat("Taxa retained:", ncol(clr_use), "\n")
cat("Stages:\n")
print(table(meta_use$collapsed_stage, useNA = "ifany"))
cat("Depths:\n")
print(table(meta_use$depth_cm, useNA = "ifany"))
cat("Localities:\n")
print(table(meta_use$locality, useNA = "ifany"))
cat("\n")

write_csv(meta_use, file.path(OUT_DIR, "003_metadata_used.csv"))

# -----------------------------
# 6. Modelo A: residualizar localidad + profundidad
# -----------------------------

cat("[3] Residualizing Model A: locality + depth_cm\n")

resid_A <- residualize_matrix(
  mat_samples_by_taxa = clr_use,
  meta = meta_use,
  formula_rhs = "locality + depth_cm"
)

write_csv(
  as.data.frame(t(resid_A)) %>%
    rownames_to_column("Original_Taxon_ID"),
  file.path(OUT_DIR, "003_residual_matrix_ModelA_locality_depth_taxa_by_samples.csv")
)

res_A <- run_residual_pca(
  resid_mat = resid_A,
  meta = meta_use,
  model_label = "ModelA_locality_depth",
  out_dir = OUT_DIR,
  plot_dir = PLOT_DIR
)

# -----------------------------
# 7. Modelo B: residualizar localidad + profundidad + ECI
# -----------------------------

res_B <- NULL

if (!is.null(eci_col)) {
  cat("[4] Residualizing Model B: locality + depth_cm + ECI\n")

  meta_B <- meta_use %>%
    mutate(ECI_used = suppressWarnings(as.numeric(.data[[eci_col]]))) %>%
    filter(is.finite(ECI_used))

  clr_B <- clr_use[meta_B$sample_id, , drop = FALSE]

  if (nrow(meta_B) >= 10 && length(unique(meta_B$ECI_used)) >= 5) {
    resid_B <- residualize_matrix(
      mat_samples_by_taxa = clr_B,
      meta = meta_B,
      formula_rhs = "locality + depth_cm + ECI_used"
    )

    write_csv(
      as.data.frame(t(resid_B)) %>%
        rownames_to_column("Original_Taxon_ID"),
      file.path(OUT_DIR, "003_residual_matrix_ModelB_locality_depth_ECI_taxa_by_samples.csv")
    )

    res_B <- run_residual_pca(
      resid_mat = resid_B,
      meta = meta_B,
      model_label = "ModelB_locality_depth_ECI",
      out_dir = OUT_DIR,
      plot_dir = PLOT_DIR
    )
  } else {
    cat("ECI detected but insufficient finite values. Skipping Model B.\n")
  }
}

# -----------------------------
# 8. Comparar modelos
# -----------------------------

cat("[5] Comparing residual models\n")

all_stats <- bind_rows(
  res_A$stats,
  if (!is.null(res_B)) res_B$stats else NULL
)

all_lm <- bind_rows(
  res_A$lm_anova,
  if (!is.null(res_B)) res_B$lm_anova else NULL
)

all_pairwise <- bind_rows(
  res_A$pairwise,
  if (!is.null(res_B)) res_B$pairwise else NULL
)

all_adonis <- bind_rows(
  res_A$adonis,
  if (!is.null(res_B)) res_B$adonis else NULL
)

write_csv(all_stats, file.path(OUT_DIR, "003_all_models_nonparametric_tests.csv"))
write_csv(all_lm, file.path(OUT_DIR, "003_all_models_lm_anova.csv"))
write_csv(all_pairwise, file.path(OUT_DIR, "003_all_models_pairwise_wilcox_stage.csv"))

if (nrow(all_adonis) > 0) {
  write_csv(all_adonis, file.path(OUT_DIR, "003_all_models_adonis2_margin.csv"))
}

# -----------------------------
# 9. Resumen
# -----------------------------

summary_tbl <- tibble(
  item = c(
    "script",
    "timestamp",
    "script002_dir",
    "n_samples_model_A",
    "n_taxa_model_A",
    "ECI_column_detected",
    "model_B_ran",
    "output_dir"
  ),
  value = as.character(c(
    SCRIPT_ID,
    TIMESTAMP,
    DIR_002,
    nrow(res_A$axis),
    ncol(resid_A),
    ifelse(is.null(eci_col), "none", eci_col),
    ifelse(is.null(res_B), "FALSE", "TRUE"),
    OUT_DIR
  ))
)

write_csv(summary_tbl, file.path(OUT_DIR, "003_run_summary.csv"))

cat("============================================================\n")
cat("[DONE] Residualized bimodal axis completed\n")
cat("Output directory:\n")
cat(OUT_DIR, "\n\n")
cat("Main outputs:\n")
cat(" - 003_all_models_nonparametric_tests.csv\n")
cat(" - 003_all_models_lm_anova.csv\n")
cat(" - 003_all_models_pairwise_wilcox_stage.csv\n")
cat(" - 003_all_models_adonis2_margin.csv\n")
cat(" - 003_ModelA_locality_depth_axis_by_sample.csv\n")
cat(" - 003_ModelB_locality_depth_ECI_axis_by_sample.csv, if ECI exists\n")
cat(" - plots/*RBMA1*.png\n")
cat("============================================================\n")

sink()

