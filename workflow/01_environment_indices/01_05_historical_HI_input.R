#!/usr/bin/env Rscript
# ============================================================
# 01_05_historical_HI_input.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Historical 1,467-taxon CLR producer; distinct from the current 1,031-taxon screen.
# Inputs (source expressions; complete list in docs/contracts/01_05_historical_HI_input.json):
#   ps <- readRDS(PS_RDS)
#   bimodal_tbl <- read_csv(
#   taxmap <- read_csv(
#   master <- read_csv(MASTER_CSV, show_col_types = FALSE)
# Outputs (source expressions; complete list in contract):
#   writeLines(OUT_DIR, file.path(BASE_OUT, paste0("LATEST_", SCRIPT_ID, ".txt")))
#   write_csv(bimodal_map, file.path(OUT_DIR, "002_bimodal_species_mapped_to_phyloseq.csv"))
#   write_csv(
#   write_csv(meta_axis, file.path(OUT_DIR, "002_meta_axis_after_join.csv"))
#   write_csv(axis_tbl, file.path(OUT_DIR, "002_bimodal_microbiome_axis_by_sample.csv"))
#   write_csv(loadings_tbl, file.path(OUT_DIR, "002_bimodal_species_axis_loadings.csv"))
#   write_csv(variance_tbl, file.path(OUT_DIR, "002_pca_variance_explained.csv"))
#   write_csv(stats_tbl, file.path(OUT_DIR, "002_bma1_nonparametric_tests.csv"))
#   write_csv(lm_anova_tbl, file.path(OUT_DIR, "002_bma1_lm_anova.csv"))
#   write_csv(pairwise_stage_tbl, file.path(OUT_DIR, "002_bma1_pairwise_wilcox_stage.csv"))
# Algorithmic provenance:
# Historical 1,467-taxon CLR input construction (pseudocount 1); ordination.
#   Aitchison (1982), doi:10.1111/j.2517-6161.1982.tb01195.x. Study-specific index calibration is not an externally validated recovery index.
# Source SHA-256: d70bceefb41553b512d8117b9b376c9f42df20fb2b37d1ec43fd0c785bd19828
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


# ============================================================
# Script 002: build_bimodal_species_axis
# Proyecto: Tipping_points
#
# Objetivo:
#   Construir un eje comunitario basado en especies bimodales
#   conservadoras detectadas en el Script 001.
#
# Fixes:
#   1. Subset robusto de taxa con intersect().
#   2. Manejo robusto de columnas duplicadas tras join.
#   3. Coerción explícita de Species_ID y Original_Taxon_ID a character.
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# -----------------------------
# 1. Parámetros
# -----------------------------

PS_RDS <- "/home/fjbalvino/R/Ciencia-Frontera/resultados/phyloseq_withConservedDegraded_20260130_1730/rank_S/phyloseq_FILTRADO_rank_S.rds"

BIMODAL_DIR <- "/home/fjbalvino/Tipping_points/Results/001_bootstrap_bimodality_species_universe_20260504_002025"

BIMODAL_CSV <- file.path(BIMODAL_DIR, "001_bimodal_species_conservative_joint.csv")
TAXMAP_CSV <- file.path(BIMODAL_DIR, "001_taxonomy_species_id_map.csv")

MASTER_CSV <- "/home/fjbalvino/R/Ciencia-Frontera/resultados/102_FINAL_VALIDATION_20260430_040736/102_final_stage_data.csv"

BASE_OUT <- "/home/fjbalvino/Tipping_points/Results"

SCRIPT_ID <- "002_build_bimodal_species_axis"
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

OUT_DIR <- file.path(BASE_OUT, paste0(SCRIPT_ID, "_", TIMESTAMP))
PLOT_DIR <- file.path(OUT_DIR, "plots")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

writeLines(OUT_DIR, file.path(BASE_OUT, paste0("LATEST_", SCRIPT_ID, ".txt")))

LOG_FILE <- file.path(OUT_DIR, "log.txt")
sink(LOG_FILE, split = TRUE)

PSEUDOCOUNT <- 1
SEED <- 123
STAGE_LEVELS <- c("Degraded", "Intermediate", "Preserved")

set.seed(SEED)

cat("============================================================\n")
cat("Script:", SCRIPT_ID, "\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("Output:", OUT_DIR, "\n")
cat("============================================================\n\n")

cat("[PARAMETERS]\n")
cat("PS_RDS:", PS_RDS, "\n")
cat("BIMODAL_CSV:", BIMODAL_CSV, "\n")
cat("TAXMAP_CSV:", TAXMAP_CSV, "\n")
cat("MASTER_CSV:", MASTER_CSV, "\n")
cat("PSEUDOCOUNT:", PSEUDOCOUNT, "\n\n")

# -----------------------------
# 2. Funciones auxiliares
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

safe_select_taxmap <- function(taxmap_df) {
  required <- c("Species_ID", "Original_Taxon_ID")

  missing_required <- setdiff(required, colnames(taxmap_df))
  if (length(missing_required) > 0) {
    stop("Taxmap no contiene columnas requeridas: ", paste(missing_required, collapse = ", "))
  }

  optional_cols <- intersect(c("Clean_Name", "taxonomy_full"), colnames(taxmap_df))

  taxmap_df %>%
    mutate(
      Species_ID = as.character(Species_ID),
      Original_Taxon_ID = as.character(Original_Taxon_ID)
    ) %>%
    select(all_of(c(required, optional_cols)))
}

# -----------------------------
# 3. Validar inputs
# -----------------------------

input_files <- c(PS_RDS, BIMODAL_CSV, TAXMAP_CSV, MASTER_CSV)
missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop("Faltan archivos:\n", paste(missing_files, collapse = "\n"))
}

# -----------------------------
# 4. Cargar objetos
# -----------------------------

cat("[1] Loading inputs\n")

ps <- readRDS(PS_RDS)

bimodal_tbl <- read_csv(
  BIMODAL_CSV,
  show_col_types = FALSE,
  col_types = cols(
    Species_ID = col_character(),
    Original_Taxon_ID = col_character(),
    .default = col_guess()
  )
)

taxmap <- read_csv(
  TAXMAP_CSV,
  show_col_types = FALSE,
  col_types = cols(
    Species_ID = col_character(),
    Original_Taxon_ID = col_character(),
    .default = col_guess()
  )
)

master <- read_csv(MASTER_CSV, show_col_types = FALSE)

bimodal_tbl <- bimodal_tbl %>%
  mutate(
    Species_ID = as.character(Species_ID),
    Original_Taxon_ID = if ("Original_Taxon_ID" %in% colnames(.)) {
      as.character(Original_Taxon_ID)
    } else {
      NA_character_
    }
  )

taxmap <- taxmap %>%
  mutate(
    Species_ID = as.character(Species_ID),
    Original_Taxon_ID = as.character(Original_Taxon_ID)
  )

cat("Phyloseq samples:", phyloseq::nsamples(ps), "\n")
cat("Phyloseq taxa:", phyloseq::ntaxa(ps), "\n")
cat("Bimodal conservative species:", nrow(bimodal_tbl), "\n")
cat("Master rows:", nrow(master), "\n\n")

required_master <- c("sample_id", "collapsed_stage", "depth_cm", "locality")
missing_master <- setdiff(required_master, colnames(master))

if (length(missing_master) > 0) {
  stop("Faltan columnas en MASTER_CSV: ", paste(missing_master, collapse = ", "))
}

if (!"Species_ID" %in% colnames(bimodal_tbl)) {
  stop("BIMODAL_CSV no contiene Species_ID.")
}

if (!all(c("Species_ID", "Original_Taxon_ID") %in% colnames(taxmap))) {
  stop("TAXMAP_CSV debe contener Species_ID y Original_Taxon_ID.")
}

# -----------------------------
# 5. Extraer matriz de conteos
# -----------------------------

cat("[2] Extracting count matrix\n")

otu <- as(phyloseq::otu_table(ps), "matrix")

if (!phyloseq::taxa_are_rows(ps)) {
  otu <- t(otu)
}

storage.mode(otu) <- "numeric"
otu[is.na(otu)] <- 0

rownames(otu) <- as.character(rownames(otu))
colnames(otu) <- as.character(colnames(otu))

cat("Count matrix taxa x samples:", nrow(otu), "x", ncol(otu), "\n")
cat("First OTU rownames:\n")
print(head(rownames(otu)))
cat("\n")

# -----------------------------
# 6. Mapear especies bimodales
# -----------------------------

cat("[3] Mapping bimodal Species_ID to original phyloseq taxa\n")

bimodal_map <- bimodal_tbl %>%
  select(Species_ID, everything()) %>%
  left_join(
    taxmap %>%
      select(Species_ID, Original_Taxon_ID),
    by = "Species_ID",
    suffix = c(".bimodal", ".taxmap")
  )

if ("Original_Taxon_ID.bimodal" %in% colnames(bimodal_map) &&
    "Original_Taxon_ID.taxmap" %in% colnames(bimodal_map)) {
  bimodal_map <- bimodal_map %>%
    mutate(
      Original_Taxon_ID = dplyr::coalesce(
        as.character(Original_Taxon_ID.bimodal),
        as.character(Original_Taxon_ID.taxmap)
      )
    ) %>%
    select(-Original_Taxon_ID.bimodal, -Original_Taxon_ID.taxmap)
}

if (!"Original_Taxon_ID" %in% colnames(bimodal_map)) {
  stop("No se pudo construir Original_Taxon_ID después del mapeo.")
}

bimodal_map <- bimodal_map %>%
  mutate(
    Species_ID = as.character(Species_ID),
    Original_Taxon_ID = as.character(Original_Taxon_ID),
    in_phyloseq_taxa_names = Original_Taxon_ID %in% as.character(phyloseq::taxa_names(ps)),
    in_otu_rownames = Original_Taxon_ID %in% rownames(otu)
  )

write_csv(bimodal_map, file.path(OUT_DIR, "002_bimodal_species_mapped_to_phyloseq.csv"))

candidate_taxa <- bimodal_map %>%
  filter(!is.na(Original_Taxon_ID)) %>%
  pull(Original_Taxon_ID) %>%
  unique() %>%
  as.character()

bimodal_taxa_original <- intersect(candidate_taxa, rownames(otu))
missing_after_intersect <- setdiff(candidate_taxa, rownames(otu))

write_csv(
  tibble(Original_Taxon_ID = missing_after_intersect),
  file.path(OUT_DIR, "002_bimodal_species_missing_from_otu_rownames.csv")
)

cat("Candidate bimodal taxa:", length(candidate_taxa), "\n")
cat("Mapped bimodal taxa in OTU rownames:", length(bimodal_taxa_original), "\n")
cat("Missing after intersect:", length(missing_after_intersect), "\n\n")

if (length(bimodal_taxa_original) < 10) {
  stop("Muy pocas especies bimodales pudieron mapearse a la matriz OTU.")
}

# -----------------------------
# 7. CLR universo completo y subset bimodal
# -----------------------------

cat("[4] CLR transform on full species universe, then subset bimodal taxa\n")

log_counts <- log(otu + PSEUDOCOUNT)
clr_all_taxa_by_samples <- sweep(log_counts, 2, colMeans(log_counts), FUN = "-")

rownames(clr_all_taxa_by_samples) <- rownames(otu)
colnames(clr_all_taxa_by_samples) <- colnames(otu)

bimodal_taxa_original <- intersect(bimodal_taxa_original, rownames(clr_all_taxa_by_samples))

if (length(bimodal_taxa_original) < 10) {
  stop("Después del CLR, muy pocos taxa bimodales coinciden con rownames.")
}

clr_bimodal_taxa_by_samples <- clr_all_taxa_by_samples[bimodal_taxa_original, , drop = FALSE]
clr_bimodal_samples_by_taxa <- t(clr_bimodal_taxa_by_samples)

write_csv(
  as.data.frame(clr_bimodal_taxa_by_samples) %>%
    rownames_to_column("Original_Taxon_ID"),
  file.path(OUT_DIR, "002_clr_bimodal_species_taxa_by_samples.csv")
)

cat("CLR bimodal matrix samples x taxa:",
    nrow(clr_bimodal_samples_by_taxa), "x",
    ncol(clr_bimodal_samples_by_taxa), "\n\n")

# -----------------------------
# 8. Metadata y alineación
# -----------------------------

cat("[5] Aligning metadata\n")

meta_ps <- as.data.frame(phyloseq::sample_data(ps), stringsAsFactors = FALSE) %>%
  rownames_to_column("sample_id_phyloseq")

meta_axis_raw <- tibble(sample_id = rownames(clr_bimodal_samples_by_taxa)) %>%
  left_join(master, by = "sample_id", suffix = c("", ".master")) %>%
  left_join(meta_ps, by = c("sample_id" = "sample_id_phyloseq"), suffix = c(".master", ".ps"))

write_csv(
  tibble(column_name = colnames(meta_axis_raw)),
  file.path(OUT_DIR, "002_meta_axis_raw_column_names.csv")
)

meta_axis <- meta_axis_raw %>%
  mutate(
    collapsed_stage_clean = coalesce_existing(
      .,
      c("collapsed_stage", "collapsed_stage.master", "collapsed_stage.ps")
    ),
    depth_cm_clean = coalesce_existing(
      .,
      c("depth_cm", "depth_cm.master", "depth_cm.ps", "Depth_cm", "depth")
    ),
    locality_clean = coalesce_existing(
      .,
      c("locality", "locality.master", "locality.ps", "Locality")
    )
  ) %>%
  mutate(
    collapsed_stage = as.character(collapsed_stage_clean),
    depth_cm = as.character(depth_cm_clean),
    locality = as.character(locality_clean)
  )

write_csv(meta_axis, file.path(OUT_DIR, "002_meta_axis_after_join.csv"))

missing_stage <- sum(is.na(meta_axis$collapsed_stage))
missing_depth <- sum(is.na(meta_axis$depth_cm))
missing_locality <- sum(is.na(meta_axis$locality))

cat("Samples missing collapsed_stage:", missing_stage, "\n")
cat("Samples missing depth_cm:", missing_depth, "\n")
cat("Samples missing locality:", missing_locality, "\n")

if (missing_stage > 0 || missing_depth > 0 || missing_locality > 0) {
  write_csv(
    meta_axis %>%
      filter(is.na(collapsed_stage) | is.na(depth_cm) | is.na(locality)) %>%
      select(sample_id, collapsed_stage, depth_cm, locality),
    file.path(OUT_DIR, "002_samples_missing_required_metadata.csv")
  )
}

keep_samples <- meta_axis %>%
  filter(!is.na(collapsed_stage), !is.na(depth_cm), !is.na(locality)) %>%
  pull(sample_id)

clr_use <- clr_bimodal_samples_by_taxa[keep_samples, , drop = FALSE]

meta_use <- meta_axis %>%
  filter(sample_id %in% keep_samples) %>%
  mutate(
    collapsed_stage = factor(collapsed_stage, levels = STAGE_LEVELS),
    depth_cm = factor(depth_cm),
    locality = factor(locality)
  ) %>%
  arrange(match(sample_id, rownames(clr_use)))

if (!all(meta_use$sample_id == rownames(clr_use))) {
  stop("Error de alineación entre metadata y matriz CLR.")
}

cat("Samples retained for axis:", nrow(meta_use), "\n")
cat("Stages:\n")
print(table(meta_use$collapsed_stage, useNA = "ifany"))
cat("Depths:\n")
print(table(meta_use$depth_cm, useNA = "ifany"))
cat("Localities:\n")
print(table(meta_use$locality, useNA = "ifany"))
cat("\n")

if (nrow(meta_use) < 10) {
  stop("Muy pocas muestras con metadata completa para construir el eje.")
}

# -----------------------------
# 9. PCA: Bimodal Microbiome Axis
# -----------------------------

cat("[6] PCA on conservative bimodal species CLR matrix\n")

pca <- prcomp(clr_use, center = TRUE, scale. = FALSE)

var_explained <- (pca$sdev^2) / sum(pca$sdev^2)

axis_tbl <- meta_use %>%
  mutate(
    BMA1_raw = pca$x[, 1],
    BMA2_raw = pca$x[, 2],
    BMA3_raw = pca$x[, 3]
  )

med_degraded <- median(axis_tbl$BMA1_raw[axis_tbl$collapsed_stage == "Degraded"], na.rm = TRUE)
med_preserved <- median(axis_tbl$BMA1_raw[axis_tbl$collapsed_stage == "Preserved"], na.rm = TRUE)

flip <- ifelse(
  is.finite(med_degraded) &&
    is.finite(med_preserved) &&
    med_preserved < med_degraded,
  -1,
  1
)

axis_tbl <- axis_tbl %>%
  mutate(
    BMA1 = flip * BMA1_raw,
    BMA2 = BMA2_raw,
    BMA3 = BMA3_raw
  )

taxmap_join <- safe_select_taxmap(taxmap)

loadings_tbl <- as.data.frame(pca$rotation[, 1:5, drop = FALSE]) %>%
  rownames_to_column("Original_Taxon_ID") %>%
  mutate(
    Original_Taxon_ID = as.character(Original_Taxon_ID),
    PC1 = flip * PC1
  ) %>%
  left_join(
    taxmap_join,
    by = "Original_Taxon_ID"
  ) %>%
  left_join(
    bimodal_tbl %>%
      mutate(Species_ID = as.character(Species_ID)) %>%
      select(
        Species_ID,
        any_of(c(
          "Bimodality_Class",
          "Stability_Score",
          "Freq_Bimodal_BC",
          "Freq_Bimodal_Dip",
          "Freq_Bimodal_Joint",
          "BC_observed",
          "Dip_P_observed"
        ))
      ),
    by = "Species_ID"
  ) %>%
  arrange(desc(abs(PC1)))

variance_tbl <- tibble(
  axis = paste0("PC", seq_along(var_explained)),
  variance_explained = var_explained,
  cumulative_variance = cumsum(var_explained)
)

write_csv(axis_tbl, file.path(OUT_DIR, "002_bimodal_microbiome_axis_by_sample.csv"))
write_csv(loadings_tbl, file.path(OUT_DIR, "002_bimodal_species_axis_loadings.csv"))
write_csv(variance_tbl, file.path(OUT_DIR, "002_pca_variance_explained.csv"))

cat("BMA1 variance explained:", round(var_explained[1], 4), "\n")
cat("BMA2 variance explained:", round(var_explained[2], 4), "\n")
cat("BMA1 orientation flip:", flip, "\n\n")

# -----------------------------
# 10. Estadística básica
# -----------------------------

cat("[7] Statistical tests\n")

kw_stage <- kruskal.test(BMA1 ~ collapsed_stage, data = axis_tbl)
kw_depth <- kruskal.test(BMA1 ~ depth_cm, data = axis_tbl)
kw_locality <- kruskal.test(BMA1 ~ locality, data = axis_tbl)

lm_fit <- lm(BMA1 ~ collapsed_stage + depth_cm + locality, data = axis_tbl)
lm_anova <- anova(lm_fit)

stats_tbl <- tibble(
  test = c(
    "Kruskal_BMA1_by_collapsed_stage",
    "Kruskal_BMA1_by_depth_cm",
    "Kruskal_BMA1_by_locality"
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
  )
)

lm_anova_tbl <- as.data.frame(lm_anova) %>%
  rownames_to_column("term")

write_csv(stats_tbl, file.path(OUT_DIR, "002_bma1_nonparametric_tests.csv"))
write_csv(lm_anova_tbl, file.path(OUT_DIR, "002_bma1_lm_anova.csv"))

cat("Kruskal BMA1 ~ collapsed_stage p:", kw_stage$p.value, "\n")
cat("Kruskal BMA1 ~ depth_cm p:", kw_depth$p.value, "\n")
cat("Kruskal BMA1 ~ locality p:", kw_locality$p.value, "\n\n")

pairwise_stage <- pairwise.wilcox.test(
  x = axis_tbl$BMA1,
  g = axis_tbl$collapsed_stage,
  p.adjust.method = "BH",
  exact = FALSE
)

pairwise_stage_tbl <- as.data.frame(as.table(pairwise_stage$p.value)) %>%
  filter(!is.na(Freq)) %>%
  rename(group1 = Var1, group2 = Var2, p_adj_BH = Freq)

write_csv(pairwise_stage_tbl, file.path(OUT_DIR, "002_bma1_pairwise_wilcox_stage.csv"))

# -----------------------------
# 11. PERMANOVA opcional
# -----------------------------

if (requireNamespace("vegan", quietly = TRUE)) {
  cat("[8] PERMANOVA on Aitchison distance of bimodal species\n")

  dist_bimodal <- dist(clr_use, method = "euclidean")

  adonis_fit <- vegan::adonis2(
    dist_bimodal ~ collapsed_stage + depth_cm + locality,
    data = axis_tbl,
    permutations = 999
  )

  adonis_tbl <- as.data.frame(adonis_fit) %>%
    rownames_to_column("term")

  write_csv(adonis_tbl, file.path(OUT_DIR, "002_adonis2_bimodal_species_aitchison.csv"))

  print(adonis_tbl)
  cat("\n")
} else {
  cat("[8] vegan not installed; skipping PERMANOVA\n\n")
}

# -----------------------------
# 12. Figuras
# -----------------------------

cat("[9] Plotting\n")

p_axis_stage <- ggplot(axis_tbl, aes(x = collapsed_stage, y = BMA1)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.65) +
  geom_jitter(aes(shape = depth_cm), width = 0.14, size = 2.8, alpha = 0.85) +
  facet_wrap(~ depth_cm, scales = "free_y") +
  theme_bw(base_size = 13) +
  labs(
    title = "Bimodal Microbiome Axis 1 across restoration stages",
    subtitle = "Axis built from conservative bimodal species",
    x = "Collapsed restoration stage",
    y = "BMA1"
  )

ggsave(file.path(PLOT_DIR, "002_BMA1_by_stage_and_depth.png"),
       p_axis_stage, width = 10, height = 6, dpi = 300)
ggsave(file.path(PLOT_DIR, "002_BMA1_by_stage_and_depth.pdf"),
       p_axis_stage, width = 10, height = 6)

p_axis_locality <- ggplot(axis_tbl, aes(x = collapsed_stage, y = BMA1)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.65) +
  geom_jitter(aes(shape = depth_cm), width = 0.14, size = 2.6, alpha = 0.85) +
  facet_wrap(~ locality, scales = "free_y") +
  theme_bw(base_size = 13) +
  labs(
    title = "BMA1 across restoration stages by locality",
    subtitle = "Testing whether the bimodal axis is geographically consistent",
    x = "Collapsed restoration stage",
    y = "BMA1"
  )

ggsave(file.path(PLOT_DIR, "002_BMA1_by_stage_and_locality.png"),
       p_axis_locality, width = 11, height = 6, dpi = 300)
ggsave(file.path(PLOT_DIR, "002_BMA1_by_stage_and_locality.pdf"),
       p_axis_locality, width = 11, height = 6)

p_pca <- ggplot(axis_tbl, aes(x = BMA1, y = BMA2)) +
  geom_point(aes(shape = depth_cm, color = collapsed_stage), size = 3.5, alpha = 0.9) +
  theme_bw(base_size = 13) +
  labs(
    title = "Bimodal species PCA",
    subtitle = paste0(
      "BMA1 = ", round(100 * var_explained[1], 1),
      "%; BMA2 = ", round(100 * var_explained[2], 1), "%"
    ),
    x = "BMA1",
    y = "BMA2",
    color = "Stage",
    shape = "Depth"
  )

if (all(table(axis_tbl$collapsed_stage) >= 3)) {
  p_pca <- p_pca +
    stat_ellipse(aes(color = collapsed_stage), type = "norm",
                 linewidth = 0.7, na.rm = TRUE)
}

ggsave(file.path(PLOT_DIR, "002_BMA1_BMA2_pca.png"),
       p_pca, width = 8, height = 6, dpi = 300)
ggsave(file.path(PLOT_DIR, "002_BMA1_BMA2_pca.pdf"),
       p_pca, width = 8, height = 6)

p_density <- ggplot(axis_tbl, aes(x = BMA1, fill = collapsed_stage)) +
  geom_density(alpha = 0.4, na.rm = TRUE) +
  facet_wrap(~ depth_cm, scales = "free_y") +
  theme_bw(base_size = 13) +
  labs(
    title = "Density of BMA1 by depth",
    subtitle = "Intermediate broadening would support transition-zone instability",
    x = "BMA1",
    y = "Density",
    fill = "Stage"
  )

ggsave(file.path(PLOT_DIR, "002_BMA1_density_by_depth.png"),
       p_density, width = 10, height = 6, dpi = 300)
ggsave(file.path(PLOT_DIR, "002_BMA1_density_by_depth.pdf"),
       p_density, width = 10, height = 6)

top_loadings <- bind_rows(
  loadings_tbl %>%
    arrange(desc(PC1)) %>%
    slice_head(n = 25) %>%
    mutate(direction = "positive"),
  loadings_tbl %>%
    arrange(PC1) %>%
    slice_head(n = 25) %>%
    mutate(direction = "negative")
) %>%
  mutate(
    label = ifelse(
      "Clean_Name" %in% colnames(.) & !is.na(Clean_Name) & Clean_Name != "",
      Clean_Name,
      Original_Taxon_ID
    ),
    label = make.unique(as.character(label)),
    label = factor(label, levels = unique(label[order(PC1)]))
  )

write_csv(top_loadings, file.path(OUT_DIR, "002_top_BMA1_species_loadings.csv"))

p_load <- ggplot(top_loadings, aes(x = label, y = PC1, fill = direction)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 10) +
  labs(
    title = "Top species loadings on BMA1",
    x = "Species",
    y = "PC1 loading"
  )

ggsave(file.path(PLOT_DIR, "002_top_BMA1_species_loadings.png"),
       p_load, width = 10, height = 12, dpi = 300)
ggsave(file.path(PLOT_DIR, "002_top_BMA1_species_loadings.pdf"),
       p_load, width = 10, height = 12)

# -----------------------------
# 13. Resumen
# -----------------------------

summary_tbl <- tibble(
  item = c(
    "script",
    "timestamp",
    "phyloseq_input",
    "bimodal_input",
    "master_input",
    "n_samples_phyloseq",
    "n_taxa_phyloseq",
    "n_bimodal_species_input",
    "n_bimodal_species_mapped",
    "n_bimodal_species_missing",
    "n_samples_axis",
    "n_bimodal_species_axis",
    "BMA1_variance_explained",
    "BMA2_variance_explained",
    "BMA1_orientation_flip",
    "kruskal_stage_p",
    "kruskal_depth_p",
    "kruskal_locality_p",
    "output_dir"
  ),
  value = as.character(c(
    SCRIPT_ID,
    TIMESTAMP,
    PS_RDS,
    BIMODAL_CSV,
    MASTER_CSV,
    phyloseq::nsamples(ps),
    phyloseq::ntaxa(ps),
    nrow(bimodal_tbl),
    length(bimodal_taxa_original),
    length(missing_after_intersect),
    nrow(axis_tbl),
    ncol(clr_use),
    var_explained[1],
    var_explained[2],
    flip,
    kw_stage$p.value,
    kw_depth$p.value,
    kw_locality$p.value,
    OUT_DIR
  ))
)

write_csv(summary_tbl, file.path(OUT_DIR, "002_run_summary.csv"))

cat("============================================================\n")
cat("[DONE] Bimodal species axis completed\n")
cat("Output directory:\n")
cat(OUT_DIR, "\n\n")
cat("Main outputs:\n")
cat(" - 002_bimodal_microbiome_axis_by_sample.csv\n")
cat(" - 002_bimodal_species_axis_loadings.csv\n")
cat(" - 002_top_BMA1_species_loadings.csv\n")
cat(" - 002_bma1_nonparametric_tests.csv\n")
cat(" - 002_bma1_lm_anova.csv\n")
cat(" - plots/002_BMA1_by_stage_and_depth.png\n")
cat(" - plots/002_BMA1_BMA2_pca.png\n")
cat(" - plots/002_BMA1_density_by_depth.png\n")
cat("============================================================\n")

sink()
