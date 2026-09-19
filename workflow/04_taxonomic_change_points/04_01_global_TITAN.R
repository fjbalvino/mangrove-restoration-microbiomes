#!/usr/bin/env Rscript
# ============================================================
# 04_01_global_TITAN.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Global profile-mean TITAN2. Archived MHI branch uses pre-003b values; five other axes match corrected metadata.
# Inputs (source expressions; complete list in docs/contracts/04_01_global_TITAN.json):
#   df <- readr::read_csv(
#   trimws(readLines(f, warn = FALSE, n = 1L))
#   ps <- readRDS(PS_RDS)
#   meta_raw <- readr::read_csv(META_CSV, show_col_types = FALSE)
# Outputs (source expressions; complete list in contract):
#   ggsave(paste0(out_base, ".png"), p, width = 8, height = 8, dpi = 600, bg = "white")
#   ggsave(paste0(out_base, ".pdf"), p, width = 8, height = 8, device = cairo_pdf, bg = "white")
#   writeLines(
#   readr::write_csv(
#   readr::write_csv(profile_audit, file.path(DIR_TABLES, "401_auditoria_perfiles.csv"))
#   readr::write_csv(meta, file.path(DIR_TABLES, "401_metadata_used.csv"))
#   readr::write_csv(taxa_used, file.path(DIR_TABLES, "401_nonviral_bimodal_taxa_used.csv"))
#   readr::write_csv(sample_counts, file.path(DIR_TABLES, "401_sample_counts.csv"))
#   saveRDS(
# Algorithmic provenance:
# TITAN2 on profile means of Hellinger abundances, per environmental gradient.
#   Baker & King (2010), doi:10.1111/j.2041-210X.2009.00007.x.
# Source SHA-256: 511ae14dbcf5cf115de939168a89ca0e9ad4e9e9a3188b5ffc5d25d37d618b8b
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


# ==============================================================================
# 401_TITAN2_global_taxa_bimodales_gradientes_ambientales.R
#
# Objetivo:
#   Detectar umbrales comunitarios coordinados de taxa bimodales candidatos
#   a lo largo de gradientes ambientales continuos usando TITAN2.
#
# Comunidad:
#   - Universo fijo de 1031 taxa candidatos del script 202
#   - Conserva los 1031 taxa; estado viral como anotacion
#
# Gradientes principales:
#   - ECI_PC1
#   - nutrients_redox_PC1
#   - physicochemical_PC1
#   - moisture_stress_PC1
#   - water_inundation_PC1
#   - vegetation_landscape_PC1
#
# Outputs:
#   tables/
#     401_titan_summary_by_gradient.csv
#     401_titan_indval_<gradient>.csv
#     401_titan_sumz_<gradient>.csv
#     401_taxa_matrix_used.csv
#     401_metadata_used.csv
#
#   plots/
#     401_titan_sumz_<gradient>.png/pdf
#     401_titan_indicator_taxa_<gradient>.png/pdf
#
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(TITAN2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(scales)
})

# ==============================================================================
# Helpers
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

get_integer_arg <- function(flag, default) {
  val <- get_arg(flag, default = NA_character_)
  if (is.na(val) || is.null(val)) return(default)
  as.integer(val)
}

get_numeric_arg <- function(flag, default) {
  val <- get_arg(flag, default = NA_character_)
  if (is.na(val) || is.null(val)) return(default)
  as.numeric(val)
}

stop_if_missing_file <- function(path, label = "file") {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
    stop(sprintf("Missing %s: %s", label, path), call. = FALSE)
  }
}

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

clean_colnames <- function(x) make.names(x, unique = TRUE)

normalize_depth <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\.0$", "", x)
  x
}

clean_tax_string <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("^[a-zA-Z]__", "", x)
  x <- gsub("^.__", "", x)
  x <- gsub("_", " ", x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "unclassified", "Unclassified", "unknown", "Unknown")] <- NA_character_
  x
}

make_tax_label <- function(tax_df, tax_id_col = "tax_id") {
  ranks <- names(tax_df)

  species_col <- intersect(c("Species", "species", "S", "Rank_S"), ranks)[1]
  genus_col   <- intersect(c("Genus", "genus", "G", "Rank_G"), ranks)[1]
  family_col  <- intersect(c("Family", "family", "F", "Rank_F"), ranks)[1]
  order_col   <- intersect(c("Order", "order", "O", "Rank_O"), ranks)[1]

  species <- if (!is.na(species_col)) clean_tax_string(tax_df[[species_col]]) else rep(NA_character_, nrow(tax_df))
  genus   <- if (!is.na(genus_col)) clean_tax_string(tax_df[[genus_col]]) else rep(NA_character_, nrow(tax_df))
  family  <- if (!is.na(family_col)) clean_tax_string(tax_df[[family_col]]) else rep(NA_character_, nrow(tax_df))
  order   <- if (!is.na(order_col)) clean_tax_string(tax_df[[order_col]]) else rep(NA_character_, nrow(tax_df))

  label <- species

  bad_species <- is.na(label) |
    label == "" |
    grepl("^sp\\.?$", label, ignore.case = TRUE) |
    grepl("uncultured|unclassified|metagenome|environmental sample", label, ignore.case = TRUE)

  label[bad_species & !is.na(genus)] <- paste0(genus[bad_species & !is.na(genus)], " sp.")

  bad_label <- is.na(label) | label == ""
  label[bad_label & !is.na(family)] <- paste0(family[bad_label & !is.na(family)], " taxon")

  bad_label <- is.na(label) | label == ""
  label[bad_label & !is.na(order)] <- paste0(order[bad_label & !is.na(order)], " taxon")

  bad_label <- is.na(label) | label == ""
  label[bad_label] <- as.character(tax_df[[tax_id_col]][bad_label])

  make.unique(label)
}

is_viral_taxon <- function(tax_df) {
  viral_regex <- paste(
    c(
      "virus", "viruses", "viral", "viridae", "virinae", "virophage", "phage",
      "caudoviricetes", "duplodnaviria", "monodnaviria", "riboviria",
      "varidnaviria", "uroviricota", "nucleocytoviricota",
      "heunggongvirae", "loebvirae", "trapavirae", "shotokuvirae",
      "orthornavirae", "pararnavirae", "emperorvirus", "bellamyvirus",
      "stormageddonvirus", "nilusvirus", "zhoulongquanvirus",
      "mimivirus", "tupanvirus", "mollivirus", "phikzvirus",
      "myovirus", "podovirus", "siphovirus"
    ),
    collapse = "|"
  )

  tax_chr <- tax_df %>%
    mutate(across(everything(), as.character))

  apply(tax_chr, 1, function(z) {
    z <- paste(z, collapse = " ")
    grepl(viral_regex, z, ignore.case = TRUE)
  })
}

get_otu_samples_by_taxa <- function(ps) {
  otu <- phyloseq::otu_table(ps)
  mat <- as(otu, "matrix")
  if (phyloseq::taxa_are_rows(ps)) mat <- t(mat)
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  mat
}

read_bimodal_taxa_from_202 <- function(out202) {
  f <- file.path(
    out202, "tables",
    "202_CLR_raw1031_muestras_por_taxa.csv"
  )
  stop_if_missing_file(f, "matriz raw1031")
  df <- readr::read_csv(
    f, n_max = 0, show_col_types = FALSE,
    name_repair = "minimal"
  )
  if (sum(names(df) == "sample_id") != 1L) {
    stop("La matriz 202 requiere una columna sample_id.", call. = FALSE)
  }
  ids <- names(df)[names(df) != "sample_id"]
  if (length(ids) != 1031L || anyNA(ids) ||
      any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("202 no contiene exactamente 1031 IDs unicos.", call. = FALSE)
  }
  ids
}

match_bimodal_taxa_to_phyloseq <- function(bimodal_taxa, ps_taxa) {
  missing <- setdiff(bimodal_taxa, ps_taxa)
  if (anyDuplicated(ps_taxa) || length(missing)) {
    stop(
      "Matching exacto fallido. IDs ausentes: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  bimodal_taxa
}

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

extract_titan_table <- function(titan_obj) {
  out <- tryCatch({
    x <- as.data.frame(TITAN2::titan.out(titan_obj))
    x %>% rownames_to_column("taxon")
  }, error = function(e) {
    tibble(error = e$message)
  })
  out
}

extract_sumz_table <- function(titan_obj) {
  out <- tryCatch({
    sx <- as.data.frame(titan_obj$sumz)
    sx %>% rownames_to_column("row")
  }, error = function(e) {
    tibble(error = e$message)
  })
  out
}

safe_titan <- function(env, taxa_mat, nboot, ncpus, min_splt, num_permutations) {
  tryCatch({
    TITAN2::titan(
      env = env,
      txa = taxa_mat,
      minSplt = min_splt,
      numPerm = num_permutations,
      boot = TRUE,
      nBoot = nboot,
      imax = FALSE,
      ivTot = FALSE,
      pur.cut = 0.95,
      rel.cut = 0.95,
      ncpus = ncpus,
      memory = FALSE
    )
  }, error = function(e) {
    msg("TITAN failed: ", e$message)
    NULL
  })
}

plot_sumz_curves <- function(titan_obj, gradient_name, out_base) {
  plot_one <- function(device_fun, file, width = 8, height = 5.5, res = NULL) {
    if (is.null(res)) {
      device_fun(file, width = width, height = height)
    } else {
      device_fun(file, width = width, height = height, units = "in", res = res)
    }

    ok <- TRUE

    # Some TITAN2 versions print sumz quantiles without opening a plot.
    # We force a graphics page first and then try plot_sumz safely.
    try(graphics::plot.new(), silent = TRUE)

    out <- try(
      TITAN2::plot_sumz(titan_obj, filter = TRUE),
      silent = TRUE
    )

    if (inherits(out, "try-error")) {
      ok <- FALSE
      graphics::plot.new()
      graphics::text(
        0.5, 0.6,
        labels = paste0("TITAN2 plot_sumz failed
", gradient_name),
        cex = 1.1
      )
      graphics::text(
        0.5, 0.45,
        labels = as.character(out),
        cex = 0.7
      )
    }

    try(
      graphics::title(main = paste0("TITAN2 community threshold: ", gradient_name)),
      silent = TRUE
    )

    dev.off()
    invisible(ok)
  }

  plot_one(
    grDevices::png,
    paste0(out_base, ".png"),
    width = 8,
    height = 5.5,
    res = 600
  )

  plot_one(
    grDevices::pdf,
    paste0(out_base, ".pdf"),
    width = 8,
    height = 5.5,
    res = NULL
  )
}

plot_indicator_taxa <- function(ind_tbl, gradient_name, out_base, top_n = 40) {
  if (!all(c("taxon") %in% names(ind_tbl))) return(invisible(NULL))

  # TITAN2 column names vary slightly across versions.
  cp_col <- intersect(c("obs.cp", "cp", "env.cp", "env.cp.obs"), names(ind_tbl))[1]
  z_col  <- intersect(c("zscore", "obs.z", "z", "maxgrp"), names(ind_tbl))[1]
  dir_col <- intersect(c("maxgrp", "grp", "direction"), names(ind_tbl))[1]
  pur_col <- intersect(c("purity", "pur"), names(ind_tbl))[1]
  rel_col <- intersect(c("reliability", "rel"), names(ind_tbl))[1]

  if (is.na(cp_col)) return(invisible(NULL))

  df <- ind_tbl %>%
    mutate(
      change_point = suppressWarnings(as.numeric(.data[[cp_col]])),
      z_value = if (!is.na(z_col)) suppressWarnings(as.numeric(.data[[z_col]])) else NA_real_,
      direction_raw = if (!is.na(dir_col)) as.character(.data[[dir_col]]) else NA_character_,
      purity = if (!is.na(pur_col)) suppressWarnings(as.numeric(.data[[pur_col]])) else NA_real_,
      reliability = if (!is.na(rel_col)) suppressWarnings(as.numeric(.data[[rel_col]])) else NA_real_
    ) %>%
    filter(is.finite(change_point)) %>%
    mutate(
      direction = case_when(
        direction_raw %in% c("1", "z+", "+", "pos", "positive") ~ "z+",
        direction_raw %in% c("2", "z-", "-", "neg", "negative") ~ "z-",
        TRUE ~ direction_raw
      ),
      score_for_rank = ifelse(is.finite(z_value), abs(z_value), 0)
    ) %>%
    arrange(desc(score_for_rank)) %>%
    slice_head(n = top_n) %>%
    mutate(
      taxon = factor(taxon, levels = rev(taxon))
    )

  p <- ggplot(df, aes(x = change_point, y = taxon, color = direction, size = score_for_rank)) +
    geom_point(alpha = 0.9) +
    labs(
      title = NULL,
      subtitle = paste0("Top indicator taxa along ", gradient_name),
      x = paste0(gradient_name, " change point"),
      y = NULL,
      color = "Direction",
      size = "|z|"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.45),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      legend.position = "right"
    )

  ggsave(paste0(out_base, ".png"), p, width = 8, height = 8, dpi = 600, bg = "white")
  ggsave(paste0(out_base, ".pdf"), p, width = 8, height = 8, device = cairo_pdf, bg = "white")
}

# ==============================================================================
# Parameters
# ==============================================================================

RESULTS_DIR <- "/home/fjbalvino/Tipping_points/resultados_finales"

latest_output <- function(root, id) {
  f <- file.path(root, paste0("LATEST_", id, ".txt"))
  if (!file.exists(f)) return(NA_character_)
  trimws(readLines(f, warn = FALSE, n = 1L))
}

RUN001 <- latest_output(RESULTS_DIR, "001_auditar_estructura_metadata")
DEFAULT_PS <- if (!is.na(RUN001)) file.path(RUN001, "rds", "001_phyloseq_canon_51.rds") else NA_character_
DEFAULT_OUT202 <- latest_output(RESULTS_DIR, "202_construir_matriz_CLR_y_eje_bimodal")
RUN101 <- latest_output(RESULTS_DIR, "101_construir_ejes_ambientales_por_bloque")
RUN003 <- latest_output(RESULTS_DIR, "003_integrar_ejes_ECI_HI_en_metadata")
DEFAULT_META <- if (!is.na(RUN003)) {
  file.path(RUN003, "tables", "003_metadata_integrada_canon_51.csv")
} else if (!is.na(RUN101)) {
  file.path(RUN101, "tables", "101_metadata_with_env_axes.csv")
} else NA_character_

PS_RDS <- get_arg("--ps", DEFAULT_PS)
OUT202 <- get_arg("--out202", DEFAULT_OUT202)
META_CSV <- get_arg("--meta_csv", DEFAULT_META)
OUT_ROOT <- get_arg("--out_root", RESULTS_DIR)

GRADIENTS_ARG <- get_arg(
  "--gradients",
  "MHI_local,nutrients_redox_PC1,physicochemical_PC1,moisture_stress_PC1,water_inundation_PC1,vegetation_landscape_PC1"
)

MIN_PREV_N <- get_integer_arg("--min_prev_n", 3)
MIN_TOTAL <- get_numeric_arg("--min_total", 1)
N_BOOT <- get_integer_arg("--n_boot", 500)
N_PERM <- get_integer_arg("--n_perm", 999)
N_CPUS <- get_integer_arg("--ncpus", 1)
MIN_SPLT <- get_integer_arg("--min_splt", 5)
SEED <- get_integer_arg("--seed", 123)

GRADIENTS <- trimws(unlist(strsplit(GRADIENTS_ARG, ",")))
GRADIENTS <- GRADIENTS[nzchar(GRADIENTS)]

SCRIPT_BASENAME <- "401_TITAN2_global_taxa_bimodales_gradientes_ambientales"

OUT_DIR <- file.path(
  OUT_ROOT,
  paste0(SCRIPT_BASENAME, "_", timestamp_now())
)

DIR_TABLES <- file.path(OUT_DIR, "tables")
DIR_PLOTS <- file.path(OUT_DIR, "plots")
DIR_OBJECTS <- file.path(OUT_DIR, "objects")
DIR_LOGS <- file.path(OUT_DIR, "logs")

safe_mkdir(OUT_DIR)
safe_mkdir(DIR_TABLES)
safe_mkdir(DIR_PLOTS)
safe_mkdir(DIR_OBJECTS)
safe_mkdir(DIR_LOGS)

LOG_FILE <- file.path(OUT_DIR, "log.txt")
sink(LOG_FILE, split = TRUE)
on.exit(sink(), add = TRUE)

writeLines(
  OUT_DIR,
  file.path(OUT_ROOT, paste0("LATEST_", SCRIPT_BASENAME, ".txt"))
)

set.seed(SEED)

msg("Starting script 401")
msg("PS_RDS: ", PS_RDS)
msg("OUT202: ", OUT202)
msg("META_CSV: ", META_CSV)
msg("OUT_DIR: ", OUT_DIR)
msg("GRADIENTS: ", paste(GRADIENTS, collapse = ", "))
msg("MIN_PREV_N: ", MIN_PREV_N)
msg("MIN_TOTAL: ", MIN_TOTAL)
msg("N_BOOT: ", N_BOOT)
msg("N_PERM: ", N_PERM)
msg("N_CPUS: ", N_CPUS)
msg("MIN_SPLT: ", MIN_SPLT)
msg("SEED: ", SEED)

stop_if_missing_file(PS_RDS, "phyloseq rds")
stop_if_missing_file(META_CSV, "metadata csv")

# ==============================================================================
# Read data
# ==============================================================================

ps <- readRDS(PS_RDS)
otu <- get_otu_samples_by_taxa(ps)

if (!is.null(phyloseq::tax_table(ps, errorIfNULL = FALSE))) {
  tax <- as.data.frame(phyloseq::tax_table(ps), stringsAsFactors = FALSE) %>%
    rownames_to_column("tax_id")
} else {
  tax <- tibble(tax_id = colnames(otu))
}

tax$tax_label <- make_tax_label(tax, tax_id_col = "tax_id")
tax$is_viral <- is_viral_taxon(tax)

readr::write_csv(
  tax %>% select(tax_id, tax_label, is_viral, everything()),
  file.path(DIR_TABLES, "401_taxonomy_with_viral_flag.csv")
)

meta_raw <- readr::read_csv(META_CSV, show_col_types = FALSE)
names(meta_raw) <- clean_colnames(names(meta_raw))

sample_col <- intersect(
  c("sample_id", "SampleID", "sample", "id_metagenomics", "run_accession"),
  names(meta_raw)
)[1]

profile_col <- intersect(
  c("profile_id", "lat_block", "latitude_block", "perfil", "core_id"),
  names(meta_raw)
)[1]

if (is.na(sample_col)) {
  stop("No pude identificar sample_id en metadata.", call. = FALSE)
}
if (is.na(profile_col)) {
  stop("No pude identificar profile_id/lat_block; TITAN2 global requiere perfiles independientes.", call. = FALSE)
}

# ECI y MHI son indices distintos: no sustituir uno por otro.
if (!"ECI_PC1" %in% names(meta_raw)) {
  eci_aliases <- intersect(c("ECI.PC1", "eci_pc1"), names(meta_raw))
  if (length(eci_aliases) == 1L) {
    meta_raw$ECI_PC1 <- meta_raw[[eci_aliases]]
  }
}

missing_grad <- setdiff(GRADIENTS, names(meta_raw))
if (length(missing_grad) > 0) {
  msg("WARNING: missing gradients will be skipped: ", paste(missing_grad, collapse = ", "))
}

GRADIENTS <- intersect(GRADIENTS, names(meta_raw))

if (length(GRADIENTS) == 0) {
  stop("No valid gradients found in metadata.", call. = FALSE)
}

required_meta <- c("restoration4", "depth_cm", "locality", GRADIENTS)
missing_meta <- setdiff(required_meta, names(meta_raw))
if (length(missing_meta) > 0) {
  stop("Faltan columnas en metadata: ", paste(missing_meta, collapse = ", "), call. = FALSE)
}

STAGE_LEVELS <- c(
  "Degraded",
  "Early restoration",
  "Intermediate restoration",
  "Advanced restoration",
  "Conserved"
)

DEPTH_LEVELS <- c("5", "20", "40")

meta <- meta_raw %>%
  transmute(
    sample_id = as.character(.data[[sample_col]]),
    restoration4 = factor(as.character(restoration4), levels = STAGE_LEVELS, ordered = TRUE),
    depth_cm = factor(normalize_depth(depth_cm), levels = DEPTH_LEVELS, ordered = TRUE),
    locality = factor(as.character(locality)),
    profile_raw = as.character(.data[[profile_col]]),
    profile_id = if_else(
      str_detect(profile_raw, fixed(as.character(locality))),
      profile_raw,
      paste(as.character(locality), profile_raw, sep = "::")
    ),
    across(all_of(GRADIENTS), ~ suppressWarnings(as.numeric(.x)))
  ) %>%
  filter(
    !is.na(sample_id),
    !is.na(restoration4),
    !is.na(depth_cm),
    !is.na(locality)
  )

# ==============================================================================
# Taxa selection: fixed raw1031 universe
# ==============================================================================

bimodal_taxa_raw <- read_bimodal_taxa_from_202(OUT202)
bimodal_taxa <- match_bimodal_taxa_to_phyloseq(bimodal_taxa_raw, colnames(otu))

viral_taxa <- tax %>%
  filter(is_viral) %>%
  pull(tax_id)

viral_taxa <- intersect(viral_taxa, colnames(otu))

# Universo acordado: los 1031 taxa de 202, incluida su fraccion viral.
# La clasificacion viral se conserva como anotacion.
selected_taxa <- bimodal_taxa

if (length(selected_taxa) < 5) {
  stop("Muy pocos taxa candidatos seleccionados: ", length(selected_taxa), call. = FALSE)
}

common_samples <- intersect(rownames(otu), meta$sample_id)

if (length(common_samples) != 51L || !setequal(rownames(otu), meta$sample_id)) {
  stop("401 requiere coincidencia exacta de 51 sample_id entre phyloseq y metadata.", call. = FALSE)
}

meta <- meta %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

counts <- otu[meta$sample_id, selected_taxa, drop = FALSE]

prev <- colSums(counts > 0, na.rm = TRUE)
tot <- colSums(counts, na.rm = TRUE)

keep_taxa <- prev >= MIN_PREV_N & tot >= MIN_TOTAL
if (!all(keep_taxa) || ncol(counts) != 1031L ||
    !setequal(colnames(counts), bimodal_taxa_raw)) {
  stop(
    "Los filtros o el matching alteran el universo de 1031 taxa.",
    call. = FALSE
  )
}
counts <- counts[, keep_taxa, drop = FALSE]

if (ncol(counts) < 5) {
  stop("Muy pocos taxa tras filtros de prevalencia/abundancia: ", ncol(counts), call. = FALSE)
}

# TITAN2 global uses the independent profile as its analysis unit. Counts are
# converted to Hellinger abundance per sample and averaged across 5/20/40 cm.
profile_audit <- meta %>%
  group_by(profile_id) %>%
  summarise(
    n_samples = n(),
    depths = paste(sort(unique(as.character(depth_cm))), collapse = "|"),
    n_locality = n_distinct(locality),
    n_stage = n_distinct(restoration4),
    valid = n_samples == 3L && setequal(as.character(depth_cm), DEPTH_LEVELS) &&
      n_locality == 1L && n_stage == 1L,
    .groups = "drop"
  )
readr::write_csv(profile_audit, file.path(DIR_TABLES, "401_auditoria_perfiles.csv"))
if (nrow(profile_audit) != 17L || !all(profile_audit$valid)) {
  stop("401 requiere exactamente 17 perfiles completos 5/20/40.", call. = FALSE)
}

# Universo raw1031: validar los 17 perfiles sin excluir muestras.
if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
  stop("Counts contiene valores invalidos.", call. = FALSE)
}
if (nrow(meta) != 51L || anyDuplicated(meta$sample_id) ||
    dplyr::n_distinct(meta$profile_id) != 17L ||
    !all(table(meta$profile_id) == 3L)) {
  stop("Se requieren 51 muestras de 17 perfiles completos.", call. = FALSE)
}
if (any(rowSums(counts) <= 0)) {
  stop("Suma cero inesperada en raw1031; revisar 202.", call. = FALSE)
}

readr::write_csv(
  profile_audit,
  file.path(DIR_TABLES, "401_auditoria_perfiles_retenidos.csv")
)
readr::write_csv(
  tibble::tibble(
    universe = "raw1031",
    selection = "observed_joint_raw",
    n_taxa = ncol(counts),
    n_samples = nrow(counts),
    n_profiles = dplyr::n_distinct(meta$profile_id),
    n_excluded_profiles = 0L
  ),
  file.path(DIR_TABLES, "401_contrato_universo_raw1031.csv")
)
msg("Universo raw1031: 51 muestras; 17 perfiles; ninguna exclusion.")

lib_size <- rowSums(counts, na.rm = TRUE)
if (any(!is.finite(lib_size) | lib_size <= 0)) {
  stop("Hay muestras sin abundancia total en los taxa seleccionados.", call. = FALSE)
}
txa_sample <- sqrt(sweep(counts, 1, lib_size, "/"))
profile_order <- unique(meta$profile_id)
txa <- rowsum(txa_sample, group = meta$profile_id, reorder = FALSE)
txa <- sweep(txa, 1, as.numeric(table(factor(meta$profile_id, levels = rownames(txa)))), "/")

meta <- meta %>%
  group_by(profile_id) %>%
  summarise(
    restoration4 = first(restoration4),
    locality = first(locality),
    across(all_of(GRADIENTS), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(match(profile_id, rownames(txa)))
stopifnot(identical(as.character(meta$profile_id), rownames(txa)))

v <- apply(txa, 2, var, na.rm = TRUE)
if (ncol(txa) != 1031L || any(!is.finite(v) | v <= 0)) {
  stop("No se conservan 1031 taxa variables a nivel de perfil.", call. = FALSE)
}
txa <- txa[, is.finite(v) & v > 0, drop = FALSE]

if (ncol(txa) < 5) {
  stop("Muy pocos taxa variables para TITAN2: ", ncol(txa), call. = FALSE)
}

msg("Independent profiles retained: ", nrow(meta))
msg("Candidate taxa retained: ", ncol(txa))
msg("Viral taxa retained within raw1031: ", sum(colnames(txa) %in% viral_taxa))

readr::write_csv(meta, file.path(DIR_TABLES, "401_metadata_used.csv"))

taxa_used <- tibble(tax_id = colnames(txa)) %>%
  left_join(tax, by = "tax_id")

readr::write_csv(taxa_used, file.path(DIR_TABLES, "401_nonviral_bimodal_taxa_used.csv"))

readr::write_csv(
  as.data.frame(txa) %>% rownames_to_column("profile_id"),
  file.path(DIR_TABLES, "401_taxa_matrix_profile_mean_hellinger_used.csv")
)

sample_counts <- bind_rows(
  meta %>% count(locality, name = "n") %>% mutate(group = "locality", level = as.character(locality)) %>% select(group, level, n),
  meta %>% count(restoration4, name = "n") %>% mutate(group = "restoration4", level = as.character(restoration4)) %>% select(group, level, n),
  tibble(group = "analysis_unit", level = "independent_profile", n = nrow(meta))
)

readr::write_csv(sample_counts, file.path(DIR_TABLES, "401_sample_counts.csv"))

# ==============================================================================
# TITAN2 global by gradient
# ==============================================================================

summary_rows <- list()

for (grad in GRADIENTS) {
  msg("Running TITAN2 for gradient: ", grad)

  grad_safe <- sanitize_name(grad)

  use_idx <- is.finite(meta[[grad]])
  n_use <- sum(use_idx)

  if (n_use < max(12, 2 * MIN_SPLT + 2)) {
    msg("Skipping ", grad, ": too few valid samples (n=", n_use, ")")
    summary_rows[[grad]] <- tibble(
      gradient = grad,
      n_profiles = n_use,
      n_taxa = ncol(txa),
      status = "skipped_too_few_samples",
      n_indicators = NA_integer_,
      n_zplus = NA_integer_,
      n_zminus = NA_integer_
    )
    next
  }

  env <- meta[[grad]][use_idx]
  txa_use <- txa[use_idx, , drop = FALSE]

  # Sort by environmental gradient, as TITAN2 expects monotonic gradient order
  ord <- order(env)
  env <- env[ord]
  txa_use <- txa_use[ord, , drop = FALSE]

  titan_obj <- safe_titan(
    env = env,
    taxa_mat = txa_use,
    nboot = N_BOOT,
    ncpus = N_CPUS,
    min_splt = MIN_SPLT,
    num_permutations = N_PERM
  )

  if (is.null(titan_obj)) {
    summary_rows[[grad]] <- tibble(
      gradient = grad,
      n_profiles = n_use,
      n_taxa = ncol(txa_use),
      status = "failed",
      n_indicators = NA_integer_,
      n_zplus = NA_integer_,
      n_zminus = NA_integer_
    )
    next
  }

  saveRDS(
    titan_obj,
    file.path(DIR_OBJECTS, paste0("401_titan_", grad_safe, ".rds"))
  )

  ind_tbl <- extract_titan_table(titan_obj)
  sumz_tbl <- extract_sumz_table(titan_obj)

  readr::write_csv(
    ind_tbl,
    file.path(DIR_TABLES, paste0("401_titan_indval_", grad_safe, ".csv"))
  )

  readr::write_csv(
    sumz_tbl,
    file.path(DIR_TABLES, paste0("401_titan_sumz_", grad_safe, ".csv"))
  )

  # Plot native sum(z)
  plot_sumz_curves(
    titan_obj,
    gradient_name = grad,
    out_base = file.path(DIR_PLOTS, paste0("401_titan_sumz_", grad_safe))
  )

  # Plot indicator taxa
  plot_indicator_taxa(
    ind_tbl,
    gradient_name = grad,
    out_base = file.path(DIR_PLOTS, paste0("401_titan_indicator_taxa_", grad_safe)),
    top_n = 40
  )

  dir_col <- intersect(c("maxgrp", "grp", "direction"), names(ind_tbl))[1]
  pur_col <- intersect(c("purity", "pur"), names(ind_tbl))[1]
  rel_col <- intersect(c("reliability", "rel"), names(ind_tbl))[1]

  indicators <- ind_tbl

  if (!is.na(pur_col) && !is.na(rel_col)) {
    indicators <- indicators %>%
      mutate(
        purity_num = suppressWarnings(as.numeric(.data[[pur_col]])),
        reliability_num = suppressWarnings(as.numeric(.data[[rel_col]]))
      ) %>%
      filter(purity_num >= 0.95, reliability_num >= 0.95)
  }

  n_ind <- nrow(indicators)

  n_zplus <- NA_integer_
  n_zminus <- NA_integer_

  if (!is.na(dir_col) && nrow(indicators) > 0) {
    d <- as.character(indicators[[dir_col]])
    n_zplus <- sum(d %in% c("1", "z+", "+", "pos", "positive"), na.rm = TRUE)
    n_zminus <- sum(d %in% c("2", "z-", "-", "neg", "negative"), na.rm = TRUE)
  }

  summary_rows[[grad]] <- tibble(
    gradient = grad,
    n_profiles = n_use,
    n_taxa = ncol(txa_use),
    status = "ok",
    n_indicators = n_ind,
    n_zplus = n_zplus,
    n_zminus = n_zminus
  )
}

summary_tbl <- bind_rows(summary_rows)

readr::write_csv(
  summary_tbl,
  file.path(DIR_TABLES, "401_titan_summary_by_gradient.csv")
)

# ==============================================================================
# Run info
# ==============================================================================

run_info <- tibble(
  key = c(
    "script",
    "timestamp",
    "ps_rds",
    "out202",
    "meta_csv",
    "out_dir",
    "gradients",
    "n_profiles",
    "n_taxa",
    "n_viral_taxa_excluded",
    "analysis_unit",
    "abundance_transform",
    "min_prev_n",
    "min_total",
    "n_boot",
    "n_perm",
    "ncpus",
    "min_splt",
    "seed"
  ),
  value = c(
    SCRIPT_BASENAME,
    timestamp_now(),
    PS_RDS,
    OUT202,
    META_CSV,
    OUT_DIR,
    paste(GRADIENTS, collapse = ","),
    as.character(nrow(meta)),
    as.character(ncol(txa)),
    as.character(length(viral_taxa)),
    "profile_mean_across_5_20_40_cm",
    "sample_Hellinger_then_profile_mean",
    as.character(MIN_PREV_N),
    as.character(MIN_TOTAL),
    as.character(N_BOOT),
    as.character(N_PERM),
    as.character(N_CPUS),
    as.character(MIN_SPLT),
    as.character(SEED)
  )
)

readr::write_csv(
  run_info,
  file.path(DIR_TABLES, "401_run_info.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(DIR_LOGS, "401_sessionInfo.txt")
)

msg("Finished script 401 successfully")
msg("Output directory: ", OUT_DIR)
