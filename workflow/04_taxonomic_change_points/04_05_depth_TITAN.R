#!/usr/bin/env Rscript
# ============================================================
# 04_05_depth_TITAN.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Depth-specific TITAN2. MHI results require corrected 003b input and rerun.
# Inputs (source expressions; complete list in docs/contracts/04_05_depth_TITAN.json):
#   trimws(readLines(f, warn = FALSE, n = 1L))
#   x <- read.csv(f, nrows = 0, check.names = FALSE)
#   ps <- readRDS(ps_file)
#   meta <- readr::read_csv(meta_csv, show_col_types = FALSE, progress = FALSE)
# Outputs (source expressions; complete list in contract):
#   writeLines(out_dir, sub("\\.txt$", "_ATTEMPT.txt", latest_file))
#   readr::write_csv(as.data.frame(x), path, na = "")
#   saveRDS(titan_obj, rds_path)
#   ggsave(
#   writeLines(out_dir, latest_file)
# Algorithmic provenance:
# TITAN2 within each sediment depth, preserving the original Hellinger denominator.
#   Baker & King (2010), doi:10.1111/j.2041-210X.2009.00007.x.
# Source SHA-256: 787b81f57d0cc8c36fbdea4ab0b187c5f91883a21cfcddf603013c16a1eca5b1
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================

# RAW1031_405_V1

# ==============================================================================
# 405_TITAN2_depth_stratified_bimodal_nonviral.R
#
# Objective:
#   Run TITAN2 depth-stratified environmental threshold analyses for the
#   raw1031 bimodal community, including viral taxa.
#
# Major corrections relative to draft:
#   - Robust matching between OUT202 taxa and phyloseq taxa.
#   - Depth normalization: 5, 20, 40 instead of possible 5.0/20.0/40.0.
#   - Viral taxonomy is annotation only; no viral exclusion.
#   - Prevalence/total filtering is performed on raw counts, before transform.
#   - TITAN2 objects are parsed using obj$sumz.cp and obj$sppmax robustly.
#   - cp_plot uses obs_cp when available, otherwise cp_50.
#   - reliable90/reliable95 require both purity and reliability cutoffs.
#   - Outputs are versioned as 405_TITAN2_depth_stratified_bimodal_nonviral_YYYYMMDD_HHMMSS.
# ==============================================================================

options(stringsAsFactors = FALSE)
options(warn = 1)

# ------------------------------------------------------------------------------
# Package checks
# ------------------------------------------------------------------------------

required_pkgs <- c(
  "phyloseq",
  "TITAN2",
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "ggplot2",
  "stringr",
  "purrr",
  "scales"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them before running this script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(stringr)
  library(purrr)
  library(scales)
})

# ------------------------------------------------------------------------------
# Argument parser
# ------------------------------------------------------------------------------

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1L

  while (i <= length(args)) {
    key <- args[[i]]

    if (!startsWith(key, "--")) {
      stop("Unexpected argument without -- prefix: ", key, call. = FALSE)
    }

    key <- sub("^--", "", key)

    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }

  out
}

args <- parse_args()

arg_chr <- function(name, default = NULL, required = FALSE) {
  val <- args[[name]]

  if (is.null(val) || identical(val, TRUE) || is.na(val) || !nzchar(as.character(val))) {
    if (required) stop("Missing required argument: --", name, call. = FALSE)
    return(default)
  }

  as.character(val)
}

arg_num <- function(name, default = NULL, required = FALSE) {
  val <- arg_chr(name, default = NULL, required = required)
  if (is.null(val)) return(default)
  as.numeric(val)
}

arg_int <- function(name, default = NULL, required = FALSE) {
  as.integer(arg_num(name, default = default, required = required))
}

arg_lgl <- function(name, default = FALSE) {
  val <- args[[name]]

  if (is.null(val)) return(default)
  if (is.logical(val)) return(val)

  val <- tolower(as.character(val))
  val %in% c("1", "true", "t", "yes", "y")
}

split_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) {
    character(0)
  } else {
    trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  }
}

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------

SCRIPT_ID <- "405"
SCRIPT_NAME <- "405_TITAN2_taxa_bimodales_por_profundidad"

results_root <- arg_chr(
  "results_root",
  default = "/home/fjbalvino/Tipping_points/resultados_finales"
)

latest_output <- function(root, id) {
  f <- file.path(root, paste0("LATEST_", id, ".txt"))
  if (!file.exists(f)) return(NA_character_)
  trimws(readLines(f, warn = FALSE, n = 1L))
}

ps_file <- arg_chr(
  "ps",
  default = {
    run001 <- latest_output(results_root, "001_auditar_estructura_metadata")
    if (!is.na(run001)) file.path(run001, "rds", "001_phyloseq_canon_51.rds") else NA_character_
  }
)

out202_dir <- arg_chr(
  "out202",
  default = latest_output(results_root, "202_construir_matriz_CLR_y_eje_bimodal")
)

run101 <- latest_output(results_root, "101_construir_ejes_ambientales_por_bloque")
run003 <- latest_output(results_root, "003_integrar_ejes_ECI_HI_en_metadata")
meta_csv <- arg_chr(
  "meta_csv",
  default = if (!is.na(run003)) {
    file.path(run003, "tables", "003_metadata_integrada_canon_51.csv")
  } else if (!is.na(run101)) {
    file.path(run101, "tables", "101_metadata_with_env_axes.csv")
  } else NA_character_
)

gradients <- split_csv(arg_chr(
  "gradients",
  default = "MHI_local,moisture_stress_PC1,nutrients_redox_PC1,physicochemical_PC1,vegetation_landscape_PC1,water_inundation_PC1"
))

depth_levels <- split_csv(arg_chr(
  "depth_levels",
  default = "5,20,40"
))

sample_id_col <- arg_chr("sample_id_col", default = "sample_id")
depth_col <- arg_chr("depth_col", default = "depth_cm")
condition_col <- arg_chr("condition_col", default = "restoration4")
locality_col <- arg_chr("locality_col", default = "locality")

min_prev_n <- arg_int("min_prev_n", default = 3L)
min_total <- arg_num("min_total", default = 1)
min_splt <- arg_int("min_splt", default = 5L)
n_boot <- arg_int("n_boot", default = 500L)
n_perm <- arg_int("n_perm", default = 999L)
ncpus <- arg_int("ncpus", default = 1L)
seed <- arg_int("seed", default = 123L)

auto_relax <- arg_lgl("auto_relax", default = FALSE)
relaxed_min_prev_n <- arg_int("relaxed_min_prev_n", default = 3L)
relaxed_min_splt <- arg_int("relaxed_min_splt", default = 5L)

abundance_mode <- arg_chr(
  "abundance_mode",
  default = "hellinger"
)

min_samples_per_depth <- arg_int("min_samples_per_depth", default = 10L)
min_gradient_unique <- arg_int("min_gradient_unique", default = 6L)

reliability_cut <- arg_num("reliability_cut", default = 0.95)
purity_cut <- arg_num("purity_cut", default = 0.95)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(results_root, paste0(SCRIPT_NAME, "_", timestamp))

tables_dir <- file.path(out_dir, "tables")
plots_dir <- file.path(out_dir, "plots_nature")
rds_dir <- file.path(out_dir, "rds")
logs_dir <- file.path(out_dir, "logs")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

latest_file <- file.path(results_root, paste0("LATEST_", SCRIPT_NAME, ".txt"))
writeLines(out_dir, sub("\\.txt$", "_ATTEMPT.txt", latest_file))

log_file <- file.path(logs_dir, "405_run.log")
sink(log_file, split = TRUE)
on.exit({
  sink()
}, add = TRUE)

cat("===== 405 TITAN2 depth-stratified raw1031 analysis =====\n")
cat("Start time: ", as.character(Sys.time()), "\n", sep = "")
cat("Output directory: ", out_dir, "\n", sep = "")

# ------------------------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------------------------

write_csv_safe <- function(x, path) {
  readr::write_csv(as.data.frame(x), path, na = "")
  invisible(path)
}

clean_taxon_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x
}

normalize_depth <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\.0$", "", x)
  x
}

clean_tax_string <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("^D_[0-9]+__", "", x)
  x <- gsub("^[a-zA-Z]__", "", x)
  x <- gsub("^.__", "", x)
  x <- gsub("_", " ", x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "unclassified", "Unclassified", "unknown", "Unknown")] <- NA_character_
  x
}

make_full_taxonomy <- function(df) {
  tax_cols <- intersect(
    c(
      "Kingdom", "Superkingdom", "Domain",
      "Phylum", "Class", "Order", "Family", "Genus", "Species",
      "kingdom", "superkingdom", "domain",
      "phylum", "class", "order", "family", "genus", "species",
      "D", "P", "C", "O", "F", "G", "S",
      "Rank_D", "Rank_P", "Rank_C", "Rank_O", "Rank_F", "Rank_G", "Rank_S"
    ),
    names(df)
  )

  if (length(tax_cols) == 0) {
    return(rep(NA_character_, nrow(df)))
  }

  apply(df[, tax_cols, drop = FALSE], 1, function(z) {
    z <- clean_tax_string(z)
    z <- z[!is.na(z) & nzchar(z)]
    if (length(z) == 0) NA_character_ else paste(z, collapse = "; ")
  })
}

make_tax_label <- function(df, tax_id_col = "taxon_id") {
  ranks <- names(df)

  species_col <- intersect(c("Species", "species", "S", "Rank_S"), ranks)[1]
  genus_col   <- intersect(c("Genus", "genus", "G", "Rank_G"), ranks)[1]
  family_col  <- intersect(c("Family", "family", "F", "Rank_F"), ranks)[1]
  order_col   <- intersect(c("Order", "order", "O", "Rank_O"), ranks)[1]
  class_col   <- intersect(c("Class", "class", "C", "Rank_C"), ranks)[1]
  phylum_col  <- intersect(c("Phylum", "phylum", "P", "Rank_P"), ranks)[1]

  species <- if (!is.na(species_col)) clean_tax_string(df[[species_col]]) else rep(NA_character_, nrow(df))
  genus   <- if (!is.na(genus_col)) clean_tax_string(df[[genus_col]]) else rep(NA_character_, nrow(df))
  family  <- if (!is.na(family_col)) clean_tax_string(df[[family_col]]) else rep(NA_character_, nrow(df))
  order   <- if (!is.na(order_col)) clean_tax_string(df[[order_col]]) else rep(NA_character_, nrow(df))
  class_  <- if (!is.na(class_col)) clean_tax_string(df[[class_col]]) else rep(NA_character_, nrow(df))
  phylum  <- if (!is.na(phylum_col)) clean_tax_string(df[[phylum_col]]) else rep(NA_character_, nrow(df))

  label <- species

  bad <- is.na(label) |
    label == "" |
    grepl("^sp\\.?$", label, ignore.case = TRUE) |
    grepl("uncultured|unclassified|metagenome|environmental sample", label, ignore.case = TRUE)

  idx <- bad & !is.na(genus)
  label[idx] <- paste0(genus[idx], " sp.")

  bad <- is.na(label) | label == ""
  idx <- bad & !is.na(family)
  label[idx] <- paste0(family[idx], " taxon")

  bad <- is.na(label) | label == ""
  idx <- bad & !is.na(order)
  label[idx] <- paste0(order[idx], " taxon")

  bad <- is.na(label) | label == ""
  idx <- bad & !is.na(class_)
  label[idx] <- paste0(class_[idx], " taxon")

  bad <- is.na(label) | label == ""
  idx <- bad & !is.na(phylum)
  label[idx] <- paste0(phylum[idx], " taxon")

  bad <- is.na(label) | label == ""
  label[bad] <- as.character(df[[tax_id_col]][bad])

  label
}

is_viral_taxonomy <- function(x) {
  x <- tolower(paste(x, collapse = " "))
  grepl(
    paste(
      c(
        "virus", "viruses", "viral", "viridae", "virinae",
        "phage", "phages", "virophage",
        "caudoviricetes", "uroviricota", "duplodnaviria",
        "monodnaviria", "varidnaviria", "riboviria",
        "nucleocytoviricota", "heunggongvirae", "loebvirae",
        "trapavirae", "shotokuvirae", "orthornavirae",
        "pararnavirae", "myovirus", "podovirus", "siphovirus"
      ),
      collapse = "|"
    ),
    x
  )
}

read_bimodal_taxa_from_out202 <- function(out202_dir) {
  f <- file.path(out202_dir, "tables",
                 "202_CLR_raw1031_muestras_por_taxa.csv")
  if (!file.exists(f)) stop("Falta matriz raw1031: ", f)
  x <- read.csv(f, nrows = 0, check.names = FALSE)
  stopifnot(names(x)[1] == "sample_id")
  ids <- names(x)[-1]
  stopifnot(length(ids) == 1031L, !anyNA(ids),
            !anyDuplicated(ids), all(nzchar(ids)))
  ids
}

match_bimodal_taxa_to_phyloseq <- function(bimodal_taxa, ps_taxa) {
  stopifnot(length(bimodal_taxa) == 1031L,
            !anyDuplicated(bimodal_taxa),
            !anyDuplicated(ps_taxa))
  missing <- setdiff(bimodal_taxa, ps_taxa)
  if (length(missing)) {
    stop("IDs raw1031 ausentes en phyloseq: ",
         paste(missing, collapse = ", "))
  }
  bimodal_taxa
}

extract_phyloseq_counts_samples_by_taxa <- function(ps) {
  otu <- phyloseq::otu_table(ps)
  mat <- as(otu, "matrix")

  if (phyloseq::taxa_are_rows(ps)) {
    mat <- t(mat)
  }

  storage.mode(mat) <- "numeric"
  mat
}

make_taxonomy_map <- function(ps) {
  if (is.null(phyloseq::tax_table(ps, errorIfNULL = FALSE))) {
    out <- tibble(
      taxon_id = phyloseq::taxa_names(ps),
      full_taxonomy = NA_character_,
      tax_label = phyloseq::taxa_names(ps),
      is_viral = FALSE
    )
    return(out)
  }

  tx <- as.data.frame(as(phyloseq::tax_table(ps), "matrix"), stringsAsFactors = FALSE)
  tx$taxon_id <- rownames(tx)
  tx$full_taxonomy <- make_full_taxonomy(tx)
  tx$tax_label <- make_tax_label(tx, tax_id_col = "taxon_id")
  tx$is_viral <- apply(tx, 1, is_viral_taxonomy)

  tibble::as_tibble(tx)
}

prepare_abundance_matrix <- function(counts, mode = "hellinger") {
  mode <- tolower(mode)

  if (mode == "counts") {
    out <- counts
  } else if (mode == "hellinger") {
    rs <- rowSums(counts, na.rm = TRUE)
    rs[rs == 0] <- NA_real_
    out <- sqrt(sweep(counts, 1, rs, "/"))
    out[is.na(out)] <- 0
  } else if (mode == "relative") {
    rs <- rowSums(counts, na.rm = TRUE)
    rs[rs == 0] <- NA_real_
    out <- sweep(counts, 1, rs, "/")
    out[is.na(out)] <- 0
  } else if (mode == "log1p_counts") {
    out <- log1p(counts)
  } else if (mode == "log1p_relative") {
    rs <- rowSums(counts, na.rm = TRUE)
    rs[rs == 0] <- NA_real_
    rel <- sweep(counts, 1, rs, "/")
    rel[is.na(rel)] <- 0
    out <- log1p(rel)
  } else {
    stop(
      "Unsupported --abundance_mode: ", mode,
      ". Valid: hellinger, counts, relative, log1p_counts, log1p_relative",
      call. = FALSE
    )
  }

  out[!is.finite(out)] <- 0
  out
}

filter_taxa_on_raw_counts <- function(counts, min_prev_n, min_total) {
  prev <- colSums(counts > 0, na.rm = TRUE)
  total <- colSums(counts, na.rm = TRUE)
  keep <- prev >= min_prev_n & total >= min_total
  counts[, keep, drop = FALSE]
}

safe_titan <- function(env, txa, min_splt, n_perm, n_boot, ncpus, seed) {
  set.seed(seed)

  ord <- order(env)
  env <- env[ord]
  txa <- txa[ord, , drop = FALSE]

  titan_fun <- get("titan", asNamespace("TITAN2"))
  fml <- names(formals(titan_fun))

  named_args <- list(
    minSplt = min_splt,
    numPerm = n_perm,
    boot = TRUE,
    nBoot = n_boot,
    imax = FALSE,
    ivTot = FALSE,
    pur.cut = purity_cut,
    rel.cut = reliability_cut,
    ncpus = ncpus,
    memory = FALSE,
    messaging = FALSE
  )

  named_args <- named_args[names(named_args) %in% fml]
  call_args <- c(list(env = env, txa = txa), named_args)

  do.call(titan_fun, call_args)
}

find_col <- function(df, candidates) {
  nms <- names(df)
  nms_norm <- tolower(gsub("[^a-z0-9]+", "", nms))
  cand_norm <- tolower(gsub("[^a-z0-9]+", "", candidates))

  idx <- match(cand_norm, nms_norm)
  idx <- idx[!is.na(idx)]

  if (length(idx) == 0) return(NA_character_)
  nms[idx[1]]
}

to_num <- function(x) suppressWarnings(as.numeric(x))

direction_from_maxgrp <- function(x) {
  x_chr <- tolower(as.character(x))

  case_when(
    x_chr %in% c("1", "z-", "-", "neg", "negative", "decrease", "decreasing") ~ "z-",
    x_chr %in% c("2", "z+", "+", "pos", "positive", "increase", "increasing") ~ "z+",
    TRUE ~ as.character(x)
  )
}

extract_sumz_from_titan <- function(titan_obj) {
  x <- as.data.frame(titan_obj$sumz.cp, check.names = FALSE)
  components <- c("sumz-", "sumz+", "fsumz-", "fsumz+")
  cols <- c("cp", "0.05", "0.10", "0.50", "0.90", "0.95")
  stopifnot(all(components %in% rownames(x)),
            all(cols %in% names(x)))
  x <- x[components, cols, drop = FALSE]
  tibble(
    component = components,
    response = ifelse(grepl("+", components, fixed = TRUE), "z+", "z-"),
    filtered = startsWith(components, "f"),
    cp_plot = as.numeric(x[["cp"]]),
    ci05 = as.numeric(x[["0.05"]]),
    ci10 = as.numeric(x[["0.10"]]),
    ci50 = as.numeric(x[["0.50"]]),
    ci90 = as.numeric(x[["0.90"]]),
    ci95 = as.numeric(x[["0.95"]])
  )
}

extract_indval_from_titan <- function(titan_obj) {
  x <- as.data.frame(titan_obj$sppmax, check.names = FALSE)
  required <- c("zenv.cp", "zscore", "maxgrp", "purity",
                "reliability", "freq", "5%", "10%", "50%",
                "90%", "95%", "filter")
  if (!all(required %in% names(x)) &&
      all(required %in% rownames(x))) {
    x <- as.data.frame(t(as.matrix(x)), check.names = FALSE)
  }
  stopifnot(all(required %in% names(x)),
            !anyDuplicated(rownames(x)))
  num <- function(n) as.numeric(x[[n]])
  out <- tibble(
    taxon_id = rownames(x),
    obs_cp = num("zenv.cp"),
    z_score = num("zscore"),
    maxgrp_raw = as.character(x$maxgrp),
    response = direction_from_maxgrp(maxgrp_raw),
    purity = num("purity"),
    reliability = num("reliability"),
    frequency = num("freq"),
    cp_05 = num("5%"),
    cp_10 = num("10%"),
    cp_50 = num("50%"),
    cp_90 = num("90%"),
    cp_95 = num("95%"),
    saved_filter = num("filter")
  ) %>%
    mutate(
      cp_plot = obs_cp,
      abs_z = abs(z_score),
      reliable90 = is.finite(purity) & is.finite(reliability) &
        purity >= 0.90 & reliability >= 0.90,
      reliable95 = is.finite(purity) & is.finite(reliability) &
        purity >= 0.95 & reliability >= 0.95
    )
  stopifnot(all(is.finite(out$obs_cp)),
            all(is.finite(out$cp_05)),
            all(is.finite(out$cp_95)),
            all(out$response %in% c("z-", "z+")))
  out
}

safe_numeric_gradient <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA_real_
  x
}

nice_gradient_label <- function(x) {
  dplyr::recode(
    x,
    "ECI_PC1" = "Integrated environmental condition",
    "MHI_local" = "Mangrove Health Index (MHI)",
    "water_inundation_PC1" = "Water / inundation",
    "vegetation_landscape_PC1" = "Vegetation / landscape",
    "physicochemical_PC1" = "Physicochemical",
    "nutrients_redox_PC1" = "Nutrients / redox",
    "moisture_stress_PC1" = "Moisture stress",
    .default = x
  )
}

save_session_info <- function(path) {
  capture.output(sessionInfo(), file = path)
}

# ------------------------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------------------------

if (!file.exists(ps_file)) stop("Phyloseq file not found: ", ps_file, call. = FALSE)
if (!dir.exists(out202_dir)) stop("OUT202 directory not found: ", out202_dir, call. = FALSE)
if (!file.exists(meta_csv)) stop("Metadata CSV not found: ", meta_csv, call. = FALSE)

cat("\n--- Reading inputs ---\n")
cat("Phyloseq: ", ps_file, "\n", sep = "")
cat("OUT202: ", out202_dir, "\n", sep = "")
cat("Metadata: ", meta_csv, "\n", sep = "")

ps <- readRDS(ps_file)
meta <- readr::read_csv(meta_csv, show_col_types = FALSE, progress = FALSE)

if (!sample_id_col %in% names(meta)) stop("sample_id column not found: ", sample_id_col, call. = FALSE)
if (!depth_col %in% names(meta)) stop("depth column not found: ", depth_col, call. = FALSE)
if (!locality_col %in% names(meta)) stop("locality column not found: ", locality_col, call. = FALSE)

profile_col <- intersect(
  c("profile_id", "lat_block", "latitude_block", "perfil", "core_id"),
  names(meta)
)[1]
if (is.na(profile_col)) {
  stop("No profile_id/lat_block column found; depth-stratified TITAN2 requires one independent observation per profile.", call. = FALSE)
}

# Robust ECI alias if needed.
if ("ECI_PC1" %in% gradients && !"ECI_PC1" %in% names(meta)) {
  eci_candidates <- intersect(
    c("ECI.PC1", "eci_pc1"),
    names(meta)
  )

  if (length(eci_candidates) > 0) {
    meta$ECI_PC1 <- meta[[eci_candidates[1]]]
    cat("Using ECI alias: ", eci_candidates[1], " -> ECI_PC1\n", sep = "")
  }
}

missing_gradients <- setdiff(gradients, names(meta))
if (length(missing_gradients) > 0) {
  stop("Missing gradient columns in metadata: ", paste(missing_gradients, collapse = ", "), call. = FALSE)
}

bimodal_taxa_raw <- read_bimodal_taxa_from_out202(out202_dir)

cat("Bimodal taxa read from OUT202: ", length(bimodal_taxa_raw), "\n", sep = "")

counts_all <- extract_phyloseq_counts_samples_by_taxa(ps)
tax_map <- make_taxonomy_map(ps)

write_csv_safe(
  tax_map,
  file.path(tables_dir, "405_taxonomy_map_from_phyloseq.csv")
)

taxa_in_ps <- colnames(counts_all)
bimodal_in_ps <- match_bimodal_taxa_to_phyloseq(bimodal_taxa_raw, taxa_in_ps)

if (length(bimodal_in_ps) == 0) {
  stop(
    "No OUT202 bimodal taxa matched taxa_names(ps). Check taxon IDs between OUT202 and phyloseq.",
    call. = FALSE
  )
}

viral_taxa <- tax_map %>%
  filter(.data$is_viral) %>%
  pull(.data$taxon_id) %>%
  unique()

selected_bimodal_taxa <- bimodal_in_ps

cat("Bimodal taxa matched in phyloseq: ", length(bimodal_in_ps), "\n", sep = "")
cat("Viral taxa excluded: 0; taxonomy is annotation only.\n")
cat("Raw1031 taxa retained before prevalence filtering: ", length(selected_bimodal_taxa), "\n", sep = "")

counts <- counts_all[, selected_bimodal_taxa, drop = FALSE]

sample_ids_ps <- rownames(counts)

meta2 <- meta %>%
  mutate(
    .sample_id = as.character(.data[[sample_id_col]]),
    .depth = normalize_depth(.data[[depth_col]]),
    .locality = as.character(.data[[locality_col]]),
    .profile_raw = as.character(.data[[profile_col]]),
    .profile_id = if_else(
      str_detect(.profile_raw, fixed(.locality)),
      .profile_raw,
      paste(.locality, .profile_raw, sep = "::")
    )
  ) %>%
  filter(.sample_id %in% sample_ids_ps)

common_samples <- intersect(sample_ids_ps, meta2$.sample_id)

if (length(common_samples) != 51L || !setequal(sample_ids_ps, meta2$.sample_id)) {
  stop("405 requiere coincidencia exacta de 51 sample_id entre phyloseq y metadata.", call. = FALSE)
}

counts <- counts[common_samples, , drop = FALSE]

meta2 <- meta2 %>%
  filter(.sample_id %in% common_samples) %>%
  arrange(match(.sample_id, common_samples))

stopifnot(identical(rownames(counts), meta2$.sample_id),
          ncol(counts) == 1031L,
          all(is.finite(counts)), all(counts >= 0),
          all(rowSums(counts) > 0),
          !anyNA(meta2$.profile_raw), all(nzchar(meta2$.profile_raw)),
          !anyNA(meta2$.locality), all(nzchar(meta2$.locality)),
          setequal(depth_levels, c("5", "20", "40")))
for (g in gradients) {
  if (any(!is.finite(safe_numeric_gradient(meta2[[g]])))) {
    stop("Gradiente incompleto: ", g,
         "; no se excluirán perfiles silenciosamente.")
  }
}
write_csv_safe(
  tibble(universe = "raw1031", selection = "observed_joint_raw",
         n_taxa = ncol(counts), n_samples = nrow(counts),
         n_profiles = n_distinct(meta2$.profile_id),
         n_viral_excluded = 0L),
  file.path(tables_dir, "405_contrato_universo_raw1031.csv")
)

cat("Matched samples: ", nrow(meta2), "\n", sep = "")
cat("Depth levels requested: ", paste(depth_levels, collapse = ", "), "\n", sep = "")
cat("Depth levels observed: ", paste(sort(unique(meta2$.depth)), collapse = ", "), "\n", sep = "")
cat("Gradients requested: ", paste(gradients, collapse = ", "), "\n", sep = "")

write_csv_safe(
  meta2,
  file.path(tables_dir, "405_metadata_used_all_samples.csv")
)

sample_counts <- meta2 %>%
  group_by(.depth) %>%
  summarise(n_samples = n(), n_profiles = n_distinct(.profile_id), .groups = "drop") %>%
  arrange(.depth)

profile_gate <- meta2 %>% group_by(.profile_id) %>%
  summarise(n = n(), depths = paste(sort(unique(.depth)), collapse = "|"),
            valid = n == 3L && setequal(.depth, c("5", "20", "40")) &&
              n_distinct(.locality) == 1L,
            .groups = "drop")
write_csv_safe(profile_gate, file.path(tables_dir, "405_auditoria_canon_51.csv"))
if (nrow(meta2) != 51L || n_distinct(meta2$.sample_id) != 51L ||
    nrow(profile_gate) != 17L || !all(profile_gate$valid) ||
    any(sample_counts$n_samples != sample_counts$n_profiles)) {
  stop("405 requiere exactamente 51 muestras en 17 perfiles completos, sin perfiles duplicados por profundidad.", call. = FALSE)
}

write_csv_safe(
  sample_counts,
  file.path(tables_dir, "405_sample_counts_by_depth.csv")
)

# ------------------------------------------------------------------------------
# Run information
# ------------------------------------------------------------------------------

run_info <- tibble(
  key = c(
    "script_id",
    "script_name",
    "timestamp",
    "ps_file",
    "out202_dir",
    "meta_csv",
    "results_root",
    "out_dir",
    "sample_id_col",
    "depth_col",
    "condition_col",
    "locality_col",
    "gradients",
    "depth_levels",
    "abundance_mode",
    "min_prev_n",
    "min_total",
    "min_splt",
    "n_boot",
    "n_perm",
    "ncpus",
    "seed",
    "auto_relax",
    "relaxed_min_prev_n",
    "relaxed_min_splt",
    "min_samples_per_depth",
    "min_gradient_unique",
    "reliability_cut",
    "purity_cut",
    "n_samples_matched",
    "n_bimodal_taxa_out202",
    "n_bimodal_taxa_matched_ps",
    "n_viral_taxa_excluded",
    "n_selected_bimodal_taxa_initial"
  ),
  value = c(
    SCRIPT_ID,
    SCRIPT_NAME,
    timestamp,
    ps_file,
    out202_dir,
    meta_csv,
    results_root,
    out_dir,
    sample_id_col,
    depth_col,
    condition_col,
    locality_col,
    paste(gradients, collapse = ","),
    paste(depth_levels, collapse = ","),
    abundance_mode,
    as.character(min_prev_n),
    as.character(min_total),
    as.character(min_splt),
    as.character(n_boot),
    as.character(n_perm),
    as.character(ncpus),
    as.character(seed),
    as.character(auto_relax),
    as.character(relaxed_min_prev_n),
    as.character(relaxed_min_splt),
    as.character(min_samples_per_depth),
    as.character(min_gradient_unique),
    as.character(reliability_cut),
    as.character(purity_cut),
    as.character(nrow(meta2)),
    as.character(length(bimodal_taxa_raw)),
    as.character(length(bimodal_in_ps)),
    as.character(0L),
    as.character(length(selected_bimodal_taxa))
  )
)

write_csv_safe(run_info, file.path(tables_dir, "405_run_info.csv"))

# ------------------------------------------------------------------------------
# TITAN2 loop
# ------------------------------------------------------------------------------

summary_rows <- list()
sumz_rows <- list()
indval_rows <- list()

analysis_index <- 1L

for (depth_i in depth_levels) {
  cat("\n==============================\n")
  cat("Depth: ", depth_i, "\n", sep = "")
  cat("==============================\n")

  depth_keep <- meta2$.depth == as.character(depth_i)

  meta_depth <- meta2[depth_keep, , drop = FALSE]
  counts_depth <- counts[depth_keep, , drop = FALSE]

  if (anyDuplicated(meta_depth$.profile_id)) {
    stop("Perfil duplicado en profundidad ", depth_i, call. = FALSE)
  }

  n_depth_samples <- nrow(meta_depth)

  for (grad in gradients) {
    cat("\n--- Running depth ", depth_i, " | gradient ", grad, " ---\n", sep = "")

    status <- "not_run"
    message <- NA_character_

    n_samples <- NA_integer_
    n_gradient_unique <- NA_integer_
    n_taxa_retained <- NA_integer_
    min_prev_use <- min_prev_n
    min_splt_use <- min_splt

    grad_vec <- safe_numeric_gradient(meta_depth[[grad]])

    complete <- !is.na(grad_vec)
    meta_g <- meta_depth[complete, , drop = FALSE]
    counts_g <- counts_depth[complete, , drop = FALSE]
    grad_vec <- grad_vec[complete]

    n_samples <- nrow(meta_g)
    n_gradient_unique <- length(unique(grad_vec))

    if (auto_relax && n_samples < 18) {
      min_prev_use <- max(3L, min(min_prev_use, relaxed_min_prev_n))
      min_splt_use <- max(5L, min(min_splt_use, relaxed_min_splt))
    }

    if (n_samples < min_samples_per_depth) {
      status <- "skipped"
      message <- paste0(
        "Too few samples after depth/gradient filtering: n=", n_samples,
        "; required >= ", min_samples_per_depth
      )
      cat("SKIP: ", message, "\n", sep = "")
    } else if (n_gradient_unique < min_gradient_unique) {
      status <- "skipped"
      message <- paste0(
        "Too few unique gradient values: unique=", n_gradient_unique,
        "; required >= ", min_gradient_unique
      )
      cat("SKIP: ", message, "\n", sep = "")
    } else if (min_splt_use * 2 >= n_samples) {
      status <- "skipped"
      message <- paste0(
        "minSplt too large for sample size: minSplt=", min_splt_use,
        ", n=", n_samples
      )
      cat("SKIP: ", message, "\n", sep = "")
    } else {
      counts_filtered <- filter_taxa_on_raw_counts(
        counts_g,
        min_prev_n = min_prev_use,
        min_total = min_total
      )

      n_taxa_retained <- ncol(counts_filtered)

      if (n_taxa_retained < 2) {
        status <- "skipped"
        message <- paste0("Too few taxa retained after raw-count filtering: ", n_taxa_retained)
        cat("SKIP: ", message, "\n", sep = "")
      } else {
        txa_filtered <- prepare_abundance_matrix(
          counts_g,
          mode = abundance_mode
        )[, colnames(counts_filtered), drop = FALSE]

        var_taxa <- apply(txa_filtered, 2, var, na.rm = TRUE)
        keep_var <- is.finite(var_taxa) & var_taxa > 0
        txa_filtered <- txa_filtered[, keep_var, drop = FALSE]
n_taxa_retained <- ncol(txa_filtered)
        write_csv_safe(
          tibble(
            taxon_id = colnames(counts_g),
            prevalence = colSums(counts_g > 0),
            total = colSums(counts_g),
            pass_raw_filter = colnames(counts_g) %in% colnames(counts_filtered),
            retained = colnames(counts_g) %in% colnames(txa_filtered)
          ),
          file.path(tables_dir,
                    paste0("405_taxa_audit_depth", depth_i, "__", grad, ".csv"))
        )

        if (n_taxa_retained < 2) {
          status <- "skipped"
          message <- paste0("Too few variable taxa retained after transform: ", n_taxa_retained)
          cat("SKIP: ", message, "\n", sep = "")
        } else {
          status <- "ok"
          message <- "TITAN2 completed"

          cat("Samples: ", n_samples, "\n", sep = "")
          cat("Unique gradient values: ", n_gradient_unique, "\n", sep = "")
          cat("Taxa retained: ", n_taxa_retained, "\n", sep = "")
          cat("min_prev_n used: ", min_prev_use, "\n", sep = "")
          cat("minSplt used: ", min_splt_use, "\n", sep = "")

          titan_id <- paste0(
            "depth", depth_i,
            "__",
            grad
          )

          rds_path <- file.path(rds_dir, paste0("405_titan_", titan_id, ".rds"))

          titan_obj <- tryCatch(
            {
              safe_titan(
                env = grad_vec,
                txa = txa_filtered,
                min_splt = min_splt_use,
                n_perm = n_perm,
                n_boot = n_boot,
                ncpus = ncpus,
                seed = seed + analysis_index
              )
            },
            error = function(e) {
              status <<- "failed"
              message <<- conditionMessage(e)
              NULL
            }
          )

          if (!is.null(titan_obj)) {
            saveRDS(titan_obj, rds_path)

            sumz_tbl <- tryCatch(
              extract_sumz_from_titan(titan_obj),
              error = function(e) {
                warning("Could not extract sumz for ", titan_id, ": ", conditionMessage(e))
                tibble()
              }
            )

            if (nrow(sumz_tbl) > 0) {
              sumz_tbl <- sumz_tbl %>%
                mutate(
                  depth_cm = depth_i,
                  gradient = grad,
                  gradient_label = nice_gradient_label(grad),
                  n_samples = n_samples,
                  n_taxa_retained = n_taxa_retained,
                  min_prev_n_used = min_prev_use,
                  min_splt_used = min_splt_use,
                  titan_rds = rds_path,
                  .before = 1
                )

              sumz_rows[[length(sumz_rows) + 1L]] <- sumz_tbl
            }

            indval_tbl <- tryCatch(
              extract_indval_from_titan(titan_obj),
              error = function(e) {
                warning("Could not extract indval for ", titan_id, ": ", conditionMessage(e))
                tibble()
              }
            )

extraction_ok <- nrow(sumz_tbl) == 4L &&
              nrow(indval_tbl) == ncol(txa_filtered) &&
              setequal(indval_tbl$taxon_id, colnames(txa_filtered))
            if (!extraction_ok) {
              status <- "failed_extraction"
              message <- "RDS guardado; extracción incompleta o IDs incompatibles"
            }
            if (nrow(indval_tbl) > 0) {
              indval_tbl <- indval_tbl %>%
                mutate(
                  depth_cm = depth_i,
                  gradient = grad,
                  gradient_label = nice_gradient_label(grad),
                  n_samples = n_samples,
                  n_taxa_retained = n_taxa_retained,
                  min_prev_n_used = min_prev_use,
                  min_splt_used = min_splt_use,
                  titan_rds = rds_path,
                  .before = 1
                )

              indval_rows[[length(indval_rows) + 1L]] <- indval_tbl
            }
          }
        }
      }
    }

    summary_rows[[length(summary_rows) + 1L]] <- tibble(
      depth_cm = depth_i,
      gradient = grad,
      gradient_label = nice_gradient_label(grad),
      status = status,
      message = message,
      n_samples_depth = n_depth_samples,
      n_samples_used = n_samples,
      n_gradient_unique = n_gradient_unique,
      n_taxa_retained = n_taxa_retained,
      min_prev_n_used = min_prev_use,
      min_splt_used = min_splt_use
    )

    analysis_index <- analysis_index + 1L
  }
}

depth_summary <- bind_rows(summary_rows)
sumz_all <- bind_rows(sumz_rows)
indval_all <- bind_rows(indval_rows)

# ------------------------------------------------------------------------------
# Annotate indicator taxa
# ------------------------------------------------------------------------------

tax_annot <- tax_map %>%
  mutate(taxon_id = as.character(taxon_id)) %>%
  select(
    taxon_id,
    tax_label,
    full_taxonomy,
    is_viral,
    everything()
  )

if (nrow(indval_all) > 0) {
  indval_annotated <- indval_all %>%
    mutate(taxon_id = as.character(taxon_id)) %>%
    left_join(tax_annot, by = "taxon_id") %>%
    mutate(
      tax_label = ifelse(is.na(tax_label) | tax_label == "", taxon_id, tax_label),
      full_taxonomy = ifelse(is.na(full_taxonomy) | full_taxonomy == "", tax_label, full_taxonomy),
      abs_z = abs(z_score),
      reliable90 = is.finite(purity) & is.finite(reliability) &
        purity >= 0.90 & reliability >= 0.90,
      reliable95 = is.finite(purity) & is.finite(reliability) &
        purity >= 0.95 & reliability >= 0.95,
      reliable_custom = is.finite(purity) & is.finite(reliability) &
        purity >= purity_cut & reliability >= reliability_cut
    )
} else {
  indval_annotated <- tibble()
}

# ------------------------------------------------------------------------------
# Write tables
# ------------------------------------------------------------------------------

write_csv_safe(
  depth_summary,
  file.path(tables_dir, "405_titan_depth_summary.csv")
)

write_csv_safe(
  sumz_all,
  file.path(tables_dir, "405_titan_sumz_by_depth_gradient.csv")
)

write_csv_safe(
  indval_annotated,
  file.path(tables_dir, "405_titan_indval_all_depths_annotated.csv")
)

if (nrow(indval_annotated) > 0) {
  reliable_summary <- indval_annotated %>%
    group_by(depth_cm, gradient, gradient_label, response) %>%
    summarise(
      n_indicators = n(),
      n_reliable90 = sum(reliable90, na.rm = TRUE),
      n_reliable95 = sum(reliable95, na.rm = TRUE),
      n_reliable_custom = sum(reliable_custom, na.rm = TRUE),
      median_cp = median(cp_plot, na.rm = TRUE),
      median_abs_z = median(abs_z, na.rm = TRUE),
      max_abs_z = max(abs_z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      median_cp = ifelse(is.nan(median_cp), NA_real_, median_cp),
      median_abs_z = ifelse(is.nan(median_abs_z), NA_real_, median_abs_z),
      max_abs_z = ifelse(is.infinite(max_abs_z), NA_real_, max_abs_z)
    )

  write_csv_safe(
    reliable_summary,
    file.path(tables_dir, "405_titan_indval_reliable_summary.csv")
  )

  top_indicators <- indval_annotated %>%
    filter(reliable_custom) %>%
    arrange(depth_cm, gradient, desc(abs_z)) %>%
    group_by(depth_cm, gradient) %>%
    slice_head(n = 30) %>%
    ungroup()

  write_csv_safe(
    top_indicators,
    file.path(tables_dir, "405_titan_top_reliable_indicators_by_depth_gradient.csv")
  )
}

# Robust filtered-sumz table.
# Some TITAN2 versions or extraction paths may not return an explicit
# `filtered` column. In that case, keep all z-/z+ community threshold rows.
if (nrow(sumz_all) == 0) {
  filtered_sumz <- tibble()
} else {
  if (!"filtered" %in% names(sumz_all)) {
    sumz_all$filtered <- TRUE
  }

  if (!"response" %in% names(sumz_all)) {
    sumz_all$response <- NA_character_
  }

  filtered_sumz <- sumz_all %>%
    filter(.data$filtered, .data$response %in% c("z-", "z+")) %>%
    arrange(.data$depth_cm, .data$gradient, .data$response)
}

write_csv_safe(
  filtered_sumz,
  file.path(tables_dir, "405_titan_filtered_sumz_by_depth_gradient.csv")
)

# ------------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------------

theme_405 <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_text(color = "black", margin = margin(b = 6)),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.35),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.25),
      plot.caption = element_text(color = "grey30", hjust = 0)
    )
}

pal_response <- c(
  "z-" = "#D55E00",
  "z+" = "#009E73"
)

if (nrow(filtered_sumz) > 0) {
  sumz_plot <- filtered_sumz %>%
    filter(is.finite(cp_plot)) %>%
    mutate(
      depth_cm = factor(depth_cm, levels = depth_levels),
      gradient_label = factor(
        gradient_label,
        levels = nice_gradient_label(gradients)
      ),
      response = factor(response, levels = c("z-", "z+"))
    )

  if (nrow(sumz_plot) > 0) {
    p_forest <- ggplot(
      sumz_plot,
      aes(
        x = cp_plot,
        y = depth_cm,
        color = response,
        fill = response,
        shape = response
      )
    ) +
      geom_vline(xintercept = 0, linewidth = 0.30, linetype = "dashed", color = "grey45") +
      geom_errorbarh(
        aes(xmin = ci05, xmax = ci95),
        height = 0.18,
        linewidth = 0.60,
        alpha = 0.30,
        na.rm = TRUE
      ) +
      geom_errorbarh(
        aes(xmin = ci10, xmax = ci90),
        height = 0.09,
        linewidth = 1.05,
        alpha = 0.90,
        na.rm = TRUE
      ) +
      geom_point(size = 3.0, stroke = 0.75, color = "black", na.rm = TRUE) +
      facet_wrap(~ gradient_label, scales = "free_x", ncol = 2) +
      scale_color_manual(values = pal_response) +
      scale_fill_manual(values = pal_response) +
      scale_shape_manual(values = c("z-" = 21, "z+" = 24)) +
      scale_x_continuous(breaks = pretty_breaks(n = 5)) +
      labs(
        title = NULL,
        subtitle = "Depth-specific TITAN2 environmental thresholds",
        x = "Environmental gradient value",
        y = "Depth (cm)",
        color = "Response",
        fill = "Response",
        shape = "Response",
        caption = "Points are filtered community change points; thick and thin intervals denote 10–90% and 5–95% bootstrap ranges."
      ) +
      theme_405(base_size = 11)

    ggsave(
      filename = file.path(plots_dir, "405_threshold_forest_by_depth.png"),
      plot = p_forest,
      width = 10.8,
      height = 7.2,
      dpi = 600,
      bg = "white"
    )

    ggsave(
      filename = file.path(plots_dir, "405_threshold_forest_by_depth.pdf"),
      plot = p_forest,
      width = 10.8,
      height = 7.2,
      device = cairo_pdf,
      bg = "white"
    )

    heat_df <- sumz_plot %>%
      mutate(
        label = ifelse(is.na(cp_plot), "", sprintf("%.2f", cp_plot)),
        tile_id = factor(
          paste0(response, " | ", depth_cm, " cm"),
          levels = rev(paste0(as.vector(outer(c("z-", "z+"), depth_levels, paste, sep = " | ")), " cm"))
        )
      )

    p_heat <- ggplot(
      heat_df,
      aes(
        x = gradient_label,
        y = tile_id,
        fill = cp_plot
      )
    ) +
      geom_tile(color = "white", linewidth = 0.35) +
      geom_text(aes(label = label), size = 3) +
      scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        name = "Change point"
      ) +
      labs(
        title = NULL,
        subtitle = "Depth-specific threshold map",
        x = NULL,
        y = "Response | depth",
        caption = "Cell values are filtered TITAN2 community change points."
      ) +
      theme_405(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1)
      )

    ggsave(
      filename = file.path(plots_dir, "405_threshold_heatmap_depth_gradient.png"),
      plot = p_heat,
      width = 9.8,
      height = 5.8,
      dpi = 600,
      bg = "white"
    )

    ggsave(
      filename = file.path(plots_dir, "405_threshold_heatmap_depth_gradient.pdf"),
      plot = p_heat,
      width = 9.8,
      height = 5.8,
      device = cairo_pdf,
      bg = "white"
    )
  }
}

if (nrow(indval_annotated) > 0) {
  density_df <- indval_annotated %>%
    filter(
      is.finite(cp_plot),
      response %in% c("z-", "z+")
    ) %>%
    mutate(
      depth_cm = factor(depth_cm, levels = depth_levels),
      gradient_label = factor(
        gradient_label,
        levels = nice_gradient_label(gradients)
      ),
      response = factor(response, levels = c("z-", "z+"))
    ) %>%
    filter(reliable_custom)

  if (nrow(density_df) >= 10) {
    p_density <- ggplot(
      density_df,
      aes(
        x = cp_plot,
        color = response,
        fill = response
      )
    ) +
      geom_density(alpha = 0.22, linewidth = 0.75, na.rm = TRUE) +
      geom_rug(alpha = 0.12, linewidth = 0.2) +
      facet_grid(depth_cm ~ gradient_label, scales = "free_x") +
      scale_color_manual(values = pal_response) +
      scale_fill_manual(values = pal_response) +
      labs(
        title = NULL,
        subtitle = "Taxon-specific TITAN2 change-point densities by depth",
        x = "Taxon-specific change point",
        y = "Density",
        color = "Response",
        fill = "Response",
        caption = paste0("Reliable indicators shown when available; purity and reliability cut = ", reliability_cut, ".")
      ) +
      theme_405(base_size = 9.5) +
      theme(
        axis.text.x = element_text(size = 7.5),
        strip.text = element_text(size = 8.5)
      )

    ggsave(
      filename = file.path(plots_dir, "405_indicator_cp_density_by_depth.png"),
      plot = p_density,
      width = 13,
      height = 8.5,
      dpi = 600,
      bg = "white"
    )

    ggsave(
      filename = file.path(plots_dir, "405_indicator_cp_density_by_depth.pdf"),
      plot = p_density,
      width = 13,
      height = 8.5,
      device = cairo_pdf,
      bg = "white"
    )
  }
}

# ------------------------------------------------------------------------------
# Figure index
# ------------------------------------------------------------------------------

fig_index <- tibble(
  figure = c(
    "405_threshold_forest_by_depth.png",
    "405_threshold_heatmap_depth_gradient.png",
    "405_indicator_cp_density_by_depth.png"
  ),
  role = c(
    "Main candidate: forest plot of z- and z+ community thresholds by depth and gradient",
    "Compact summary: threshold map by depth, gradient and response",
    "Supplementary/Extended Data: taxon-specific change point density by depth"
  ),
  path = file.path(plots_dir, figure),
  exists = file.exists(path)
)

write_csv_safe(
  fig_index,
  file.path(tables_dir, "405_figures_index.csv")
)

# ------------------------------------------------------------------------------
# Session info and final report
# ------------------------------------------------------------------------------

save_session_info(file.path(logs_dir, "405_sessionInfo.txt"))
if (nrow(depth_summary) != length(depth_levels) * length(gradients) ||
    any(depth_summary$status != "ok")) {
  print(depth_summary)
  stop("405 incompleto: revisar resumen y RDS; LATEST no actualizado.")
}
writeLines(out_dir, latest_file)

cat("\n===== Finished 405 TITAN2 depth-stratified analysis =====\n")
cat("End time: ", as.character(Sys.time()), "\n", sep = "")
cat("Output directory: ", out_dir, "\n", sep = "")
cat("LATEST file: ", latest_file, "\n", sep = "")

cat("\nStatus summary:\n")
print(depth_summary)

cat("\nFiles written:\n")
cat("  - ", file.path(tables_dir, "405_titan_depth_summary.csv"), "\n", sep = "")
cat("  - ", file.path(tables_dir, "405_titan_sumz_by_depth_gradient.csv"), "\n", sep = "")
cat("  - ", file.path(tables_dir, "405_titan_filtered_sumz_by_depth_gradient.csv"), "\n", sep = "")
cat("  - ", file.path(tables_dir, "405_titan_indval_all_depths_annotated.csv"), "\n", sep = "")
cat("  - ", file.path(plots_dir, "405_threshold_forest_by_depth.png"), "\n", sep = "")
cat("  - ", file.path(plots_dir, "405_threshold_heatmap_depth_gradient.png"), "\n", sep = "")
cat("  - ", file.path(plots_dir, "405_indicator_cp_density_by_depth.png"), "\n", sep = "")
