#!/usr/bin/env Rscript
# ============================================================
# 04_03_annotate_indicators.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Joins taxonomic annotation to TITAN indicators.
# Inputs (source expressions; complete list in docs/contracts/04_03_annotate_indicators.json):
#   RUN001 <- if (file.exists(LATEST001)) trimws(readLines(LATEST001, warn = FALSE, n = 1L)) else NA_character_
#   DEFAULT401 <- if (file.exists(LATEST401)) trimws(readLines(LATEST401, warn = FALSE, n = 1L)) else NA_character_
#   ps <- readRDS(PS_RDS)
#   taxa_used <- readr::read_csv(taxa_used_file, show_col_types = FALSE) %>%
#   ind <- readr::read_csv(f, show_col_types = FALSE)
# Outputs (source expressions; complete list in contract):
#   readr::write_csv(
#   readr::write_csv(ann, out_file)
# Algorithmic provenance:
# Join phyloseq taxonomy to the saved taxon-level TITAN2 indicator tables.
#   Baker & King (2010), doi:10.1111/j.2041-210X.2009.00007.x.
# Source SHA-256: 5be068821821c00be25ecf54c83fbaec338cc468c802ca8c551c32cf7341c300
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


# ==============================================================================
# 403_annotate_TITAN2_indval_taxonomy.R
#
# Objetivo:
#   Anotar taxonomía completa en las tablas TITAN2:
#     401_titan_indval_<gradient>.csv
#
#   El campo "taxon" suele ser un ID numérico/tax_id. Este script cruza esos IDs
#   contra la tax_table del phyloseq y contra 031_nonviral_bimodal_taxa_used.csv
#   si existe.
#
# Outputs:
#   tables/
#     403_titan_indval_<gradient>_annotated.csv
#     403_titan_indval_all_gradients_annotated.csv
#     403_titan_indval_reliable95_annotated.csv
#     403_titan_indval_reliable90_annotated.csv
#     403_titan_taxonomy_match_audit.csv
#     403_titan_taxonomic_response_summary.csv
#     403_taxonomy_map_from_phyloseq.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

# ==============================================================================
# Helpers
# ==============================================================================

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[[hit + 1]]
}

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
  cat("\n")
  flush.console()
}

safe_mkdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

stop_if_missing_file <- function(path, label = "file") {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
    stop(sprintf("Missing %s: %s", label, path), call. = FALSE)
  }
}

clean_tax_string <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""

  x <- gsub("^D_[0-9]+__", "", x)
  x <- gsub("^[a-zA-Z]__", "", x)
  x <- gsub("^.__", "", x)
  x <- gsub("_", " ", x)
  x <- trimws(x)

  x[x %in% c(
    "",
    "NA",
    "NaN",
    "unclassified",
    "Unclassified",
    "unknown",
    "Unknown",
    "not assigned",
    "Not assigned"
  )] <- NA_character_

  x
}

make_tax_label <- function(tax_df, tax_id_col = "tax_id") {
  ranks <- names(tax_df)

  species_col <- intersect(c("Species", "species", "S", "Rank_S"), ranks)[1]
  genus_col   <- intersect(c("Genus", "genus", "G", "Rank_G"), ranks)[1]
  family_col  <- intersect(c("Family", "family", "F", "Rank_F"), ranks)[1]
  order_col   <- intersect(c("Order", "order", "O", "Rank_O"), ranks)[1]
  class_col   <- intersect(c("Class", "class", "C", "Rank_C"), ranks)[1]
  phylum_col  <- intersect(c("Phylum", "phylum", "P", "Rank_P"), ranks)[1]

  species <- if (!is.na(species_col)) clean_tax_string(tax_df[[species_col]]) else rep(NA_character_, nrow(tax_df))
  genus   <- if (!is.na(genus_col)) clean_tax_string(tax_df[[genus_col]]) else rep(NA_character_, nrow(tax_df))
  family  <- if (!is.na(family_col)) clean_tax_string(tax_df[[family_col]]) else rep(NA_character_, nrow(tax_df))
  order   <- if (!is.na(order_col)) clean_tax_string(tax_df[[order_col]]) else rep(NA_character_, nrow(tax_df))
  class_  <- if (!is.na(class_col)) clean_tax_string(tax_df[[class_col]]) else rep(NA_character_, nrow(tax_df))
  phylum  <- if (!is.na(phylum_col)) clean_tax_string(tax_df[[phylum_col]]) else rep(NA_character_, nrow(tax_df))

  label <- species

  bad_species <- is.na(label) |
    label == "" |
    grepl("^sp\\.?$", label, ignore.case = TRUE) |
    grepl("uncultured|unclassified|metagenome|environmental sample", label, ignore.case = TRUE)

  idx <- bad_species & !is.na(genus)
  label[idx] <- paste0(genus[idx], " sp.")

  bad_label <- is.na(label) | label == ""
  idx <- bad_label & !is.na(family)
  label[idx] <- paste0(family[idx], " taxon")

  bad_label <- is.na(label) | label == ""
  idx <- bad_label & !is.na(order)
  label[idx] <- paste0(order[idx], " taxon")

  bad_label <- is.na(label) | label == ""
  idx <- bad_label & !is.na(class_)
  label[idx] <- paste0(class_[idx], " taxon")

  bad_label <- is.na(label) | label == ""
  idx <- bad_label & !is.na(phylum)
  label[idx] <- paste0(phylum[idx], " taxon")

  bad_label <- is.na(label) | label == ""
  label[bad_label] <- as.character(tax_df[[tax_id_col]][bad_label])

  label
}

make_taxonomy_full <- function(tax_df) {
  rank_priority <- c(
    "Superkingdom", "Kingdom", "Domain",
    "Phylum", "Class", "Order", "Family", "Genus", "Species",
    "superkingdom", "kingdom", "domain",
    "phylum", "class", "order", "family", "genus", "species",
    "D", "P", "C", "O", "F", "G", "S",
    "Rank_D", "Rank_P", "Rank_C", "Rank_O", "Rank_F", "Rank_G", "Rank_S"
  )

  cols <- intersect(rank_priority, names(tax_df))

  if (length(cols) == 0) {
    return(rep(NA_character_, nrow(tax_df)))
  }

  apply(tax_df[, cols, drop = FALSE], 1, function(z) {
    z <- clean_tax_string(z)
    z <- z[!is.na(z) & nzchar(z)]
    if (length(z) == 0) return(NA_character_)
    paste(z, collapse = "; ")
  })
}

standardize_tax_table <- function(ps) {
  if (!is.null(phyloseq::tax_table(ps, errorIfNULL = FALSE))) {
    tax <- as.data.frame(phyloseq::tax_table(ps), stringsAsFactors = FALSE) %>%
      rownames_to_column("tax_id")
  } else {
    tax <- tibble(tax_id = taxa_names(ps))
  }

  tax <- tax %>%
    mutate(
      tax_id = as.character(tax_id),
      taxon = tax_id,
      tax_id_clean = make.names(tax_id, unique = FALSE),
      tax_id_clean2 = gsub(
        "^X(?=[0-9])",
        "",
        make.names(tax_id, unique = FALSE),
        perl = TRUE
      )
    )

  tax$tax_label <- make_tax_label(tax, tax_id_col = "tax_id")
  tax$taxonomy_full <- make_taxonomy_full(tax)

  tax
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

add_missing_col <- function(df, col, value = NA_character_) {
  if (!col %in% names(df)) {
    df[[col]] <- value
  }
  df
}

add_taxonomy <- function(ind, tax_map) {
  ind <- ind %>%
    mutate(
      taxon = as.character(taxon),
      taxon_clean = make.names(taxon, unique = FALSE),
      taxon_clean2 = gsub(
        "^X(?=[0-9])",
        "",
        make.names(taxon, unique = FALSE),
        perl = TRUE
      )
    )

  tax_small <- tax_map %>%
    select(
      tax_id,
      tax_id_clean,
      tax_id_clean2,
      tax_label,
      taxonomy_full,
      everything()
    )

  # Match 1: exact taxon == tax_id
  out <- ind %>%
    left_join(
      tax_small,
      by = c("taxon" = "tax_id"),
      suffix = c("", ".tax_exact")
    )

  out$taxonomy_match_method <- ifelse(!is.na(out$tax_label), "exact_tax_id", NA_character_)

  # Match 2: taxon_clean == tax_id_clean
  missing_idx <- which(is.na(out$tax_label))

  if (length(missing_idx) > 0) {
    m2 <- ind[missing_idx, , drop = FALSE] %>%
      left_join(
        tax_small,
        by = c("taxon_clean" = "tax_id_clean"),
        suffix = c("", ".tax_clean")
      )

    replace_cols <- setdiff(names(m2), names(ind))

    for (cc in replace_cols) {
      if (!cc %in% names(out)) out[[cc]] <- NA
      out[[cc]][missing_idx] <- m2[[cc]]
    }

    good <- !is.na(m2$tax_label)
    out$taxonomy_match_method[missing_idx[good]] <- "clean_name"
  }

  # Match 3: taxon_clean2 == tax_id_clean2
  missing_idx <- which(is.na(out$tax_label))

  if (length(missing_idx) > 0) {
    m3 <- ind[missing_idx, , drop = FALSE] %>%
      left_join(
        tax_small,
        by = c("taxon_clean2" = "tax_id_clean2"),
        suffix = c("", ".tax_clean2")
      )

    replace_cols <- setdiff(names(m3), names(ind))

    for (cc in replace_cols) {
      if (!cc %in% names(out)) out[[cc]] <- NA
      out[[cc]][missing_idx] <- m3[[cc]]
    }

    good <- !is.na(m3$tax_label)
    out$taxonomy_match_method[missing_idx[good]] <- "clean_name_no_X"
  }

  out <- out %>%
    mutate(
      taxonomy_match_method = ifelse(
        is.na(taxonomy_match_method),
        "unmatched",
        taxonomy_match_method
      ),
      tax_label = ifelse(is.na(tax_label) | tax_label == "", taxon, tax_label),
      taxonomy_full = ifelse(is.na(taxonomy_full) | taxonomy_full == "", tax_label, taxonomy_full)
    )

  out <- out %>%
    select(
      any_of(c(
        "gradient",
        "taxon",
        "tax_label",
        "taxonomy_full",
        "taxonomy_match_method",
        "obs_cp",
        "z_score",
        "maxgrp_raw",
        "direction",
        "purity",
        "reliability",
        "frequency",
        "cp_05",
        "cp_10",
        "cp_50",
        "cp_90",
        "cp_95",
        "is_reliable_95",
        "is_reliable_90"
      )),
      everything(),
      -any_of(c("taxon_clean", "taxon_clean2"))
    )

  out
}

# ==============================================================================
# Main
# ==============================================================================

RESULTS_DIR <- "/home/fjbalvino/Tipping_points/resultados_finales"
LATEST001 <- file.path(RESULTS_DIR, "LATEST_001_auditar_estructura_metadata.txt")
RUN001 <- if (file.exists(LATEST001)) trimws(readLines(LATEST001, warn = FALSE, n = 1L)) else NA_character_
DEFAULT_PS <- file.path(RUN001, "rds", "001_phyloseq_canon_51.rds")
LATEST401 <- file.path(RESULTS_DIR, "LATEST_401_TITAN2_global_taxa_bimodales_gradientes_ambientales.txt")
DEFAULT401 <- if (file.exists(LATEST401)) trimws(readLines(LATEST401, warn = FALSE, n = 1L)) else NA_character_
OUT401 <- get_arg("--out401", Sys.getenv("OUT401", unset = DEFAULT401))
PS_RDS <- get_arg("--ps", DEFAULT_PS)

if (is.na(OUT401) || !nzchar(OUT401) || !dir.exists(OUT401)) {
  stop("Usa --out401 /ruta/al/output/401", call. = FALSE)
}

stop_if_missing_file(PS_RDS, "phyloseq rds")

DIR_TABLES <- file.path(OUT401, "tables")
DIR_LOGS <- file.path(OUT401, "logs")

safe_mkdir(DIR_TABLES)
safe_mkdir(DIR_LOGS)

msg("OUT401: ", OUT401)
msg("PS_RDS: ", PS_RDS)

ps <- readRDS(PS_RDS)
if (phyloseq::nsamples(ps) != 51L) stop("403 requiere el phyloseq canonico de 51 muestras.", call. = FALSE)
tax_map <- standardize_tax_table(ps)

# Si existe la tabla de taxa usados del 031, la usamos para reforzar labels.
taxa_used_file <- file.path(DIR_TABLES, "401_nonviral_bimodal_taxa_used.csv")

if (file.exists(taxa_used_file)) {
  msg("Reading taxa-used table: ", taxa_used_file)

  taxa_used <- readr::read_csv(taxa_used_file, show_col_types = FALSE) %>%
    mutate(tax_id = as.character(tax_id))

  if (!"tax_label" %in% names(taxa_used)) {
    taxa_used$tax_label <- NA_character_
  }

  if (!"taxonomy_full" %in% names(taxa_used)) {
    taxa_used$taxonomy_full <- NA_character_
  }

  taxa_used_small <- taxa_used %>%
    select(tax_id, tax_label, taxonomy_full) %>%
    rename(
      tax_label_used = tax_label,
      taxonomy_full_used = taxonomy_full
    )

  tax_map <- tax_map %>%
    select(-any_of(c("tax_label_used", "taxonomy_full_used"))) %>%
    left_join(taxa_used_small, by = "tax_id") %>%
    mutate(
      tax_label = ifelse(
        !is.na(tax_label_used) & tax_label_used != "",
        tax_label_used,
        tax_label
      ),
      taxonomy_full = ifelse(
        !is.na(taxonomy_full_used) & taxonomy_full_used != "",
        taxonomy_full_used,
        taxonomy_full
      )
    ) %>%
    select(-any_of(c("tax_label_used", "taxonomy_full_used")))
}

readr::write_csv(
  tax_map,
  file.path(DIR_TABLES, "403_taxonomy_map_from_phyloseq.csv")
)

ind_files <- list.files(
  DIR_TABLES,
  pattern = "^401_titan_indval_.*\\.csv$",
  full.names = TRUE
)

ind_files <- ind_files[
  !grepl("402_|403_|all_gradients|reliable_summary|annotated", basename(ind_files))
]

if (length(ind_files) == 0) {
  stop("No encontré archivos 401_titan_indval_<gradient>.csv en: ", DIR_TABLES, call. = FALSE)
}

all_annotated <- list()

for (f in ind_files) {
  grad <- basename(f)
  grad <- gsub("^401_titan_indval_", "", grad)
  grad <- gsub("\\.csv$", "", grad)

  msg("Annotating: ", grad)

  ind <- readr::read_csv(f, show_col_types = FALSE)

  if (!"taxon" %in% names(ind)) {
    msg("WARNING: file without taxon column: ", f)
    next
  }

  if (!"gradient" %in% names(ind)) {
    ind <- ind %>% mutate(gradient = grad, .before = 1)
  }

  ind <- ind %>%
    mutate(
      taxon = as.character(taxon),
      obs_cp = if ("obs_cp" %in% names(.)) safe_numeric(obs_cp) else NA_real_,
      z_score = if ("z_score" %in% names(.)) safe_numeric(z_score) else NA_real_,
      purity = if ("purity" %in% names(.)) safe_numeric(purity) else NA_real_,
      reliability = if ("reliability" %in% names(.)) safe_numeric(reliability) else NA_real_,
      cp_05 = if ("cp_05" %in% names(.)) safe_numeric(cp_05) else NA_real_,
      cp_10 = if ("cp_10" %in% names(.)) safe_numeric(cp_10) else NA_real_,
      cp_50 = if ("cp_50" %in% names(.)) safe_numeric(cp_50) else NA_real_,
      cp_90 = if ("cp_90" %in% names(.)) safe_numeric(cp_90) else NA_real_,
      cp_95 = if ("cp_95" %in% names(.)) safe_numeric(cp_95) else NA_real_
    )

  ind <- add_missing_col(ind, "direction", NA_character_)
  ind <- add_missing_col(ind, "is_reliable_95", FALSE)
  ind <- add_missing_col(ind, "is_reliable_90", FALSE)

  ann <- add_taxonomy(ind, tax_map)

  out_file <- file.path(
    DIR_TABLES,
    paste0("403_titan_indval_", grad, "_annotated.csv")
  )

  readr::write_csv(ann, out_file)

  msg("Wrote: ", out_file)
  msg("Matched: ", sum(ann$taxonomy_match_method != "unmatched"), "/", nrow(ann))

  all_annotated[[grad]] <- ann
}

all_tbl <- bind_rows(all_annotated)

readr::write_csv(
  all_tbl,
  file.path(DIR_TABLES, "403_titan_indval_all_gradients_annotated.csv")
)

if ("is_reliable_95" %in% names(all_tbl)) {
  reliable95 <- all_tbl %>%
    filter(is_reliable_95 %in% TRUE) %>%
    arrange(gradient, direction, obs_cp, desc(abs(z_score)))

  readr::write_csv(
    reliable95,
    file.path(DIR_TABLES, "403_titan_indval_reliable95_annotated.csv")
  )
}

if ("is_reliable_90" %in% names(all_tbl)) {
  reliable90 <- all_tbl %>%
    filter(is_reliable_90 %in% TRUE) %>%
    arrange(gradient, direction, obs_cp, desc(abs(z_score)))

  readr::write_csv(
    reliable90,
    file.path(DIR_TABLES, "403_titan_indval_reliable90_annotated.csv")
  )
}

audit <- all_tbl %>%
  count(gradient, taxonomy_match_method, name = "n") %>%
  group_by(gradient) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  arrange(gradient, desc(n))

readr::write_csv(
  audit,
  file.path(DIR_TABLES, "403_titan_taxonomy_match_audit.csv")
)

summary_tax <- all_tbl %>%
  mutate(
    direction = ifelse(is.na(direction), "NA", as.character(direction)),
    tax_label = ifelse(is.na(tax_label), taxon, tax_label),
    taxonomy_full = ifelse(is.na(taxonomy_full), tax_label, taxonomy_full)
  ) %>%
  group_by(gradient, direction, tax_label, taxonomy_full) %>%
  summarise(
    n_rows = n(),
    best_abs_z = suppressWarnings(max(abs(z_score), na.rm = TRUE)),
    median_cp = suppressWarnings(median(obs_cp, na.rm = TRUE)),
    reliable95 = any(is_reliable_95 %in% TRUE, na.rm = TRUE),
    reliable90 = any(is_reliable_90 %in% TRUE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    best_abs_z = ifelse(is.infinite(best_abs_z), NA_real_, best_abs_z),
    median_cp = ifelse(is.nan(median_cp), NA_real_, median_cp)
  ) %>%
  arrange(gradient, direction, desc(reliable95), desc(reliable90), desc(best_abs_z))

readr::write_csv(
  summary_tax,
  file.path(DIR_TABLES, "403_titan_taxonomic_response_summary.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(DIR_LOGS, "403_sessionInfo.txt")
)

msg("Done.")
msg("Main annotated file: ", file.path(DIR_TABLES, "403_titan_indval_all_gradients_annotated.csv"))
msg("Reliable95 file: ", file.path(DIR_TABLES, "403_titan_indval_reliable95_annotated.csv"))
msg("Audit file: ", file.path(DIR_TABLES, "403_titan_taxonomy_match_audit.csv"))
