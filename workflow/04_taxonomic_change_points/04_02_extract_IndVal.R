#!/usr/bin/env Rscript
# ============================================================
# 04_02_extract_IndVal.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Extracts complete indicator statistics; the 401 one-row indicator summary is not a taxon count.
# Inputs (source expressions; complete list in docs/contracts/04_02_extract_IndVal.json):
#   DEFAULT401 <- if (file.exists(LATEST401)) trimws(readLines(LATEST401, warn = FALSE, n = 1L)) else NA_character_
#   obj <- readRDS(rds)
# Outputs (source expressions; complete list in contract):
#   readr::write_csv(std, out_file)
#   readr::write_csv(
# Algorithmic provenance:
# Extract full taxon-level IndVal and sum-z thresholds from saved TITAN2 RDS objects.
#   Baker & King (2010), doi:10.1111/j.2041-210X.2009.00007.x.
# Source SHA-256: e5bfa1d62fbe36dccf033d9853177655ae211155df8b638d589047a5f4ab0853
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


# ==============================================================================
# 402_regenerate_TITAN2_indval_tables_from_rds.R
#
# Objetivo:
#   Regenerar correctamente los archivos:
#     401_titan_indval_<gradient>.csv
#
#   a partir de los objetos TITAN2 .rds ya generados por el script 401.
#
# Motivo:
#   TITAN2::titan.out() no está exportada en algunas versiones del paquete.
#   Este script intenta:
#     1) getFromNamespace("titan.out", "TITAN2")
#     2) obj$sppmax
#     3) obj$ivz / estructuras alternativas
#
# Outputs:
#   tables/
#     401_titan_indval_<gradient>.csv
#     402_titan_indval_all_gradients.csv
#     402_titan_indval_reliable_summary.csv
#     402_titan_object_structure_<gradient>.txt
# ==============================================================================

suppressPackageStartupMessages({
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

safe_mkdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
  cat("\n")
  flush.console()
}

sanitize_gradient <- function(x) {
  x <- gsub("^401_titan_", "", x)
  x <- gsub("\\.rds$", "", x)
  x
}

clean_names_soft <- function(x) {
  x <- as.character(x)
  x <- gsub("^X", "", x)
  x <- gsub("\\.+", ".", x)
  x <- gsub("\\s+", "_", x)
  x
}

find_col <- function(df, candidates) {
  nms <- names(df)
  hit <- intersect(candidates, nms)
  if (length(hit) > 0) return(hit[1])

  # relaxed matching
  nms_low <- tolower(gsub("[^a-z0-9]+", "", nms))
  cand_low <- tolower(gsub("[^a-z0-9]+", "", candidates))

  idx <- match(cand_low, nms_low)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) return(nms[idx[1]])

  NA_character_
}

num_or_na <- function(x) {
  suppressWarnings(as.numeric(x))
}

direction_from_maxgrp <- function(x) {
  x_chr <- tolower(as.character(x))

  case_when(
    x_chr %in% c("1", "z-", "-", "neg", "negative", "decrease", "decreasing") ~ "z-",
    x_chr %in% c("2", "z+", "+", "pos", "positive", "increase", "increasing") ~ "z+",
    TRUE ~ as.character(x)
  )
}

as_table_with_taxon <- function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  if (!"taxon" %in% names(x)) {
    rn <- rownames(x)
    if (!is.null(rn) && length(rn) == nrow(x) && !all(rn == seq_len(nrow(x)))) {
      x <- x %>% rownames_to_column("taxon")
    } else {
      x <- x %>% mutate(taxon = paste0("taxon_", row_number()), .before = 1)
    }
  }

  names(x) <- clean_names_soft(names(x))
  as_tibble(x)
}

try_internal_titan_out <- function(obj) {
  out <- tryCatch({
    fn <- getFromNamespace("titan.out", "TITAN2")
    as_table_with_taxon(fn(obj))
  }, error = function(e) {
    NULL
  })

  out
}

try_sppmax <- function(obj) {
  if (is.null(obj$sppmax)) return(NULL)

  x <- obj$sppmax

  # Some objects store sppmax as transposed taxa x metrics or metrics x taxa.
  df <- as.data.frame(x, stringsAsFactors = FALSE)

  # If there are very few rows and many columns, likely transposed.
  if (nrow(df) < ncol(df) && nrow(df) < 20) {
    df_t <- as.data.frame(t(df), stringsAsFactors = FALSE)
    df_t <- df_t %>% rownames_to_column("taxon")
    names(df_t) <- clean_names_soft(names(df_t))
    return(as_tibble(df_t))
  }

  out <- as_table_with_taxon(df)
  out
}

try_alternative_slots <- function(obj) {
  candidates <- c(
    "ivz",
    "indval",
    "taxa",
    "taxaout",
    "spp",
    "sppmax.boot",
    "maxspp"
  )

  for (cc in candidates) {
    if (!is.null(obj[[cc]])) {
      out <- tryCatch(as_table_with_taxon(obj[[cc]]), error = function(e) NULL)
      if (!is.null(out) && nrow(out) > 0) return(out)
    }
  }

  NULL
}

standardize_indval_table <- function(raw, gradient) {
  if (is.null(raw) || nrow(raw) == 0) {
    return(tibble(
      gradient = gradient,
      taxon = NA_character_,
      status = "empty_extraction"
    ))
  }

  names(raw) <- clean_names_soft(names(raw))

  tax_col <- find_col(raw, c("taxon", "species", "spp", "otu", "tax_id"))
  if (is.na(tax_col)) {
    raw <- raw %>% mutate(taxon = paste0("taxon_", row_number()), .before = 1)
    tax_col <- "taxon"
  }

  cp_col <- find_col(raw, c("zenv.cp", "zenv_cp",
    "obs.cp", "obs_cp", "env.cp", "env_cp", "cp", "change.point",
    "change_point", "env", "Env"
  ))

  z_col <- find_col(raw, c(
    "zscore", "z.score", "z_score", "obs.z", "obs_z", "z", "Z"
  ))

  maxgrp_col <- find_col(raw, c(
    "maxgrp", "max_grp", "grp", "group", "direction", "response"
  ))

  pur_col <- find_col(raw, c(
    "purity", "pur", "pure"
  ))

  rel_col <- find_col(raw, c(
    "reliability", "rel", "reliable"
  ))

  freq_col <- find_col(raw, c(
    "freq", "frequency"
  ))

  cp05_col <- find_col(raw, c("X5.","0.05", "X0.05", "q05", "cp.0.05", "cp05", "5.", "5%"))
  cp10_col <- find_col(raw, c("X10.","0.10", "X0.10", "q10", "cp.0.10", "cp10", "10.", "10%"))
  cp50_col <- find_col(raw, c("X50.","0.50", "X0.50", "q50", "cp.0.50", "cp50", "50.", "50%"))
  cp90_col <- find_col(raw, c("X90.","0.90", "X0.90", "q90", "cp.0.90", "cp90", "90.", "90%"))
  cp95_col <- find_col(raw, c("X95.","0.95", "X0.95", "q95", "cp.0.95", "cp95", "95.", "95%"))

  out <- raw %>%
    transmute(
      gradient = gradient,
      taxon = as.character(.data[[tax_col]]),
      obs_cp = if (!is.na(cp_col)) num_or_na(.data[[cp_col]]) else NA_real_,
      z_score = if (!is.na(z_col)) num_or_na(.data[[z_col]]) else NA_real_,
      maxgrp_raw = if (!is.na(maxgrp_col)) as.character(.data[[maxgrp_col]]) else NA_character_,
      direction = direction_from_maxgrp(maxgrp_raw),
      purity = if (!is.na(pur_col)) num_or_na(.data[[pur_col]]) else NA_real_,
      reliability = if (!is.na(rel_col)) num_or_na(.data[[rel_col]]) else NA_real_,
      frequency = if (!is.na(freq_col)) num_or_na(.data[[freq_col]]) else NA_real_,
      cp_05 = if (!is.na(cp05_col)) num_or_na(.data[[cp05_col]]) else NA_real_,
      cp_10 = if (!is.na(cp10_col)) num_or_na(.data[[cp10_col]]) else NA_real_,
      cp_50 = if (!is.na(cp50_col)) num_or_na(.data[[cp50_col]]) else NA_real_,
      cp_90 = if (!is.na(cp90_col)) num_or_na(.data[[cp90_col]]) else NA_real_,
      cp_95 = if (!is.na(cp95_col)) num_or_na(.data[[cp95_col]]) else NA_real_
    ) %>%
    mutate(
      is_reliable_95 = ifelse(
        is.finite(purity) & is.finite(reliability),
        purity >= 0.95 & reliability >= 0.95,
        FALSE
      ),
      is_reliable_90 = ifelse(
        is.finite(purity) & is.finite(reliability),
        purity >= 0.90 & reliability >= 0.90,
        FALSE
      )
    ) %>%
    arrange(direction, obs_cp, desc(abs(z_score)))

  # If cp quantiles are absent but obs_cp is present, keep the core table anyway.
  out
}

write_structure_file <- function(obj, gradient, file) {
  capture.output({
    cat("Gradient:", gradient, "\n\n")
    cat("Top-level names:\n")
    print(names(obj))
    cat("\nObject str(max.level=2):\n")
    str(obj, max.level = 2)
    cat("\nobj$sppmax names/dim preview:\n")
    if (!is.null(obj$sppmax)) {
      print(dim(obj$sppmax))
      print(head(as.data.frame(obj$sppmax)))
    } else {
      cat("obj$sppmax is NULL\n")
    }
  }, file = file)
}

extract_sumz_thresholds <- function(obj, gradient) {
  if (is.null(obj$sumz.cp)) stop("Falta sumz.cp en el RDS.", call. = FALSE)
  x <- as.data.frame(obj$sumz.cp, check.names = FALSE) %>% rownames_to_column("component")
  names(x) <- gsub("^X", "", names(x))
  get_num <- function(candidates) {
    hit <- intersect(candidates, names(x))
    if (!length(hit)) return(rep(NA_real_, nrow(x)))
    suppressWarnings(as.numeric(x[[hit[[1]]]]))
  }
  tibble(
    gradient = gradient,
    component = as.character(x$component),
    cp = get_num(c("cp", "CP")),
    q05 = get_num(c("0.05", "5.", "X0.05")),
    q10 = get_num(c("0.10", "10.", "X0.10")),
    q50 = get_num(c("0.50", "50.", "X0.50")),
    q90 = get_num(c("0.90", "90.", "X0.90")),
    q95 = get_num(c("0.95", "95.", "X0.95"))
  )
}

# ==============================================================================
# Main
# ==============================================================================

RESULTS_DIR <- "/home/fjbalvino/Tipping_points/resultados_finales"
LATEST401 <- file.path(RESULTS_DIR, "LATEST_401_TITAN2_global_taxa_bimodales_gradientes_ambientales.txt")
DEFAULT401 <- if (file.exists(LATEST401)) trimws(readLines(LATEST401, warn = FALSE, n = 1L)) else NA_character_
OUT401 <- get_arg("--out401", Sys.getenv("OUT401", unset = DEFAULT401))

if (is.na(OUT401) || !nzchar(OUT401) || !dir.exists(OUT401)) {
  stop("Usa --out401 /ruta/al/output/401", call. = FALSE)
}

DIR_OBJECTS <- file.path(OUT401, "objects")
DIR_TABLES <- file.path(OUT401, "tables")
DIR_LOGS <- file.path(OUT401, "logs")

safe_mkdir(DIR_TABLES)
safe_mkdir(DIR_LOGS)

rds_files <- list.files(
  DIR_OBJECTS,
  pattern = "^401_titan_.*\\.rds$",
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("No encontré objetos TITAN2 en: ", DIR_OBJECTS, call. = FALSE)
}

all_tables <- list()
all_sumz <- list()

for (rds in rds_files) {
  gradient <- sanitize_gradient(basename(rds))

  msg("Processing gradient: ", gradient)
  msg("Reading: ", rds)

  obj <- readRDS(rds)
  all_sumz[[gradient]] <- extract_sumz_thresholds(obj, gradient)

  write_structure_file(
    obj,
    gradient,
    file.path(DIR_LOGS, paste0("402_titan_object_structure_", gradient, ".txt"))
  )

  if (is.null(obj$sppmax) || nrow(obj$sppmax) != 1031L ||
      anyDuplicated(rownames(obj$sppmax))) {
    stop("El RDS no contiene 1031 taxa unicos en sppmax.", call. = FALSE)
  }
  raw <- as.data.frame(obj$sppmax, check.names = FALSE) %>%
    rownames_to_column("taxon")
  extraction_method <- "obj$sppmax"

  if (is.null(raw) || nrow(raw) == 0) {
    raw <- try_sppmax(obj)
    extraction_method <- "obj$sppmax"
  }

  if (is.null(raw) || nrow(raw) == 0) {
    raw <- try_alternative_slots(obj)
    extraction_method <- "alternative_slot"
  }

  if (is.null(raw) || nrow(raw) == 0) {
    msg("WARNING: Could not extract indval table for ", gradient)
    std <- tibble(
      gradient = gradient,
      taxon = NA_character_,
      status = "failed_extraction"
    )
  } else {
    std <- standardize_indval_table(raw, gradient) %>%
      mutate(extraction_method = extraction_method, .after = gradient)
  }

  if (nrow(std) != 1031L ||
      any(!is.finite(std$obs_cp)) ||
      any(!is.finite(std$cp_05)) ||
      any(!is.finite(std$cp_95)) ||
      anyNA(std$direction)) {
    stop("Extraccion incompleta de IDs, puntos, intervalos o direccion.",
         call. = FALSE)
  }

  out_file <- file.path(DIR_TABLES, paste0("401_titan_indval_", gradient, ".csv"))
  if (file.exists(out_file)) {
    file.copy(out_file, paste0(out_file, ".bak_",
      format(Sys.time(), "%Y%m%d_%H%M%S")), overwrite = FALSE)
  }

  readr::write_csv(std, out_file)

  msg("Wrote: ", out_file)
  msg("Rows: ", nrow(std))

  all_tables[[gradient]] <- std
}

all_tbl <- bind_rows(all_tables)
sumz_tbl <- bind_rows(all_sumz)

readr::write_csv(
  all_tbl,
  file.path(DIR_TABLES, "402_titan_indval_all_gradients.csv")
)

readr::write_csv(
  sumz_tbl,
  file.path(DIR_TABLES, "402_titan_sumz_thresholds_all_gradients.csv")
)

summary_tbl <- all_tbl %>%
  filter(!is.na(taxon)) %>%
  group_by(gradient, direction) %>%
  summarise(
    n_taxa = n(),
    n_reliable_95 = sum(is_reliable_95, na.rm = TRUE),
    n_reliable_90 = sum(is_reliable_90, na.rm = TRUE),
    median_cp = median(obs_cp, na.rm = TRUE),
    min_cp = min(obs_cp, na.rm = TRUE),
    max_cp = max(obs_cp, na.rm = TRUE),
    median_abs_z = median(abs(z_score), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(gradient, direction)

readr::write_csv(
  summary_tbl,
  file.path(DIR_TABLES, "402_titan_indval_reliable_summary.csv")
)

msg("Done.")
msg("Combined table: ", file.path(DIR_TABLES, "402_titan_indval_all_gradients.csv"))
msg("Summary table: ", file.path(DIR_TABLES, "402_titan_indval_reliable_summary.csv"))
msg("sum(z) thresholds: ", file.path(DIR_TABLES, "402_titan_sumz_thresholds_all_gradients.csv"))
