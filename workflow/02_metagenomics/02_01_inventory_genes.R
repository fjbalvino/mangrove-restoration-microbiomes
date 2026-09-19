#!/usr/bin/env Rscript
# ============================================================
# 02_01_inventory_genes.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Inventory of the already assembled/annotated abundance table; raw-read producer not recovered.
# Inputs (source expressions; complete list in docs/contracts/02_01_inventory_genes.json):
#   run003 <- readLines(latest003, warn = FALSE)[1]
#   header_line <- readLines(input_tsv, n = 1, warn = FALSE)
#   meta_raw <- read_csv(meta_csv, show_col_types = FALSE)
# Outputs (source expressions; complete list in contract):
#   write.table(
#   ggsave(out_file, p, width = 10, height = 5, dpi = 300, bg = "white")
#   write_tsv(
#   write_tsv(profile_gate, file.path(out_dir, "tables", "501_auditoria_canon_51.tsv"))
#   ggsave(
#   writeLines(out_dir, latest_file)
# Algorithmic provenance:
# Inventory the existing annotated gene-abundance table and sample columns.
#   Cantalapiedra et al. (2021), eggNOG-mapper v2, doi:10.1093/molbev/msab293; this is annotation provenance, not evidence that this script runs eggNOG.
# Source SHA-256: b5db501f8a2efbc9ea9cf1d33ea464e235e91029768c1800864369bceda779c4
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


# ==============================================================================
# 501_inventariar_universo_proteico.R
#
# Objetivo:
#   Inventariar una matriz masiva proteína/ORF x muestra sin cargarla completa
#   en memoria.
#
# Input esperado:
#   FINAL_ANNOTATED_ABUNDANCE_noMZ_bacteria_fran.tsv
#
# Estructura esperada:
#   columna 1 = #query / protein_id / ORF_id
#   columnas restantes = muestras
#
# Outputs:
#   Results/501_inventariar_universo_proteico_YYYYMMDD_HHMMSS/
#     tables/
#       501_inventory_summary.tsv
#       501_header_columns.tsv
#       501_sample_sums.tsv
#       501_sample_detection.tsv
#       501_prevalence_distribution.tsv
#       501_total_abundance_bins.tsv
#       501_origin_sample_feature_counts.tsv
#       501_top_features_by_total.tsv
#       501_feature_summary.tsv.gz
#       501_metadata_sample_matching.tsv
#     figures/
#       501_sample_total_abundance.png
#       501_sample_detected_features.png
#       501_prevalence_distribution.png
#       501_total_abundance_bins.png
#     logs/
#       501_log.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})

timestamp_now <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

message2 <- function(...) {
  message(sprintf(...))
}

stop2 <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

clean_colnames <- function(x) {
  x %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("[^A-Za-z0-9_\\.]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    tolower()
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  defaults <- list(
    input_tsv = "/home/fjbalvino/Tipping_points/FINAL_ANNOTATED_ABUNDANCE_noMZ_bacteria_fran.tsv",
    meta_csv = NA_character_,
    out_root = "/home/fjbalvino/Tipping_points/resultados_finales",
    chunk_size = "100000",
    write_feature_summary = "1",
    top_n = "5000"
  )

  if (length(args) == 0) {
    return(defaults)
  }

  if (length(args) %% 2 != 0) {
    stop2("Argumentos inválidos. Usa pares tipo --input_tsv VALUE --meta_csv VALUE.")
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

detect_origin_sample <- function(query_id) {
  query_id <- as.character(query_id)

  origin <- dplyr::case_when(
    str_detect(query_id, "\\.contigs_") ~ str_replace(query_id, "\\.contigs_.*$", ""),
    str_detect(query_id, "\\.contig_") ~ str_replace(query_id, "\\.contig_.*$", ""),
    str_detect(query_id, "_contigs_") ~ str_replace(query_id, "_contigs_.*$", ""),
    TRUE ~ str_replace(query_id, "\\..*$", "")
  )

  origin
}

row_sd_fast <- function(mat) {
  n <- ncol(mat)

  if (n <= 1) {
    return(rep(NA_real_, nrow(mat)))
  }

  rs <- rowSums(mat, na.rm = TRUE)
  rs2 <- rowSums(mat * mat, na.rm = TRUE)

  variance <- (rs2 - (rs * rs / n)) / (n - 1)
  variance[variance < 0 & variance > -1e-8] <- 0

  sqrt(variance)
}

bin_total_abundance <- function(x) {
  cut(
    x,
    breaks = c(-Inf, 0, 1, 5, 10, 20, 50, 100, 500, 1000, 10000, Inf),
    labels = c(
      "0",
      "1",
      "2-5",
      "6-10",
      "11-20",
      "21-50",
      "51-100",
      "101-500",
      "501-1000",
      "1001-10000",
      ">10000"
    ),
    right = TRUE
  )
}

write_tsv_connection <- function(df, con, col_names = FALSE) {
  write.table(
    df,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = col_names,
    na = ""
  )
}

make_basic_plot <- function(df, x_col, y_col, out_file, title, xlab, ylab) {
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_col(width = 0.8) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    ) +
    labs(
      title = title,
      x = xlab,
      y = ylab
    )

  ggsave(out_file, p, width = 10, height = 5, dpi = 300, bg = "white")
}

# ------------------------------------------------------------------------------
# Args
# ------------------------------------------------------------------------------

args <- parse_args()

input_tsv <- args$input_tsv
meta_csv <- args$meta_csv
out_root <- args$out_root
chunk_size <- as.integer(args$chunk_size)
write_feature_summary <- as.integer(args$write_feature_summary) == 1
top_n <- as.integer(args$top_n)

if (!file.exists(input_tsv)) {
  stop2("No existe input_tsv: %s", input_tsv)
}

if (is.na(meta_csv) || !nzchar(meta_csv)) {
  latest003 <- file.path(out_root, "LATEST_003_integrar_ejes_ECI_HI_en_metadata.txt")
  if (!file.exists(latest003)) stop2("No existe latest de metadata integrada: %s", latest003)
  run003 <- readLines(latest003, warn = FALSE)[1]
  meta_csv <- file.path(run003, "tables", "003_metadata_integrada_canon_51.csv")
}
if (!file.exists(meta_csv)) stop2("No existe meta_csv: %s", meta_csv)

if (!is.finite(chunk_size) || chunk_size < 1000) {
  stop2("chunk_size inválido: %s", args$chunk_size)
}

if (!is.finite(top_n) || top_n < 100) {
  stop2("top_n inválido: %s", args$top_n)
}

script_id <- "501"
script_base <- "inventariar_universo_proteico"

out_dir <- file.path(
  out_root,
  paste0(script_id, "_", script_base, "_", timestamp_now())
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "logs", "501_log.txt")
sink(log_file, split = TRUE)

on.exit({
  sink()
}, add = TRUE)

message2("===== 501 inventariar universo proteico =====")
message2("Start: %s", as.character(Sys.time()))
message2("Input: %s", input_tsv)
message2("Metadata: %s", meta_csv)
message2("Out dir: %s", out_dir)
message2("Chunk size: %s", chunk_size)
message2("Write feature summary: %s", write_feature_summary)
message2("Top N: %s", top_n)

# ------------------------------------------------------------------------------
# Header inspection
# ------------------------------------------------------------------------------

header_line <- readLines(input_tsv, n = 1, warn = FALSE)
header_cols <- strsplit(header_line, "\t", fixed = TRUE)[[1]]

if (length(header_cols) < 2) {
  stop2("El archivo no parece tabular o tiene menos de 2 columnas.")
}

query_col <- header_cols[1]
sample_cols <- header_cols[-1]

header_tbl <- tibble(
  column_index = seq_along(header_cols),
  column_name = header_cols,
  inferred_role = if_else(column_index == 1, "feature_id", "sample_abundance")
)

write_tsv(
  header_tbl,
  file.path(out_dir, "tables", "501_header_columns.tsv")
)

message2("Columns detected: %s", length(header_cols))
message2("Feature/query column: %s", query_col)
message2("Sample columns: %s", length(sample_cols))
message2("First sample columns: %s", paste(head(sample_cols, 10), collapse = ", "))

# ------------------------------------------------------------------------------
# Metadata matching
# ------------------------------------------------------------------------------

metadata_matching <- tibble()

if (!is.na(meta_csv) && nzchar(meta_csv) && file.exists(meta_csv)) {
  meta_raw <- read_csv(meta_csv, show_col_types = FALSE)
  original_meta_names <- names(meta_raw)
  names(meta_raw) <- clean_colnames(names(meta_raw))

  candidate_id_cols <- intersect(
    c("sample_id", ".sample_id", "id_metagenomics", "run_accession", "sample", "sample_name"),
    names(meta_raw)
  )

  if (length(candidate_id_cols) > 0) {
    metadata_matching <- map_dfr(candidate_id_cols, function(id_col) {
      meta_ids <- unique(as.character(meta_raw[[id_col]]))
      meta_ids <- meta_ids[!is.na(meta_ids) & nzchar(meta_ids)]

      tibble(
        metadata_id_column = id_col,
        n_metadata_ids = length(meta_ids),
        n_matrix_samples = length(sample_cols),
        n_common = length(intersect(sample_cols, meta_ids)),
        matrix_not_in_metadata = paste(setdiff(sample_cols, meta_ids), collapse = ";"),
        metadata_not_in_matrix = paste(setdiff(meta_ids, sample_cols), collapse = ";")
      )
    })

    write_tsv(
      metadata_matching,
      file.path(out_dir, "tables", "501_metadata_sample_matching.tsv")
    )

    message2("Metadata ID columns checked: %s", paste(candidate_id_cols, collapse = ", "))
    print(metadata_matching)

    best_match <- metadata_matching %>% arrange(desc(n_common)) %>% slice_head(n = 1L)
    best_id_col <- best_match$metadata_id_column[[1L]]
    canonical_ids <- trimws(as.character(meta_raw[[best_id_col]]))
    canonical_ids <- canonical_ids[!is.na(canonical_ids) & nzchar(canonical_ids)]
    if (length(canonical_ids) != 51L || anyDuplicated(canonical_ids) ||
        best_match$n_common[[1L]] != 51L) {
      stop2("501 requiere una metadata canonica de 51 IDs unicos, todos presentes en la matriz proteica.")
    }
    required_design <- c("lat_block", "depth_cm", "locality", "restoration4")
    if (length(setdiff(required_design, names(meta_raw)))) {
      stop2("501 requiere lat_block, depth_cm, locality y restoration4 en la metadata canonica.")
    }
    depth_num <- suppressWarnings(as.numeric(meta_raw$depth_cm))
    profile_key <- paste(meta_raw$locality, meta_raw$lat_block, sep = "::")
    profile_gate <- tibble(profile_id = profile_key, depth_cm = depth_num) %>%
      group_by(profile_id) %>%
      summarise(n = n(), valid = n == 3L && setequal(depth_cm, c(5, 20, 40)), .groups = "drop")
    write_tsv(profile_gate, file.path(out_dir, "tables", "501_auditoria_canon_51.tsv"))
    if (nrow(profile_gate) != 17L || !all(profile_gate$valid)) {
      stop2("501 requiere exactamente 17 perfiles completos 5/20/40.")
    }
  } else {
    metadata_matching <- tibble(
      note = "No candidate sample ID columns found in metadata",
      metadata_columns = paste(original_meta_names, collapse = ";")
    )

    write_tsv(
      metadata_matching,
      file.path(out_dir, "tables", "501_metadata_sample_matching.tsv")
    )

    message2("No candidate sample ID columns found in metadata.")
  }
} else {
  metadata_matching <- tibble(
    note = "metadata file not found or not provided",
    meta_csv = meta_csv
  )

  write_tsv(
    metadata_matching,
    file.path(out_dir, "tables", "501_metadata_sample_matching.tsv")
  )

  message2("Metadata file not found or not provided.")
}

# ------------------------------------------------------------------------------
# Streaming inventory
# ------------------------------------------------------------------------------

n_samples <- length(sample_cols)

col_sums <- setNames(rep(0, n_samples), sample_cols)
col_detected_features <- setNames(rep(0, n_samples), sample_cols)
col_max <- setNames(rep(0, n_samples), sample_cols)

prevalence_counts <- setNames(rep(0, n_samples + 1), as.character(0:n_samples))

total_bins_levels <- levels(bin_total_abundance(c(0, 1, 2, 6, 11, 21, 51, 101, 501, 1001, 10001)))
total_bins_counts <- setNames(rep(0, length(total_bins_levels)), total_bins_levels)

origin_counts <- new.env(parent = emptyenv())
top_accumulator <- tibble(
  feature_id = character(),
  origin_sample = character(),
  total_abundance = numeric(),
  prevalence = integer(),
  mean_abundance = numeric(),
  sd_abundance = numeric(),
  max_abundance = numeric()
)

n_features_total <- 0
n_features_nonzero <- 0
n_features_zero <- 0
n_negative_values <- 0
n_nonfinite_values <- 0
chunk_counter <- 0

feature_summary_file <- file.path(out_dir, "tables", "501_feature_summary.tsv.gz")
feature_con <- NULL

if (write_feature_summary) {
  feature_con <- gzfile(feature_summary_file, open = "wt")

  write_tsv_connection(
    tibble(
      feature_id = character(),
      origin_sample = character(),
      total_abundance = numeric(),
      prevalence = integer(),
      mean_abundance = numeric(),
      sd_abundance = numeric(),
      max_abundance = numeric()
    ),
    feature_con,
    col_names = TRUE
  )
}

col_type_string <- paste0("c", paste(rep("d", n_samples), collapse = ""))

callback_fun <- function(x, pos) {
  chunk_counter <<- chunk_counter + 1

  names(x)[1] <- query_col

  query_id <- as.character(x[[query_col]])

  mat <- as.matrix(x[, sample_cols])
  storage.mode(mat) <- "numeric"

  nonfinite_mask <- !is.finite(mat)
  n_nonfinite_values <<- n_nonfinite_values + sum(nonfinite_mask, na.rm = TRUE)

  if (any(nonfinite_mask, na.rm = TRUE)) {
    mat[nonfinite_mask] <- 0
  }

  n_negative_values <<- n_negative_values + sum(mat < 0, na.rm = TRUE)

  row_total <- rowSums(mat, na.rm = TRUE)
  row_prev <- rowSums(mat > 0, na.rm = TRUE)
  row_mean <- row_total / n_samples
  row_sd <- row_sd_fast(mat)
  row_max <- apply(mat, 1, max, na.rm = TRUE)

  n_chunk <- length(query_id)

  n_features_total <<- n_features_total + n_chunk
  n_features_nonzero <<- n_features_nonzero + sum(row_total > 0, na.rm = TRUE)
  n_features_zero <<- n_features_zero + sum(row_total == 0, na.rm = TRUE)

  col_sums <<- col_sums + colSums(mat, na.rm = TRUE)
  col_detected_features <<- col_detected_features + colSums(mat > 0, na.rm = TRUE)
  col_max <<- pmax(col_max, apply(mat, 2, max, na.rm = TRUE))

  prev_tab <- table(factor(row_prev, levels = 0:n_samples))
  prevalence_counts <<- prevalence_counts + as.numeric(prev_tab)

  total_bin <- bin_total_abundance(row_total)
  total_bin_tab <- table(factor(total_bin, levels = total_bins_levels))
  total_bins_counts <<- total_bins_counts + as.numeric(total_bin_tab)

  origin_sample <- detect_origin_sample(query_id)
  origin_tab <- table(origin_sample)

  for (nm in names(origin_tab)) {
    old <- if (exists(nm, envir = origin_counts, inherits = FALSE)) {
      get(nm, envir = origin_counts)
    } else {
      0
    }

    assign(nm, old + as.integer(origin_tab[[nm]]), envir = origin_counts)
  }

  feature_df <- tibble(
    feature_id = query_id,
    origin_sample = origin_sample,
    total_abundance = row_total,
    prevalence = as.integer(row_prev),
    mean_abundance = row_mean,
    sd_abundance = row_sd,
    max_abundance = row_max
  )

  if (write_feature_summary) {
    write_tsv_connection(feature_df, feature_con, col_names = FALSE)
  }

  top_take <- min(top_n, nrow(feature_df))

  chunk_top <- feature_df %>%
    arrange(desc(total_abundance), desc(prevalence)) %>%
    slice_head(n = top_take)

  top_accumulator <<- bind_rows(top_accumulator, chunk_top) %>%
    arrange(desc(total_abundance), desc(prevalence)) %>%
    distinct(feature_id, .keep_all = TRUE) %>%
    slice_head(n = top_n)

  if (chunk_counter %% 10 == 0) {
    message2(
      "Processed chunk %s | rows so far: %s | nonzero features: %s | current top retained: %s",
      chunk_counter,
      format(n_features_total, big.mark = ","),
      format(n_features_nonzero, big.mark = ","),
      nrow(top_accumulator)
    )
  }

  invisible()
}

read_tsv_chunked(
  file = input_tsv,
  callback = SideEffectChunkCallback$new(callback_fun),
  chunk_size = chunk_size,
  col_types = col_type_string,
  progress = TRUE
)

if (!is.null(feature_con)) {
  close(feature_con)
}

message2("Finished streaming matrix.")
message2("Total features: %s", format(n_features_total, big.mark = ","))
message2("Nonzero features: %s", format(n_features_nonzero, big.mark = ","))
message2("Zero-total features: %s", format(n_features_zero, big.mark = ","))

# ------------------------------------------------------------------------------
# Output summaries
# ------------------------------------------------------------------------------

sample_sums <- tibble(
  sample_id = names(col_sums),
  total_abundance = as.numeric(col_sums),
  detected_features = as.numeric(col_detected_features),
  max_abundance = as.numeric(col_max)
) %>%
  arrange(desc(total_abundance))

write_tsv(
  sample_sums,
  file.path(out_dir, "tables", "501_sample_sums.tsv")
)

sample_detection <- sample_sums %>%
  mutate(
    detected_feature_fraction = detected_features / n_features_total
  )

write_tsv(
  sample_detection,
  file.path(out_dir, "tables", "501_sample_detection.tsv")
)

prevalence_distribution <- tibble(
  prevalence = as.integer(names(prevalence_counts)),
  n_features = as.numeric(prevalence_counts),
  fraction_features = n_features / sum(n_features)
)

write_tsv(
  prevalence_distribution,
  file.path(out_dir, "tables", "501_prevalence_distribution.tsv")
)

total_abundance_bins <- tibble(
  total_abundance_bin = names(total_bins_counts),
  n_features = as.numeric(total_bins_counts),
  fraction_features = n_features / sum(n_features)
)

write_tsv(
  total_abundance_bins,
  file.path(out_dir, "tables", "501_total_abundance_bins.tsv")
)

origin_sample_feature_counts <- tibble(
  origin_sample = ls(origin_counts),
  n_features = map_int(ls(origin_counts), ~ get(.x, envir = origin_counts))
) %>%
  arrange(desc(n_features))

write_tsv(
  origin_sample_feature_counts,
  file.path(out_dir, "tables", "501_origin_sample_feature_counts.tsv")
)

top_features_by_total <- top_accumulator %>%
  arrange(desc(total_abundance), desc(prevalence)) %>%
  distinct(feature_id, .keep_all = TRUE) %>%
  slice_head(n = top_n)

write_tsv(
  top_features_by_total,
  file.path(out_dir, "tables", "501_top_features_by_total.tsv")
)

inventory_summary <- tibble(
  key = c(
    "script",
    "timestamp",
    "input_tsv",
    "input_size_bytes",
    "out_dir",
    "query_col",
    "n_columns_total",
    "n_sample_columns",
    "n_features_total",
    "n_features_nonzero",
    "n_features_zero_total",
    "n_negative_values",
    "n_nonfinite_values",
    "chunk_size",
    "write_feature_summary",
    "feature_summary_file",
    "top_n",
    "metadata_csv"
  ),
  value = c(
    "501_inventariar_universo_proteico.R",
    as.character(Sys.time()),
    normalizePath(input_tsv, mustWork = TRUE),
    as.character(file.info(input_tsv)$size),
    out_dir,
    query_col,
    as.character(length(header_cols)),
    as.character(n_samples),
    as.character(n_features_total),
    as.character(n_features_nonzero),
    as.character(n_features_zero),
    as.character(n_negative_values),
    as.character(n_nonfinite_values),
    as.character(chunk_size),
    as.character(write_feature_summary),
    ifelse(write_feature_summary, feature_summary_file, "not_written"),
    as.character(top_n),
    meta_csv
  )
)

write_tsv(
  inventory_summary,
  file.path(out_dir, "tables", "501_inventory_summary.tsv")
)

# ------------------------------------------------------------------------------
# Simple QC figures
# ------------------------------------------------------------------------------

make_basic_plot(
  sample_sums,
  x_col = "sample_id",
  y_col = "total_abundance",
  out_file = file.path(out_dir, "figures", "501_sample_total_abundance.png"),
  title = "Total protein/ORF abundance per sample",
  xlab = "Sample",
  ylab = "Total abundance"
)

make_basic_plot(
  sample_detection,
  x_col = "sample_id",
  y_col = "detected_features",
  out_file = file.path(out_dir, "figures", "501_sample_detected_features.png"),
  title = "Detected protein/ORF features per sample",
  xlab = "Sample",
  ylab = "Detected features"
)

p_prev <- ggplot(prevalence_distribution, aes(x = prevalence, y = n_features)) +
  geom_col(width = 0.8) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Protein/ORF prevalence distribution",
    x = "Number of samples where feature is detected",
    y = "Number of features"
  )

ggsave(
  file.path(out_dir, "figures", "501_prevalence_distribution.png"),
  p_prev,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

p_total_bins <- ggplot(total_abundance_bins, aes(x = total_abundance_bin, y = n_features)) +
  geom_col(width = 0.8) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Protein/ORF total abundance bins",
    x = "Total abundance bin",
    y = "Number of features"
  )

ggsave(
  file.path(out_dir, "figures", "501_total_abundance_bins.png"),
  p_total_bins,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------------------------
# Latest pointer
# ------------------------------------------------------------------------------

latest_file <- file.path(out_root, "LATEST_501_inventariar_universo_proteico.txt")
writeLines(out_dir, latest_file)

message2("===== DONE 501 =====")
message2("Output: %s", out_dir)
message2("Latest: %s", latest_file)
message2("Main tables:")
message2("  %s", file.path(out_dir, "tables", "501_inventory_summary.tsv"))
message2("  %s", file.path(out_dir, "tables", "501_metadata_sample_matching.tsv"))
message2("  %s", file.path(out_dir, "tables", "501_sample_sums.tsv"))
message2("  %s", file.path(out_dir, "tables", "501_prevalence_distribution.tsv"))
message2("  %s", file.path(out_dir, "tables", "501_total_abundance_bins.tsv"))
message2("  %s", file.path(out_dir, "tables", "501_origin_sample_feature_counts.tsv"))
message2("  %s", file.path(out_dir, "tables", "501_top_features_by_total.tsv"))
message2("  %s", file.path(out_dir, "tables", "501_feature_summary.tsv.gz"))
