#!/usr/bin/env Rscript
# ============================================================
# 03_01_distributional_screen.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Contains raw and residual screens. Only raw CLR joint criteria define the current 1,031 taxa; residual/bootstrap candidates are sensitivity results.
# Inputs (source expressions; complete list in docs/contracts/03_01_distributional_screen.json):
#   out <- trimws(readLines(f, n = 1L, warn = FALSE))
#   ps <- readRDS(PS_RDS)
#   ext_meta <- read_csv(META_CSV, show_col_types = FALSE)
# Outputs (source expressions; complete list in contract):
#   writeLines(OUT_DIR, LATEST_ATTEMPT_FILE)
#   write_csv(parameter_gate, file.path(TABLE_DIR, "201_compuerta_parametros.csv"))
#   write_csv(input_gate, file.path(TABLE_DIR, "201_compuerta_inputs.csv"))
#   write_csv(profile_gate, file.path(TABLE_DIR, "201_auditoria_perfiles_bootstrap.csv"))
#   write_csv(design_gate, file.path(TABLE_DIR, "201_compuerta_diseno.csv"))
#   write_csv(
#   write_csv(filter_audit, file.path(TABLE_DIR, "201_auditoria_filtrado_taxa.csv"))
#   write_csv(tax, file.path(TABLE_DIR, "201_mapa_taxonomia_especies.csv"))
#   saveRDS(
#   write_csv(final_gate, file.path(TABLE_DIR, "201_compuerta_seleccion_bimodal.csv"))
# Algorithmic provenance:
# CLR transformation; Sarle bimodality and Hartigan dip tests with BH adjustment; sensitivity screens.
#   McArdle & Anderson (2001), doi:10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2; Peres-Neto et al. (2006), doi:10.1890/0012-9658(2006)87[2614:VPOESD]2.0.CO;2.
# Source SHA-256: 01b56607a94a4d2eab758e274d55cbebcf46113879d64e8692e252f17e3773bc
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


# 201_detectar_especies_bimodales_bootstrap.R
# Descubrimiento exploratorio y conservador de especies con distribuciones
# compatibles con multimodalidad/bimodalidad.
#
# Principios:
#   1. La etapa de restauracion NO participa en el filtrado ni en la seleccion.
#   2. La matriz primaria se residualiza por localidad + profundidad.
#   3. locality x depth se usa como sensibilidad para gradientes verticales
#      distintos entre localidades; no redefine por si sola la seleccion.
#   4. El bootstrap remuestrea perfiles completos dentro de localidad.
#   5. La seleccion observada y el soporte bootstrap primario usan Dip con BH.
#   6. La bimodalidad transversal identifica centinelas candidatos; no prueba
#      estados estables alternativos, causalidad ni histeresis.

suppressPackageStartupMessages({
  library(phyloseq)
  library(diptest)
  library(moments)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# Utilidades
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (!length(i)) return(default)
  if (i == length(args)) stop("Argumento sin valor: ", flag, call. = FALSE)
  args[[i + 1L]]
}

parse_int_arg <- function(x, flag) {
  z <- suppressWarnings(as.integer(x))
  if (length(z) != 1L || is.na(z)) {
    stop("Valor entero invalido para ", flag, ": ", x, call. = FALSE)
  }
  z
}

parse_num_arg <- function(x, flag) {
  z <- suppressWarnings(as.numeric(x))
  if (length(z) != 1L || !is.finite(z)) {
    stop("Valor numerico invalido para ", flag, ": ", x, call. = FALSE)
  }
  z
}

stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S", tz = "UTC")

read_latest <- function(root, id) {
  f <- file.path(root, paste0("LATEST_", id, ".txt"))
  if (!file.exists(f)) {
    stop("No existe ", f, "; ejecute primero 001", call. = FALSE)
  }
  out <- trimws(readLines(f, n = 1L, warn = FALSE))
  if (!nzchar(out) || !dir.exists(out)) {
    stop("LATEST invalido o inexistente: ", f, call. = FALSE)
  }
  out
}

safe_bc <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4L || length(unique(x)) < 4L) return(NA_real_)
  sk <- suppressWarnings(moments::skewness(x))
  ku <- suppressWarnings(moments::kurtosis(x))
  if (!is.finite(sk) || !is.finite(ku) || ku <= 0) return(NA_real_)
  as.numeric((sk^2 + 1) / ku)
}

safe_dip <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4L || length(unique(x)) < 4L) {
    return(c(stat = NA_real_, p = NA_real_))
  }
  z <- tryCatch(diptest::dip.test(x), error = function(e) NULL)
  if (is.null(z)) return(c(stat = NA_real_, p = NA_real_))
  c(stat = unname(z$statistic), p = z$p.value)
}

metrics_matrix <- function(m) {
  if (!is.matrix(m) || ncol(m) < 1L) {
    stop("metrics_matrix requiere una matriz con al menos un taxon", call. = FALSE)
  }
  dip <- t(vapply(seq_len(ncol(m)), function(j) safe_dip(m[, j]), numeric(2)))
  tibble(
    Original_Taxon_ID = colnames(m),
    N = colSums(is.finite(m)),
    SD_CLR = apply(m, 2, sd, na.rm = TRUE),
    Sarle_BC = vapply(seq_len(ncol(m)), function(j) safe_bc(m[, j]), numeric(1)),
    Dip_statistic = dip[, "stat"],
    Dip_P = dip[, "p"]
  ) %>%
    mutate(Dip_Q_BH = p.adjust(Dip_P, method = "BH"))
}

normalize_depth <- function(x) {
  suppressWarnings(
    as.numeric(str_extract(as.character(x), "[0-9]+(?:\\.[0-9]+)?"))
  )
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

residualize <- function(m, meta, interaction = FALSE) {
  if (nrow(m) != nrow(meta)) {
    stop("Matriz y metadata difieren en numero de filas", call. = FALSE)
  }
  f <- if (isTRUE(interaction)) {
    ~ locality * factor(depth_cm)
  } else {
    ~ locality + factor(depth_cm)
  }
  X <- model.matrix(f, data = meta)
  qx <- qr(X)
  if (qx$rank < ncol(X)) {
    stop(
      "Matriz de residualizacion sin rango completo: ",
      paste(deparse(f), collapse = " "),
      call. = FALSE
    )
  }
  out <- qr.resid(qx, m)
  colnames(out) <- colnames(m)
  out
}

bootstrap_profile_indices <- function(meta) {
  unlist(
    lapply(sort(unique(meta$locality)), function(loc) {
      profiles <- unique(meta$profile_id[meta$locality == loc])
      sampled <- sample(profiles, length(profiles), replace = TRUE)
      unlist(
        lapply(sampled, function(pr) {
          idx <- which(meta$profile_id == pr)
          idx[order(meta$depth_cm[idx])]
        }),
        use.names = FALSE
      )
    }),
    use.names = FALSE
  )
}

wilson_interval <- function(k, n, conf_level = 0.95) {
  if (!is.finite(k) || !is.finite(n) || n <= 0L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  z <- qnorm(1 - (1 - conf_level) / 2)
  phat <- k / n
  denom <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denom
  half <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / denom
  c(lower = max(0, center - half), upper = min(1, center + half))
}

# ------------------------------------------------------------------------------
# Parametros y rutas
# ------------------------------------------------------------------------------

OUT_ROOT <- get_arg(
  "--out_root",
  "/home/fjbalvino/Tipping_points/resultados_finales"
)
RUN001 <- get_arg(
  "--run001",
  read_latest(OUT_ROOT, "001_auditar_estructura_metadata")
)
PS_RDS <- get_arg(
  "--ps_rds",
  file.path(RUN001, "rds", "001_phyloseq_canon_51.rds")
)
META_CSV <- get_arg(
  "--meta_csv",
  file.path(RUN001, "tables", "001_metadata_canon_51.csv")
)
PROFILE_COL <- get_arg("--profile_col", "lat_block")
LOCALITY_COL <- get_arg("--locality_col", "locality")
DEPTH_COL <- get_arg("--depth_col", "depth_cm")
SEED <- parse_int_arg(get_arg("--seed", "123"), "--seed")
N_BOOT <- parse_int_arg(get_arg("--n_boot", "1000"), "--n_boot")
PSEUDOCOUNT <- parse_num_arg(get_arg("--pseudocount", "1"), "--pseudocount")
MIN_PREV_N <- parse_int_arg(get_arg("--min_prev_n", "5"), "--min_prev_n")
MIN_TOTAL <- parse_num_arg(get_arg("--min_total", "20"), "--min_total")
MIN_SD_CLR <- parse_num_arg(get_arg("--min_sd_clr", "0.05"), "--min_sd_clr")
BC_CUTOFF <- parse_num_arg(get_arg("--bc_cutoff", "0.555"), "--bc_cutoff")
DIP_FDR <- parse_num_arg(get_arg("--dip_fdr", "0.05"), "--dip_fdr")
STABILITY_CUTOFF <- parse_num_arg(
  get_arg("--stability_cutoff", "0.80"),
  "--stability_cutoff"
)

SCRIPT_ID <- "201_detectar_especies_bimodales_bootstrap"
OUT_DIR <- file.path(OUT_ROOT, paste0(SCRIPT_ID, "_", stamp()))
TABLE_DIR <- file.path(OUT_DIR, "tables")
PLOT_DIR <- file.path(OUT_DIR, "plots")
LOG_DIR <- file.path(OUT_DIR, "logs")
RDS_DIR <- file.path(OUT_DIR, "rds")
LATEST_FILE <- file.path(OUT_ROOT, paste0("LATEST_", SCRIPT_ID, ".txt"))
LATEST_ATTEMPT_FILE <- file.path(
  OUT_ROOT,
  paste0("LATEST_ATTEMPT_", SCRIPT_ID, ".txt")
)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
writeLines(OUT_DIR, LATEST_ATTEMPT_FILE)

log_con <- file(file.path(LOG_DIR, "201_log.txt"), "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0L) sink(type = "message")
  while (sink.number() > 0L) sink()
  close(log_con)
}, add = TRUE)

parameter_gate <- tibble(
  gate = c(
    "n_boot_at_least_100",
    "pseudocount_positive",
    "min_prev_n_valid_for_canon51",
    "min_total_nonnegative",
    "min_sd_clr_nonnegative",
    "bc_cutoff_in_unit_interval",
    "dip_fdr_in_unit_interval",
    "stability_cutoff_in_unit_interval"
  ),
  pass = c(
    N_BOOT >= 100L,
    PSEUDOCOUNT > 0,
    MIN_PREV_N >= 1L && MIN_PREV_N <= 51L,
    MIN_TOTAL >= 0,
    MIN_SD_CLR >= 0,
    BC_CUTOFF > 0 && BC_CUTOFF < 1,
    DIP_FDR > 0 && DIP_FDR <= 1,
    STABILITY_CUTOFF >= 0 && STABILITY_CUTOFF <= 1
  ),
  detail = c(
    as.character(N_BOOT),
    as.character(PSEUDOCOUNT),
    as.character(MIN_PREV_N),
    as.character(MIN_TOTAL),
    as.character(MIN_SD_CLR),
    as.character(BC_CUTOFF),
    as.character(DIP_FDR),
    as.character(STABILITY_CUTOFF)
  )
)
write_csv(parameter_gate, file.path(TABLE_DIR, "201_compuerta_parametros.csv"))
if (!all(parameter_gate$pass)) {
  stop("Fallo la compuerta de parametros de 201", call. = FALSE)
}

if (!file.exists(PS_RDS)) {
  stop("No existe phyloseq canonico de 51 muestras: ", PS_RDS, call. = FALSE)
}
if (!file.exists(META_CSV)) {
  stop("No existe metadata canonica de 51 muestras: ", META_CSV, call. = FALSE)
}

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)

cat("===== 201 Conservative bimodality discovery =====\n")
cat("Start UTC:", format(Sys.time(), tz = "UTC"), "\n")
cat("Phyloseq:", PS_RDS, "\n")
cat("Metadata:", META_CSV, "\n")
cat("Output:", OUT_DIR, "\n")
cat("Bootstrap iterations:", N_BOOT, "\n")

# ------------------------------------------------------------------------------
# Lectura, alineacion y compuerta del canon
# ------------------------------------------------------------------------------

ps <- readRDS(PS_RDS)
otu <- as(phyloseq::otu_table(ps), "matrix")
if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)
storage.mode(otu) <- "numeric"

ps_meta_df <- as.data.frame(
  phyloseq::sample_data(ps),
  stringsAsFactors = FALSE
)
# as.data.frame() puede conservar la clase S3 "sample_data" de phyloseq.
# Se fuerza una data.frame base antes de usar cualquier operación tabular.
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

# Los nombres de fila de sample_data son los identificadores canónicos de
# phyloseq y deben alinearse con las columnas de la OTU. Si sample_data ya
# contiene sample_id, se audita su concordancia antes de retirar la copia para
# evitar que rownames_to_column() cree un nombre duplicado.
if (ps_meta_has_sample_id) {
  ps_meta_df[["sample_id"]] <- NULL
}
ps_meta <- data.frame(
  sample_id = ps_meta_row_ids,
  ps_meta_df,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = NULL
)

ext_meta <- read_csv(META_CSV, show_col_types = FALSE)
if (!"sample_id" %in% names(ext_meta)) {
  stop("meta_csv requiere sample_id", call. = FALSE)
}
ext_meta <- ext_meta %>% mutate(sample_id = as.character(sample_id))

input_gate <- tibble(
  gate = c(
    "phyloseq_has_51_samples",
    "sample_data_existing_id_agrees_with_rownames",
    "otu_and_sample_data_ids_match",
    "external_metadata_has_51_unique_ids",
    "phyloseq_and_external_metadata_ids_match",
    "taxon_ids_unique",
    "otu_values_finite",
    "otu_values_nonnegative",
    "otu_values_integer_like",
    "all_samples_have_positive_library_size"
  ),
  pass = c(
    isTRUE(ncol(otu) == 51L && nrow(ps_meta) == 51L),
    isTRUE(ps_meta_id_agreement),
    isTRUE(!anyDuplicated(ps_meta$sample_id) &&
      !anyDuplicated(colnames(otu)) &&
      setequal(ps_meta$sample_id, colnames(otu))),
    isTRUE(nrow(ext_meta) == 51L &&
      !anyNA(ext_meta$sample_id) &&
      !anyDuplicated(ext_meta$sample_id)),
    isTRUE(setequal(ps_meta$sample_id, ext_meta$sample_id)),
    isTRUE(!is.null(rownames(otu)) && !anyDuplicated(rownames(otu))),
    isTRUE(all(is.finite(otu))),
    isTRUE(all(is.finite(otu)) && all(otu >= 0)),
    isTRUE(all(is.finite(otu)) && all(abs(otu - round(otu)) < 1e-8)),
    isTRUE(all(is.finite(otu)) && all(colSums(otu) > 0))
  ),
  detail = c(
    paste0("otu=", ncol(otu), "; sample_data=", nrow(ps_meta)),
    if (ps_meta_has_sample_id) {
      paste0("column_present; exact_ordered_match=", ps_meta_id_agreement)
    } else {
      "column_absent; phyloseq rownames used"
    },
    paste0("otu_ids=", length(colnames(otu)), "; metadata_ids=", nrow(ps_meta)),
    paste0("n=", nrow(ext_meta), "; duplicated=", anyDuplicated(ext_meta$sample_id)),
    paste0("shared=", length(intersect(ps_meta$sample_id, ext_meta$sample_id))),
    paste0("n_taxa=", nrow(otu)),
    paste0("nonfinite=", sum(!is.finite(otu))),
    paste0("negative=", sum(otu < 0, na.rm = TRUE)),
    "raw count matrix required",
    paste0("empty_samples=", sum(colSums(otu) <= 0))
  )
)
write_csv(input_gate, file.path(TABLE_DIR, "201_compuerta_inputs.csv"))
if (!all(input_gate$pass)) {
  stop("Fallo la compuerta de inputs de 201", call. = FALSE)
}

meta0 <- inner_join(
  ps_meta,
  ext_meta,
  by = "sample_id",
  suffix = c(".ps", ".ext")
)

meta <- meta0 %>%
  mutate(
    locality = str_squish(
      coalesce_col(
        meta0,
        c(
          LOCALITY_COL,
          paste0(LOCALITY_COL, ".ext"),
          paste0(LOCALITY_COL, ".ps")
        )
      )
    ),
    profile_raw = str_squish(
      coalesce_col(
        meta0,
        c(
          PROFILE_COL,
          paste0(PROFILE_COL, ".ext"),
          paste0(PROFILE_COL, ".ps")
        )
      )
    ),
    depth_cm = normalize_depth(
      coalesce_col(
        meta0,
        c(
          DEPTH_COL,
          paste0(DEPTH_COL, ".ext"),
          paste0(DEPTH_COL, ".ps")
        )
      )
    ),
    profile_id = if_else(
      str_detect(profile_raw, fixed(locality)),
      profile_raw,
      paste(locality, profile_raw, sep = "::")
    )
  ) %>%
  filter(sample_id %in% colnames(otu)) %>%
  arrange(locality, profile_id, depth_cm)

if (anyNA(meta[, c("sample_id", "locality", "profile_id", "depth_cm")])) {
  stop("Metadata incompleta para sample_id/locality/profile/depth", call. = FALSE)
}
if (anyDuplicated(meta$sample_id)) {
  stop("sample_id duplicado despues de integrar metadata", call. = FALSE)
}

profile_gate <- meta %>%
  group_by(profile_id) %>%
  summarise(
    n = n(),
    n_localities = n_distinct(locality),
    locality_values = paste(sort(unique(locality)), collapse = "|"),
    depths = paste(sort(unique(depth_cm)), collapse = "|"),
    valid = n == 3L &&
      setequal(depth_cm, c(5, 20, 40)) &&
      n_localities == 1L,
    .groups = "drop"
  )
write_csv(profile_gate, file.path(TABLE_DIR, "201_auditoria_perfiles_bootstrap.csv"))

X_additive <- model.matrix(~ locality + factor(depth_cm), data = meta)
X_interaction <- model.matrix(~ locality * factor(depth_cm), data = meta)

design_gate <- tibble(
  gate = c(
    "canonical_51_samples",
    "canonical_17_complete_profiles",
    "canonical_three_depths",
    "canonical_three_localities",
    "residualization_additive_full_rank",
    "residualization_interaction_full_rank"
  ),
  pass = c(
    nrow(meta) == 51L && n_distinct(meta$sample_id) == 51L,
    nrow(profile_gate) == 17L && all(profile_gate$valid),
    setequal(meta$depth_cm, c(5, 20, 40)),
    n_distinct(meta$locality) == 3L,
    qr(X_additive)$rank == ncol(X_additive),
    qr(X_interaction)$rank == ncol(X_interaction)
  ),
  detail = c(
    paste0("n=", nrow(meta)),
    paste0("profiles=", nrow(profile_gate), "; valid=", sum(profile_gate$valid)),
    paste(sort(unique(meta$depth_cm)), collapse = ","),
    paste(sort(unique(meta$locality)), collapse = ","),
    "~ locality + factor(depth_cm)",
    "~ locality * factor(depth_cm)"
  )
)
write_csv(design_gate, file.path(TABLE_DIR, "201_compuerta_diseno.csv"))
if (!all(design_gate$pass)) {
  stop("Fallo la compuerta de diseno de 201", call. = FALSE)
}

write_csv(
  tibble(
    role = c("primary", "sensitivity"),
    residualization_formula = c(
      "CLR ~ locality + factor(depth_cm)",
      "CLR ~ locality * factor(depth_cm)"
    ),
    restoration_used = c(FALSE, FALSE),
    interpretation = c(
      "primary candidate discovery and profile bootstrap",
      "diagnostic for locality-specific vertical structure"
    )
  ),
  file.path(TABLE_DIR, "201_modelos_residualizacion.csv")
)

otu <- otu[, meta$sample_id, drop = FALSE]

# ------------------------------------------------------------------------------
# Filtrado independiente del resultado y transformacion CLR
# ------------------------------------------------------------------------------

prev <- rowSums(otu > 0)
total <- rowSums(otu)
keep_abundance <- prev >= MIN_PREV_N & total >= MIN_TOTAL
if (sum(keep_abundance) < 10L) {
  stop("Menos de 10 taxa pasaron prevalencia y abundancia", call. = FALSE)
}

otu_abundance <- otu[keep_abundance, , drop = FALSE]
log_counts <- log(t(otu_abundance) + PSEUDOCOUNT)
clr_abundance <- sweep(log_counts, 1, rowMeans(log_counts), "-")
sd_clr <- apply(clr_abundance, 2, sd)
keep_sd <- is.finite(sd_clr) & sd_clr >= MIN_SD_CLR
if (sum(keep_sd) < 10L) {
  stop("Menos de 10 taxa pasaron el filtro de variacion CLR", call. = FALSE)
}

filter_audit <- tibble(
  Original_Taxon_ID = rownames(otu),
  prevalence_n = as.integer(prev),
  total_count = as.numeric(total),
  pass_prevalence_total = keep_abundance,
  SD_CLR_after_abundance_filter = NA_real_,
  pass_SD_CLR = FALSE,
  retained_for_screen = FALSE
)
filter_audit$SD_CLR_after_abundance_filter[
  match(colnames(clr_abundance), filter_audit$Original_Taxon_ID)
] <- sd_clr
filter_audit$pass_SD_CLR[
  match(colnames(clr_abundance), filter_audit$Original_Taxon_ID)
] <- keep_sd
filter_audit$retained_for_screen <-
  filter_audit$pass_prevalence_total & filter_audit$pass_SD_CLR
write_csv(filter_audit, file.path(TABLE_DIR, "201_auditoria_filtrado_taxa.csv"))

otu_f <- otu_abundance[keep_sd, , drop = FALSE]
clr <- clr_abundance[, keep_sd, drop = FALSE]
resid_clr <- residualize(clr, meta, interaction = FALSE)
resid_clr_interaction <- residualize(clr, meta, interaction = TRUE)

# ------------------------------------------------------------------------------
# Metricas observadas y taxonomia
# ------------------------------------------------------------------------------

raw_metrics <- metrics_matrix(clr) %>%
  rename_with(~ paste0("raw_", .x), -Original_Taxon_ID)
resid_metrics <- metrics_matrix(resid_clr) %>%
  rename_with(~ paste0("resid_", .x), -Original_Taxon_ID)
interaction_metrics <- metrics_matrix(resid_clr_interaction) %>%
  rename_with(~ paste0("resid_interaction_", .x), -Original_Taxon_ID)

tax_obj <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
if (is.null(tax_obj)) {
  stop("El phyloseq no contiene tax_table", call. = FALSE)
}
tax <- as.data.frame(tax_obj, stringsAsFactors = FALSE) %>%
  rownames_to_column("Original_Taxon_ID") %>%
  as_tibble() %>%
  arrange(Original_Taxon_ID)
rank_cols <- setdiff(names(tax), "Original_Taxon_ID")
tax <- tax %>%
  mutate(Species_ID = sprintf("Species_%06d", row_number()))
tax$taxonomy_full <- apply(
  as.data.frame(tax[, rank_cols, drop = FALSE]),
  1,
  function(z) paste(z[!is.na(z) & nzchar(as.character(z))], collapse = ";")
)
tax$Clean_Name <- coalesce_col(
  tax,
  c("Species", "species", "S", "Taxon")
)
write_csv(tax, file.path(TABLE_DIR, "201_mapa_taxonomia_especies.csv"))

all_metrics <- raw_metrics %>%
  full_join(resid_metrics, by = "Original_Taxon_ID") %>%
  full_join(interaction_metrics, by = "Original_Taxon_ID") %>%
  left_join(tax, by = "Original_Taxon_ID") %>%
  mutate(
    prevalence_n = prev[Original_Taxon_ID],
    total_count = total[Original_Taxon_ID],
    observed_joint_raw =
      raw_Sarle_BC > BC_CUTOFF & raw_Dip_Q_BH <= DIP_FDR,
    observed_joint_residual =
      resid_Sarle_BC > BC_CUTOFF & resid_Dip_Q_BH <= DIP_FDR,
    observed_joint_residual_interaction =
      resid_interaction_Sarle_BC > BC_CUTOFF &
      resid_interaction_Dip_Q_BH <= DIP_FDR
  ) %>%
  arrange(resid_Dip_Q_BH, desc(resid_Sarle_BC))

write_csv(
  all_metrics,
  file.path(TABLE_DIR, "201_metricas_bimodalidad_todas_especies.csv")
)

candidates <- all_metrics %>%
  filter(observed_joint_residual) %>%
  pull(Original_Taxon_ID)

cat("Taxa retained for screen:", ncol(clr), "\n")
cat("Observed additive-residual candidates:", length(candidates), "\n")
cat(
  "Observed interaction-residual sensitivity candidates:",
  sum(all_metrics$observed_joint_residual_interaction, na.rm = TRUE),
  "\n"
)

# ------------------------------------------------------------------------------
# Bootstrap de perfiles completos dentro de localidad
# ------------------------------------------------------------------------------

if (length(candidates)) {
  boot_hits_nominal <- matrix(
    0L,
    nrow = N_BOOT,
    ncol = length(candidates),
    dimnames = list(NULL, candidates)
  )
  boot_hits_bh <- boot_hits_nominal
  boot_bc <- matrix(
    NA_real_,
    nrow = N_BOOT,
    ncol = length(candidates),
    dimnames = list(NULL, candidates)
  )
  boot_dip_p <- boot_bc
  boot_dip_q <- boot_bc

  for (b in seq_len(N_BOOT)) {
    idx <- bootstrap_profile_indices(meta)
    mb <- meta[idx, , drop = FALSE]
    rb <- residualize(
      clr[idx, candidates, drop = FALSE],
      mb,
      interaction = FALSE
    )

    bc_values <- vapply(
      seq_along(candidates),
      function(j) safe_bc(rb[, j]),
      numeric(1)
    )
    dip_p_values <- vapply(
      seq_along(candidates),
      function(j) safe_dip(rb[, j])[["p"]],
      numeric(1)
    )
    dip_q_values <- p.adjust(dip_p_values, method = "BH")

    boot_bc[b, ] <- bc_values
    boot_dip_p[b, ] <- dip_p_values
    boot_dip_q[b, ] <- dip_q_values
    boot_hits_nominal[b, ] <-
      is.finite(bc_values) &
      is.finite(dip_p_values) &
      bc_values > BC_CUTOFF &
      dip_p_values <= DIP_FDR
    boot_hits_bh[b, ] <-
      is.finite(bc_values) &
      is.finite(dip_q_values) &
      bc_values > BC_CUTOFF &
      dip_q_values <= DIP_FDR

    if (b %% 25L == 0L || b == N_BOOT) {
      cat("Bootstrap", b, "of", N_BOOT, "\n")
    }
  }

  hits_bh <- colSums(boot_hits_bh)
  support_bh <- hits_bh / N_BOOT
  hits_nominal <- colSums(boot_hits_nominal)
  support_nominal <- hits_nominal / N_BOOT
  support_ci <- t(
    vapply(
      hits_bh,
      function(k) wilson_interval(k, N_BOOT),
      numeric(2)
    )
  )

  stability <- tibble(
    Original_Taxon_ID = candidates,
    bootstrap_joint_hits_BH = as.integer(hits_bh),
    profile_bootstrap_joint_support = as.numeric(support_bh),
    Stability_Score = as.numeric(support_bh),
    bootstrap_support_MC_SE = sqrt(support_bh * (1 - support_bh) / N_BOOT),
    bootstrap_support_Wilson95_lower = support_ci[, "lower"],
    bootstrap_support_Wilson95_upper = support_ci[, "upper"],
    bootstrap_joint_hits_nominal = as.integer(hits_nominal),
    profile_bootstrap_joint_support_nominal = as.numeric(support_nominal),
    bootstrap_median_BC = apply(boot_bc, 2, median, na.rm = TRUE),
    bootstrap_q025_BC = apply(
      boot_bc,
      2,
      quantile,
      probs = 0.025,
      na.rm = TRUE
    ),
    bootstrap_q975_BC = apply(
      boot_bc,
      2,
      quantile,
      probs = 0.975,
      na.rm = TRUE
    ),
    bootstrap_median_Dip_P = apply(boot_dip_p, 2, median, na.rm = TRUE),
    bootstrap_median_Dip_Q_BH = apply(boot_dip_q, 2, median, na.rm = TRUE),
    bootstrap_iterations = N_BOOT,
    bootstrap_unit = "complete_profile_within_locality",
    bootstrap_FDR_family = "observed_additive_residual_candidates_within_iteration",
    note = paste(
      "Primary support requires BC above cutoff and Dip Q-BH at or below FDR",
      "within each bootstrap iteration; nominal Dip-P support is diagnostic only"
    )
  )

  saveRDS(
    list(
      candidates = candidates,
      hits_BH = boot_hits_bh,
      hits_nominal = boot_hits_nominal,
      Sarle_BC = boot_bc,
      Dip_P = boot_dip_p,
      Dip_Q_BH = boot_dip_q,
      seed = SEED,
      RNGkind = RNGkind()
    ),
    file.path(RDS_DIR, "201_bootstrap_diagnostics.rds")
  )
} else {
  stability <- tibble(
    Original_Taxon_ID = character(),
    bootstrap_joint_hits_BH = integer(),
    profile_bootstrap_joint_support = numeric(),
    Stability_Score = numeric(),
    bootstrap_support_MC_SE = numeric(),
    bootstrap_support_Wilson95_lower = numeric(),
    bootstrap_support_Wilson95_upper = numeric(),
    bootstrap_joint_hits_nominal = integer(),
    profile_bootstrap_joint_support_nominal = numeric(),
    bootstrap_median_BC = numeric(),
    bootstrap_q025_BC = numeric(),
    bootstrap_q975_BC = numeric(),
    bootstrap_median_Dip_P = numeric(),
    bootstrap_median_Dip_Q_BH = numeric(),
    bootstrap_iterations = integer(),
    bootstrap_unit = character(),
    bootstrap_FDR_family = character(),
    note = character()
  )
  saveRDS(
    list(
      candidates = character(),
      seed = SEED,
      RNGkind = RNGkind(),
      note = "No observed additive-residual candidates"
    ),
    file.path(RDS_DIR, "201_bootstrap_diagnostics.rds")
  )
}

write_csv(
  stability,
  file.path(TABLE_DIR, "201_estabilidad_bootstrap_por_perfil.csv")
)

# ------------------------------------------------------------------------------
# Seleccion conservadora y sensibilidad de residualizacion
# ------------------------------------------------------------------------------

final <- all_metrics %>%
  filter(observed_joint_residual) %>%
  left_join(stability, by = "Original_Taxon_ID") %>%
  mutate(
    conservative_bimodal =
      profile_bootstrap_joint_support >= STABILITY_CUTOFF,
    interaction_residual_sensitivity_pass =
      observed_joint_residual_interaction,
    Bimodality_Class = if_else(
      interaction_residual_sensitivity_pass,
      "conservative_primary_plus_interaction_sensitivity",
      "conservative_primary_only"
    ),
    Stability_Score_definition =
      "profile bootstrap support using Sarle BC plus Dip Q-BH",
    selection_rule = paste0(
      "additive_residual_BC>", BC_CUTOFF,
      "; additive_residual_Dip_Q_BH<=", DIP_FDR,
      "; profile_bootstrap_BC_and_Dip_Q_BH_support>=",
      STABILITY_CUTOFF
    ),
    sensitivity_rule = paste0(
      "locality_x_depth_residual_BC>", BC_CUTOFF,
      "; locality_x_depth_residual_Dip_Q_BH<=", DIP_FDR
    )
  ) %>%
  filter(conservative_bimodal) %>%
  arrange(desc(profile_bootstrap_joint_support), resid_Dip_Q_BH)

write_csv(
  final,
  file.path(TABLE_DIR, "201_especies_bimodales_conservadoras.csv")
)

sensitivity_summary <- tibble(
  criterion = c(
    "observed_raw_BC_plus_Dip_Q_BH",
    "observed_additive_residual_BC_plus_Dip_Q_BH",
    "observed_interaction_residual_BC_plus_Dip_Q_BH",
    "final_bootstrap_BH_support",
    "final_and_interaction_residual_sensitivity"
  ),
  n_taxa = c(
    sum(all_metrics$observed_joint_raw, na.rm = TRUE),
    sum(all_metrics$observed_joint_residual, na.rm = TRUE),
    sum(all_metrics$observed_joint_residual_interaction, na.rm = TRUE),
    nrow(final),
    sum(final$interaction_residual_sensitivity_pass, na.rm = TRUE)
  ),
  role = c(
    "diagnostic",
    "primary_observed_discovery",
    "sensitivity_only",
    "primary_final_selection",
    "stronger_cross_residualization_subset"
  )
)
write_csv(
  sensitivity_summary,
  file.path(TABLE_DIR, "201_resumen_sensibilidad_residualizacion.csv")
)

# ------------------------------------------------------------------------------
# Compuerta final
# ------------------------------------------------------------------------------

probability_columns <- c(
  "raw_Dip_P",
  "raw_Dip_Q_BH",
  "resid_Dip_P",
  "resid_Dip_Q_BH",
  "resid_interaction_Dip_P",
  "resid_interaction_Dip_Q_BH"
)
probability_values <- unlist(
  all_metrics[probability_columns],
  use.names = FALSE
)
probability_values <- probability_values[is.finite(probability_values)]

expected_final_ids <- stability %>%
  filter(profile_bootstrap_joint_support >= STABILITY_CUTOFF) %>%
  pull(Original_Taxon_ID)

final_gate <- tibble(
  gate = c(
    "all_screened_taxa_map_to_taxonomy",
    "reported_probabilities_in_unit_interval",
    "bootstrap_completed_for_all_candidates",
    "bootstrap_support_in_unit_interval",
    "final_ids_unique",
    "final_selection_matches_prespecified_rule",
    "final_support_meets_cutoff"
  ),
  pass = c(
    all(all_metrics$Original_Taxon_ID %in% tax$Original_Taxon_ID),
    length(probability_values) > 0L &&
      all(probability_values >= 0 & probability_values <= 1),
    nrow(stability) == length(candidates) &&
      (nrow(stability) == 0L ||
         all(stability$bootstrap_iterations == N_BOOT)),
    nrow(stability) == 0L ||
      all(
        stability$profile_bootstrap_joint_support >= 0 &
          stability$profile_bootstrap_joint_support <= 1 &
          stability$profile_bootstrap_joint_support_nominal >= 0 &
          stability$profile_bootstrap_joint_support_nominal <= 1
      ),
    !anyDuplicated(final$Original_Taxon_ID),
    setequal(final$Original_Taxon_ID, expected_final_ids),
    nrow(final) == 0L ||
      all(final$profile_bootstrap_joint_support >= STABILITY_CUTOFF)
  ),
  detail = c(
    paste0(
      "mapped=",
      sum(all_metrics$Original_Taxon_ID %in% tax$Original_Taxon_ID),
      "/",
      nrow(all_metrics)
    ),
    "finite Dip P/Q values constrained to [0,1]",
    paste0("candidates=", length(candidates), "; bootstrapped=", nrow(stability)),
    "BH-adjusted and nominal bootstrap supports constrained to [0,1]",
    paste0("n_final=", nrow(final)),
    "observed additive residual criteria plus BH-adjusted bootstrap support",
    paste0("cutoff=", STABILITY_CUTOFF)
  )
)
write_csv(final_gate, file.path(TABLE_DIR, "201_compuerta_seleccion_bimodal.csv"))
if (!all(final_gate$pass)) {
  stop("Fallo la compuerta final de seleccion bimodal", call. = FALSE)
}

# ------------------------------------------------------------------------------
# Figuras
# ------------------------------------------------------------------------------

plot_df <- all_metrics %>%
  mutate(selected = Original_Taxon_ID %in% final$Original_Taxon_ID)

p_screen <- ggplot(
  plot_df,
  aes(
    resid_Sarle_BC,
    -log10(pmax(resid_Dip_Q_BH, .Machine$double.xmin)),
    color = selected
  )
) +
  geom_point(alpha = 0.65, size = 1.4) +
  geom_vline(xintercept = BC_CUTOFF, linetype = 2) +
  geom_hline(yintercept = -log10(DIP_FDR), linetype = 2) +
  scale_color_manual(values = c(`FALSE` = "grey65", `TRUE` = "#0072B2")) +
  labs(
    x = "Sarle BC after locality + depth residualization",
    y = expression(-log[10]("Dip q (BH)")),
    color = "Retained"
  ) +
  theme_bw(base_size = 10)

ggsave(
  file.path(PLOT_DIR, "201_screen_bimodalidad_residualizada.png"),
  p_screen,
  width = 7,
  height = 5,
  dpi = 300
)
ggsave(
  file.path(PLOT_DIR, "201_screen_bimodalidad_residualizada.pdf"),
  p_screen,
  width = 7,
  height = 5
)

if (nrow(stability) > 0L) {
  support_plot_df <- stability %>%
    left_join(
      all_metrics %>%
        select(
          Original_Taxon_ID,
          observed_joint_residual_interaction
        ),
      by = "Original_Taxon_ID"
    )

  p_support <- ggplot(
    support_plot_df,
    aes(
      profile_bootstrap_joint_support_nominal,
      profile_bootstrap_joint_support,
      color = observed_joint_residual_interaction
    )
  ) +
    geom_abline(slope = 1, intercept = 0, color = "grey70", linetype = 2) +
    geom_hline(yintercept = STABILITY_CUTOFF, linetype = 3) +
    geom_point(alpha = 0.8, size = 1.8) +
    scale_color_manual(
      values = c(`FALSE` = "grey60", `TRUE` = "#D55E00"),
      name = "Passes locality x depth\nsensitivity"
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      x = "Bootstrap support using nominal Dip P",
      y = "Bootstrap support using Dip Q (BH)"
    ) +
    theme_bw(base_size = 10)

  ggsave(
    file.path(PLOT_DIR, "201_soporte_bootstrap_nominal_vs_BH.png"),
    p_support,
    width = 6.4,
    height = 5.5,
    dpi = 300
  )
  ggsave(
    file.path(PLOT_DIR, "201_soporte_bootstrap_nominal_vs_BH.pdf"),
    p_support,
    width = 6.4,
    height = 5.5
  )
}

# ------------------------------------------------------------------------------
# Procedencia, resumen y cierre transaccional de LATEST
# ------------------------------------------------------------------------------

input_md5 <- tibble(
  input = c("ps_rds", "meta_csv"),
  path = c(PS_RDS, META_CSV),
  md5 = unname(as.character(tools::md5sum(c(PS_RDS, META_CSV))))
)
write_csv(input_md5, file.path(TABLE_DIR, "201_input_md5.csv"))

run_info <- tibble(
  parameter = c(
    "script",
    "run001",
    "ps_rds",
    "meta_csv",
    "profile_col",
    "locality_col",
    "depth_col",
    "seed",
    "rng_kind",
    "n_boot",
    "pseudocount",
    "min_prev_n",
    "min_total",
    "min_sd_clr",
    "bc_cutoff",
    "bc_definition",
    "dip_fdr",
    "dip_interpretation",
    "stability_cutoff",
    "primary_residualization",
    "sensitivity_residualization",
    "bootstrap_FDR_family",
    "n_samples",
    "n_profiles",
    "n_taxa_before_filter",
    "n_taxa_screened",
    "n_observed_raw_candidates",
    "n_observed_residual_candidates",
    "n_observed_interaction_sensitivity_candidates",
    "n_final_conservative",
    "n_final_passing_interaction_sensitivity",
    "latest_attempt_file",
    "latest_file"
  ),
  value = as.character(c(
    SCRIPT_ID,
    RUN001,
    PS_RDS,
    META_CSV,
    PROFILE_COL,
    LOCALITY_COL,
    DEPTH_COL,
    SEED,
    paste(RNGkind(), collapse = ";"),
    N_BOOT,
    PSEUDOCOUNT,
    MIN_PREV_N,
    MIN_TOTAL,
    MIN_SD_CLR,
    BC_CUTOFF,
    "(sample_skewness^2 + 1) / sample_raw_kurtosis",
    DIP_FDR,
    "evidence_against_unimodality_not_proof_of_exactly_two_modes",
    STABILITY_CUTOFF,
    "CLR ~ locality + factor(depth_cm)",
    "CLR ~ locality * factor(depth_cm)",
    "observed_additive_residual_candidates_within_iteration",
    nrow(meta),
    n_distinct(meta$profile_id),
    nrow(otu),
    ncol(clr),
    sum(all_metrics$observed_joint_raw, na.rm = TRUE),
    length(candidates),
    sum(all_metrics$observed_joint_residual_interaction, na.rm = TRUE),
    nrow(final),
    sum(final$interaction_residual_sensitivity_pass, na.rm = TRUE),
    LATEST_ATTEMPT_FILE,
    LATEST_FILE
  ))
)
write_csv(run_info, file.path(TABLE_DIR, "201_parametros_y_resumen.csv"))

report <- c(
  "201 conservative bimodality discovery",
  "======================================",
  "",
  paste0("Phyloseq: ", PS_RDS),
  paste0("Metadata: ", META_CSV),
  paste0("Output: ", OUT_DIR),
  "",
  "Primary selection:",
  "  - prevalence, total abundance and CLR-SD filters are outcome-independent",
  "  - restoration stage is never used",
  "  - CLR is residualized by locality + factor(depth_cm)",
  "  - observed Sarle BC and Dip Q-BH define candidates",
  "  - complete profiles are bootstrapped within locality",
  "  - each iteration adjusts Dip P by BH across observed candidates",
  "  - final retention uses BH-adjusted bootstrap joint support",
  "",
  "Sensitivity:",
  "  - observed CLR is additionally residualized by locality * factor(depth_cm)",
  "  - this tests whether locality-specific vertical gradients explain candidates",
  "  - sensitivity status is reported but does not redefine primary selection",
  "",
  "Interpretation limit:",
  "  Dip rejects unimodality; combined Dip and Sarle evidence is compatible with",
  "  bimodality or multimodality but does not establish alternative stable states,",
  "  hysteresis, causality or temporal switching.",
  "",
  paste0("Taxa screened: ", ncol(clr)),
  paste0("Observed additive-residual candidates: ", length(candidates)),
  paste0("Final conservative taxa: ", nrow(final)),
  paste0(
    "Final taxa passing locality x depth sensitivity: ",
    sum(final$interaction_residual_sensitivity_pass, na.rm = TRUE)
  )
)
writeLines(report, file.path(OUT_DIR, "201_report.txt"))

writeLines(
  capture.output(sessionInfo()),
  file.path(LOG_DIR, "201_sessionInfo.txt")
)

# LATEST solo se actualiza despues de pasar todas las compuertas y escribir
# todos los productos principales.
writeLines(OUT_DIR, LATEST_FILE)

cat("===== DONE 201 =====\n")
cat("Final conservative taxa:", nrow(final), "\n")
cat(
  "Passing locality x depth sensitivity:",
  sum(final$interaction_residual_sensitivity_pass, na.rm = TRUE),
  "\n"
)
cat("Output:", OUT_DIR, "\n")
cat("Latest:", LATEST_FILE, "\n")
cat("Finished UTC:", format(Sys.time(), tz = "UTC"), "\n")
