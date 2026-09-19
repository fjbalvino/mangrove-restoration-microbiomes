#!/usr/bin/env Rscript
# ============================================================
# 05_01_composition_reference.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# 7,126 KO and 6,681 PFAM; composition, dispersion and matched-reference distances. 998 retained composition permutations.
# Inputs (source expressions; complete list in docs/contracts/05_01_composition_reference.json):
#   trimws(readLines(hit[[1L]], n = 1L, warn = FALSE))
#   run202 <- trimws(readLines(pointer202, n = 1L, warn = FALSE))
#   x <- read.delim(
#   md <- read.csv(
# Outputs (source expressions; complete list in contract):
#   writeLines(
#   write.table(
#   ggplot2::ggsave(
#   write_tsv(cell_table, file.path(DIRS["tables"], "602_locality_stage_design_cells.tsv"))
#   write_tsv(
#   saveRDS(
# Algorithmic provenance:
# Technical KO/PFAM filters; CLR/Aitchison; profile-blocked composition and reference-distance tests.
#   McArdle & Anderson (2001), doi:10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2.
# Source SHA-256: bb4462cce0636527c7a39316738c1fcb272eb1412d91cec18e9a2de14c1e0218
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


# 602: independent global analysis of functional composition.
# KO is primary; PFAM is cross-ontology confirmation.
# Complete three-depth profiles are permuted within locality.

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (!length(i)) return(default)
  if (i == length(args)) stop("Argument without value: ", flag)
  args[[i + 1L]]
}
read_latest_any <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) {
    stop("No LATEST file for 504. Run 504 or provide --run504.")
  }
  trimws(readLines(hit[[1L]], n = 1L, warn = FALSE))
}

# INPUTS_082A_CANON202
BASE <- get_arg("--base", "/home/fjbalvino/Tipping_points")
RESULTS <- get_arg("--out_root", file.path(BASE, "resultados_finales"))
RUN082 <- get_arg(
  "--run082",
  file.path(
    BASE, "Results",
    "082A_build_complete_functional_matrices_KO_PFAM_20260726_011523"
  )
)

META_FILE <- get_arg("--meta", NULL)
if (is.null(META_FILE)) {
  pointer202 <- file.path(
    BASE, "resultados_finales",
    "LATEST_202_construir_matriz_CLR_y_eje_bimodal.txt"
  )
  if (!file.exists(pointer202)) stop("Missing pointer: ", pointer202)
  run202 <- trimws(readLines(pointer202, n = 1L, warn = FALSE))
  META_FILE <- file.path(
    run202, "tables", "202_metadata_raw1031_alineada.csv"
  )
}

KO_FILE <- file.path(
  RUN082, "tables", "082A_KEGG_ko_abundance_samples_x_functions.tsv.gz"
)
PFAM_FILE <- file.path(
  RUN082, "tables", "082A_PFAM_abundance_samples_x_functions.tsv.gz"
)
input_files <- c(metadata = META_FILE, KO = KO_FILE, PFAM = PFAM_FILE)
if (!all(file.exists(input_files))) {
  stop("Missing inputs: ", paste(input_files[!file.exists(input_files)], collapse = ", "))
}

N_PERM <- as.integer(get_arg("--n_perm", "999"))
SEED <- as.integer(get_arg("--seed", "123"))
MIN_PREV_PRIMARY <- 26L
MIN_PREV_PERMISSIVE <- 5L
MIN_TOTAL <- 10
PC_PRIMARY <- 1
PC_SENSITIVITY <- 0.5
DEPTHS <- c(5, 20, 40)

STAGES <- c(
  "Degraded",
  "Early restoration",
  "Intermediate restoration",
  "Advanced restoration",
  "Conserved"
)

needed <- c("vegan", "permute", "ggplot2", "dplyr")
missing <- needed[
  !vapply(needed, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "))
}

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S", tz = "UTC")

OUT <- file.path(
  RESULTS,
  paste0(
    "602_validar_composicion_funcional_splitplot_",
    stamp
  )
)

DIRS <- file.path(
  OUT,
  c("tables", "figures", "rds", "logs", "provenance")
)
names(DIRS) <- c(
  "tables", "figures", "rds", "logs", "provenance"
)

invisible(
  lapply(
    c(OUT, DIRS),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)

writeLines(
  OUT,
  file.path(
    RESULTS,
    "LATEST_ATTEMPT_602_validar_composicion_funcional_splitplot.txt"
  )
)

log_con <- file(
  file.path(DIRS["logs"], "602_log.txt"),
  "wt"
)

sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat("602 functional global validation\n")
cat("Started UTC:", format(Sys.time(), tz = "UTC"), "\n")
cat("Output:", OUT, "\n\n")

write_tsv <- function(x, path) {
  write.table(
    x,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

read_matrix <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " input not found: ", path)
  }

  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)

  x <- read.delim(
    con,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character",
    quote = "",
    comment.char = ""
  )

  if (ncol(x) < 2L) {
    stop(label, " matrix is empty")
  }

  ids <- trimws(as.character(x[[1]]))
  features <- names(x)[-1]

  if (any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(label, " contains invalid sample IDs")
  }

  if (any(!nzchar(features)) || anyDuplicated(features)) {
    stop(label, " contains invalid feature IDs")
  }

  m <- matrix(
    suppressWarnings(as.numeric(as.matrix(x[-1L]))),
    nrow = nrow(x),
    ncol = ncol(x) - 1L
  )
  rownames(m) <- ids
  colnames(m) <- features

  if (any(!is.finite(m)) || any(m < 0)) {
    stop(label, " contains invalid abundances")
  }

  if (any(rowSums(m) <= 0)) {
    stop(label, " contains an empty sample")
  }

  m
}

feature_qc <- function(m) {
  data.frame(
    feature = colnames(m),
    prevalence_n = colSums(m > 0),
    total_abundance = colSums(m),
    variance_log1p = apply(log1p(m), 2, stats::var),
    stringsAsFactors = FALSE
  )
}

filter_matrix <- function(m, qc, min_prev) {
  keep <- (
    qc$prevalence_n >= min_prev &
    qc$total_abundance >= MIN_TOTAL &
    is.finite(qc$variance_log1p) &
    qc$variance_log1p > 0
  )

  if (!any(keep)) {
    stop(
      "No features passed filter prevalence >= ",
      min_prev
    )
  }

  out <- m[, qc$feature[keep], drop = FALSE]
  if (any(rowSums(out) <= 0)) {
    stop("Empty samples after filtering at prevalence >= ", min_prev)
  }
  out
}

clr <- function(m, pc) {
  z <- log(m + pc)
  sweep(z, 1, rowMeans(z), "-")
}

make_distance <- function(m, method, pc) {
  if (method == "Aitchison") {
    return(stats::dist(clr(m, pc)))
  }

  if (method == "Bray-Hellinger") {
    h <- sqrt(m / rowSums(m))
    return(vegan::vegdist(h, method = "bray"))
  }

  stop("Unknown distance: ", method)
}

fit_permanova <- function(
  d,
  md,
  permutations,
  stage_model
) {
  if (stage_model == "factor") {
    fit <- vegan::adonis2(
      d ~ locality + depth_f + restoration4,
      data = md,
      permutations = permutations,
      by = "margin"
    )

    term <- "restoration4"
  } else {
    fit <- vegan::adonis2(
      d ~ locality + depth_f + restoration_ord,
      data = md,
      permutations = permutations,
      by = "margin"
    )

    term <- "restoration_ord"
  }

  z <- fit[term, , drop = FALSE]

  data.frame(
    stage_model = stage_model,
    term = term,
    df = unname(z$Df),
    sum_of_squares = unname(z$SumOfSqs),
    r2 = unname(z$R2),
    pseudo_f = unname(z$F),
    p_value = unname(z$`Pr(>F)`),
    n_permutations = nrow(permutations),
    minimum_attainable_p = 1 / (nrow(permutations) + 1),
    seed = SEED,
    permutation_scheme = "whole_profiles_within_locality_depth_order_preserved"
  )
}

fit_dispersion <- function(
  d,
  md,
  permutations
) {
  bd <- vegan::betadisper(
    d,
    group = md$restoration4,
    type = "centroid",
    bias.adjust = TRUE
  )

  pt <- vegan::permutest(
    bd,
    permutations = permutations
  )

  tb <- as.data.frame(
    pt$tab,
    check.names = FALSE
  )

  p_col <- grep("^Pr", names(tb), value = TRUE)[1]
  sum_col <- grep("^Sum", names(tb), value = TRUE)[1]
  mean_col <- grep("^Mean", names(tb), value = TRUE)[1]

  data.frame(
    df = tb$Df[1],
    sum_of_squares = tb[[sum_col]][1],
    mean_square = tb[[mean_col]][1],
    pseudo_f = tb$F[1],
    p_value = tb[[p_col]][1],
    n_permutations = nrow(permutations),
    minimum_attainable_p = 1 / (nrow(permutations) + 1),
    seed = SEED,
    permutation_scheme = "whole_profiles_within_locality_depth_order_preserved"
  )
}

pcoa_scores <- function(
  d,
  md,
  ontology
) {
  fit <- stats::cmdscale(
    d,
    k = 2,
    eig = TRUE
  )

  positive <- fit$eig[fit$eig > 0]
  explained <- 100 * fit$eig[1:2] / sum(positive)

  out <- data.frame(
    sample_id = rownames(fit$points),
    axis1 = fit$points[, 1],
    axis2 = fit$points[, 2],
    axis1_percent = explained[1],
    axis2_percent = explained[2],
    ontology = ontology
  )

  cbind(
    out,
    md[
      match(out$sample_id, md$sample_id),
      c(
        "restoration4",
        "restoration_ord",
        "locality",
        "depth_cm",
        "depth_f",
        "profile_id"
      )
    ]
  )
}

plot_pcoa <- function(x, ontology) {
  colors <- c(
    "Degraded" = "#B84A3C",
    "Early restoration" = "#E69F00",
    "Intermediate restoration" = "#8C6BB1",
    "Advanced restoration" = "#2A9D8F",
    "Conserved" = "#276FBF"
  )

  p1 <- round(unique(x$axis1_percent)[1], 1)
  p2 <- round(unique(x$axis2_percent)[1], 1)

  p <- ggplot2::ggplot(
    x,
    ggplot2::aes(
      axis1,
      axis2,
      color = restoration4,
      shape = depth_f
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "grey85",
      linewidth = 0.25
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      color = "grey85",
      linewidth = 0.25
    ) +
    ggplot2::geom_point(
      size = 3,
      alpha = 0.88
    ) +
    ggplot2::facet_wrap(
      ~ locality,
      nrow = 1
    ) +
    ggplot2::scale_color_manual(
      values = colors,
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "5" = 16,
        "20" = 17,
        "40" = 15
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      x = paste0("PCoA1 (", p1, "%)"),
      y = paste0("PCoA2 (", p2, "%)"),
      color = "Restoration stage",
      shape = "Depth (cm)"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(
        fill = "grey95"
      ),
      legend.position = "bottom"
    )

  stem <- file.path(
    DIRS["figures"],
    paste0(
      "602_",
      ontology,
      "_primary_Aitchison_PCoA"
    )
  )

  ggplot2::ggsave(
    paste0(stem, ".png"),
    p,
    width = 12,
    height = 5.8,
    dpi = 400,
    bg = "white"
  )

  ggplot2::ggsave(
    paste0(stem, ".pdf"),
    p,
    width = 12,
    height = 5.8,
    bg = "white"
  )
}

cat("[1/7] Reading inputs\n")

md <- read.csv(
  META_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required <- c(
  "sample_id",
  "restoration4",
  "depth_cm",
  "locality",
  "lat_block"
)

missing_columns <- setdiff(required, names(md))

if (length(missing_columns)) {
  stop(
    "Missing metadata columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

md$sample_id <- trimws(as.character(md$sample_id))
md$restoration4 <- trimws(
  as.character(md$restoration4)
)
md$locality <- trimws(as.character(md$locality))
md$lat_block_chr <- trimws(
  as.character(md$lat_block)
)
md$depth_cm <- suppressWarnings(
  as.numeric(md$depth_cm)
)

if (
  anyDuplicated(md$sample_id) ||
  any(!nzchar(md$sample_id))
) {
  stop("Metadata sample IDs are invalid")
}

if (any(!md$restoration4 %in% STAGES)) {
  stop("Unexpected restoration stage")
}

if (anyNA(md$depth_cm)) {
  stop("Invalid depth_cm")
}

md$restoration4 <- factor(
  md$restoration4,
  levels = STAGES
)

md$restoration_ord <- (
  match(as.character(md$restoration4), STAGES) - 1L
)

md$locality <- factor(md$locality)

md$depth_f <- factor(
  md$depth_cm,
  levels = DEPTHS
)

md$profile_id <- interaction(
  md$locality,
  md$lat_block_chr,
  drop = TRUE,
  sep = "::"
)

raw <- list(
  KEGG_ko = read_matrix(
    KO_FILE,
    "KEGG_ko"
  ),
  PFAM = read_matrix(
    PFAM_FILE,
    "PFAM"
  )
)

for (nm in names(raw)) {
  if (!setequal(md$sample_id, rownames(raw[[nm]]))) {
    stop(
      nm,
      " sample IDs do not match metadata"
    )
  }
}

cat("[2/7] Auditing repeated profiles\n")

md <- md[
  order(
    md$locality,
    md$profile_id,
    md$depth_cm,
    md$sample_id
  ),
]

rownames(md) <- md$sample_id

raw <- lapply(
  raw,
  function(m) {
    m[md$sample_id, , drop = FALSE]
  }
)

idx <- split(
  seq_len(nrow(md)),
  md$profile_id
)

profile_ok <- vapply(
  idx,
  function(i) {
    length(i) == 3L &&
      identical(
        sort(md$depth_cm[i]),
        DEPTHS
      ) &&
      length(unique(md$restoration4[i])) == 1L &&
      length(unique(md$locality[i])) == 1L
  },
  logical(1)
)

if (
  nrow(md) != 51L ||
  length(idx) != 17L ||
  !all(profile_ok)
) {
  stop(
    paste(
      "The expected 51-sample/17-profile/",
      "three-depth design was not recovered"
    )
  )
}

cell_table <- as.data.frame(table(md$locality, md$restoration4))
names(cell_table) <- c("locality", "restoration4", "n_samples")
cell_table$n_profiles <- mapply(
  function(loc, stage) {
    length(unique(md$profile_id[md$locality == loc & md$restoration4 == stage]))
  },
  cell_table$locality,
  cell_table$restoration4
)
cell_table$empty <- cell_table$n_profiles == 0L
cell_table$single_profile <- cell_table$n_profiles == 1L
write_tsv(cell_table, file.path(DIRS["tables"], "602_locality_stage_design_cells.tsv"))
if (any(cell_table$empty)) {
  warning("Some locality x stage cells are empty; model rank and estimability are checked explicitly.")
}

X_primary <- model.matrix(~ locality + depth_f + restoration4, data = md)
if (qr(X_primary)$rank < ncol(X_primary)) {
  stop("Primary factor model is not estimable: design matrix is rank deficient")
}

cat("[3/7] Applying outcome-independent filters\n")

qc <- lapply(raw, feature_qc)

filtered <- lapply(
  names(raw),
  function(nm) {
    list(
      common = filter_matrix(
        raw[[nm]],
        qc[[nm]],
        MIN_PREV_PRIMARY
      ),
      permissive = filter_matrix(
        raw[[nm]],
        qc[[nm]],
        MIN_PREV_PERMISSIVE
      )
    )
  }
)

names(filtered) <- names(raw)

for (nm in names(qc)) {
  qc[[nm]]$keep_primary <- (
    qc[[nm]]$prevalence_n >= MIN_PREV_PRIMARY &
    qc[[nm]]$total_abundance >= MIN_TOTAL &
    qc[[nm]]$variance_log1p > 0
  )

  qc[[nm]]$keep_permissive <- (
    qc[[nm]]$prevalence_n >= MIN_PREV_PERMISSIVE &
    qc[[nm]]$total_abundance >= MIN_TOTAL &
    qc[[nm]]$variance_log1p > 0
  )

  con <- gzfile(
    file.path(
      DIRS["tables"],
      paste0(
        "602_",
        nm,
        "_feature_QC.tsv.gz"
      )
    ),
    "wt"
  )

  write.table(
    qc[[nm]],
    con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  close(con)
}

filter_summary <- do.call(
  rbind,
  lapply(
    names(raw),
    function(nm) {
      data.frame(
        ontology = nm,
        complete_features = ncol(raw[[nm]]),
        primary_features = ncol(
          filtered[[nm]]$common
        ),
        permissive_features = ncol(
          filtered[[nm]]$permissive
        ),
        primary_prevalence_n = MIN_PREV_PRIMARY,
        permissive_prevalence_n =
          MIN_PREV_PERMISSIVE,
        minimum_total_abundance = MIN_TOTAL,
        selection_uses_HI_or_stage = FALSE
      )
    }
  )
)

write_tsv(
  filter_summary,
  file.path(
    DIRS["tables"],
    "602_feature_filter_summary.tsv"
  )
)

print(filter_summary)

cat("[4/7] Constructing profile-restricted permutations\n")

control <- permute::how(
  nperm = N_PERM,
  blocks = md$locality,
  plots = permute::Plots(
    strata = md$profile_id,
    type = "free"
  ),
  within = permute::Within(type = "none")
)

set.seed(SEED)

perms_raw <- permute::shuffleSet(
  nrow(md),
  nset = N_PERM,
  control = control
)
perms <- as.matrix(unique(as.data.frame(perms_raw)))

audit_one <- function(p) {
  c(
    locality = all(
      as.character(md$locality[p]) ==
        as.character(md$locality)
    ),
    depth = all(
      md$depth_cm[p] == md$depth_cm
    ),
    whole_profile = all(
      vapply(
        idx,
        function(i) {
          length(
            unique(md$profile_id[p[i]])
          ) == 1L
        },
        logical(1)
      )
    )
  )
}

perm_audit <- t(
  apply(perms, 1, audit_one)
)

perm_audit <- data.frame(
  permutation = seq_len(nrow(perms)),
  perm_audit
)

perm_audit$valid <- with(
  perm_audit,
  locality & depth & whole_profile
)

if (!all(perm_audit$valid)) {
  stop("Restricted permutation audit failed")
}

write_tsv(
  perm_audit,
  file.path(
    DIRS["tables"],
    "602_permutation_audit.tsv"
  )
)

write_tsv(
  data.frame(
    requested = N_PERM,
    generated = nrow(perms),
    unique = nrow(
      unique(as.data.frame(perms))
    ),
    valid = sum(perm_audit$valid),
    invalid = sum(!perm_audit$valid),
    scheme = paste(
      "complete profiles within locality;",
      "depth order preserved"
    )
  ),
  file.path(
    DIRS["tables"],
    "602_permutation_summary.tsv"
  )
)

cat("[5/7] Running PERMANOVA and BETADISPER\n")

specs <- data.frame(
  analysis_id = c(
    "primary_common_CLR_pc1",
    "sensitivity_common_CLR_pc0.5",
    "sensitivity_permissive_CLR_pc1",
    "sensitivity_common_Hellinger_Bray"
  ),
  role = c(
    "primary",
    "pseudocount",
    "prevalence",
    "distance"
  ),
  filter = c(
    "common",
    "common",
    "permissive",
    "common"
  ),
  distance = c(
    "Aitchison",
    "Aitchison",
    "Aitchison",
    "Bray-Hellinger"
  ),
  pseudocount = c(
    PC_PRIMARY,
    PC_SENSITIVITY,
    PC_PRIMARY,
    NA_real_
  )
)

permanova <- list()
dispersion <- list()
pcoa <- list()
primary_data <- list(metadata = md)

k <- 1L
j <- 1L

for (ontology in names(filtered)) {
  for (s in seq_len(nrow(specs))) {
    sp <- specs[s, ]

    m <- filtered[[ontology]][[sp$filter]]

    d <- make_distance(
      m,
      sp$distance,
      sp$pseudocount
    )

    # Factor is the only inferential stage model. An ordinal trend is not used
    # because the five cross-sectional stages are not equally spaced in time.
    for (stage_model in "factor") {
      z <- fit_permanova(
        d,
        md,
        perms,
        stage_model
      )

      z$ontology <- ontology
      z$analysis_id <- sp$analysis_id
      z$analysis_role <- sp$role
      z$filter <- sp$filter
      z$n_features <- ncol(m)
      z$distance <- sp$distance
      z$pseudocount <- sp$pseudocount

      permanova[[k]] <- z
      k <- k + 1L
    }

    z <- fit_dispersion(
      d,
      md,
      perms
    )

    z$ontology <- ontology
    z$analysis_id <- sp$analysis_id
    z$analysis_role <- sp$role
    z$filter <- sp$filter
    z$n_features <- ncol(m)
    z$distance <- sp$distance
    z$pseudocount <- sp$pseudocount

    dispersion[[j]] <- z
    j <- j + 1L

    if (sp$role == "primary") {
      pcoa[[ontology]] <- pcoa_scores(
        d,
        md,
        ontology
      )

      primary_data[[ontology]] <- list(
        abundance = m,
        clr = clr(m, sp$pseudocount),
        distance = d,
        features = colnames(m),
        pseudocount = sp$pseudocount,
        min_prevalence_n = MIN_PREV_PRIMARY
      )
    }
  }
}

permanova <- do.call(rbind, permanova)
dispersion <- do.call(rbind, dispersion)
pcoa <- do.call(rbind, pcoa)

front <- c(
  "ontology",
  "analysis_id",
  "analysis_role",
  "filter",
  "n_features",
  "distance",
  "pseudocount"
)

permanova <- permanova[
  ,
  c(
    front,
    setdiff(names(permanova), front)
  )
]

dispersion <- dispersion[
  ,
  c(
    front,
    setdiff(names(dispersion), front)
  )
]

write_tsv(
  permanova,
  file.path(
    DIRS["tables"],
    "602_PERMANOVA_results.tsv"
  )
)

write_tsv(
  dispersion,
  file.path(
    DIRS["tables"],
    "602_BETADISPER_results.tsv"
  )
)

write_tsv(
  pcoa,
  file.path(
    DIRS["tables"],
    "602_primary_PCoA_scores.tsv"
  )
)

primary <- permanova[
  permanova$analysis_role == "primary" &
    permanova$stage_model == "factor",
]

disp_primary <- dispersion[
  dispersion$analysis_role == "primary",
  c("ontology", "p_value")
]

names(disp_primary)[2] <- "dispersion_p_value"

primary <- merge(
  primary,
  disp_primary,
  by = "ontology",
  sort = FALSE
)

primary$interpretation_code <- ifelse(
  primary$p_value < 0.05 &
    primary$dispersion_p_value >= 0.05,
  paste(
    "stage_location_signal_without_detected",
    "dispersion_difference",
    sep = "_"
  ),
  ifelse(
    primary$p_value < 0.05,
    "stage_signal_with_dispersion_caution",
    "no_global_stage_signal_detected"
  )
)

primary$scope <- paste(
  "internal functional gene-potential",
  "composition"
)

write_tsv(
  primary,
  file.path(
    DIRS["tables"],
    "602_primary_summary.tsv"
  )
)

print(
  primary[
    ,
    c(
      "ontology",
      "n_features",
      "r2",
      "pseudo_f",
      "p_value",
      "dispersion_p_value",
      "interpretation_code"
    )
  ]
)

# Independent direction analysis: distance to the Conserved reference of the
# same locality and depth. Conserved anchors are excluded from inference when
# leave-one-out is impossible. No HI/HFR or stage-selected functions are used.
partial_f_stage <- function(y, dat) {
  reduced <- stats::lm(y ~ locality, data = dat)
  full <- stats::lm(y ~ locality + restoration4, data = dat)
  rss0 <- sum(stats::resid(reduced)^2)
  rss1 <- sum(stats::resid(full)^2)
  df_num <- stats::df.residual(reduced) - stats::df.residual(full)
  df_den <- stats::df.residual(full)
  c(
    f = ((rss0 - rss1) / df_num) / (rss1 / df_den),
    df_num = df_num,
    df_den = df_den
  )
}

reference_samples <- list()
reference_profiles <- list()
reference_tests <- list()

for (ontology in names(primary_data)[names(primary_data) != "metadata"]) {
  z <- primary_data[[ontology]]$clr
  z <- z[md$sample_id, , drop = FALSE]
  dist_to_reference <- rep(NA_real_, nrow(md))

  for (i in seq_len(nrow(md))) {
    ref <- which(
      md$locality == md$locality[i] &
        md$depth_cm == md$depth_cm[i] &
        md$restoration4 == "Conserved"
    )
    if (md$restoration4[i] == "Conserved") {
      ref <- setdiff(ref, i)
    }
    if (length(ref)) {
      centroid <- colMeans(z[ref, , drop = FALSE])
      dist_to_reference[i] <- sqrt(sum((z[i, ] - centroid)^2))
    }
  }

  sample_out <- data.frame(
    sample_id = md$sample_id,
    profile_id = md$profile_id,
    locality = md$locality,
    depth_cm = md$depth_cm,
    restoration4 = md$restoration4,
    ontology = ontology,
    distance_to_local_depth_Conserved = dist_to_reference,
    used_for_inference = md$restoration4 != "Conserved" & is.finite(dist_to_reference),
    reference_rule = "same_locality_and_depth; Conserved self excluded"
  )
  reference_samples[[ontology]] <- sample_out

  missing_anchor <- with(
    md,
    expand.grid(locality = levels(locality), depth_cm = DEPTHS)
  )
  missing_anchor$n_conserved <- mapply(
    function(loc, dep) {
      sum(md$locality == loc & md$depth_cm == dep & md$restoration4 == "Conserved")
    },
    missing_anchor$locality,
    missing_anchor$depth_cm
  )
  if (any(missing_anchor$n_conserved < 1L)) {
    stop("Missing local x depth Conserved reference for ", ontology)
  }

  profile_out <- sample_out[sample_out$used_for_inference, ] |>
    dplyr::group_by(profile_id, locality, restoration4, ontology) |>
    dplyr::summarise(
      mean_distance_to_Conserved = mean(distance_to_local_depth_Conserved),
      n_depths = dplyr::n_distinct(depth_cm),
      .groups = "drop"
    ) |>
    dplyr::arrange(locality, profile_id)
  if (any(profile_out$n_depths != 3L)) {
    stop("Reference distance profile mean lacks a depth for ", ontology)
  }
  profile_out$locality <- factor(profile_out$locality)
  profile_out$restoration4 <- droplevels(factor(profile_out$restoration4, levels = STAGES))
  xref <- model.matrix(~ locality + restoration4, data = profile_out)
  if (qr(xref)$rank < ncol(xref)) {
    stop("Reference-anchored stage model is rank deficient for ", ontology)
  }
  reference_profiles[[ontology]] <- profile_out

  control_profile <- permute::how(
    nperm = N_PERM,
    blocks = profile_out$locality,
    within = permute::Within(type = "free")
  )
  set.seed(SEED + match(ontology, names(primary_data)))
  profile_perm_raw <- permute::shuffleSet(
    nrow(profile_out), nset = N_PERM, control = control_profile
  )
  profile_perm <- as.matrix(unique(as.data.frame(profile_perm_raw)))
  observed <- partial_f_stage(profile_out$mean_distance_to_Conserved, profile_out)
  perm_f <- vapply(seq_len(nrow(profile_perm)), function(ii) {
    partial_f_stage(
      profile_out$mean_distance_to_Conserved[profile_perm[ii, ]],
      profile_out
    )[["f"]]
  }, numeric(1))
  reference_tests[[ontology]] <- data.frame(
    ontology = ontology,
    model = "profile_mean_distance_to_Conserved ~ locality + restoration4",
    stage_model = "factor_nonconserved_only",
    observed_f = observed[["f"]],
    df_num = observed[["df_num"]],
    df_den = observed[["df_den"]],
    p_value = (1 + sum(perm_f >= observed[["f"]])) / (1 + length(perm_f)),
    n_profiles = nrow(profile_out),
    n_permutations = nrow(profile_perm),
    minimum_attainable_p = 1 / (nrow(profile_perm) + 1),
    seed = SEED + match(ontology, names(primary_data)),
    anchor = "Conserved_same_locality_and_depth",
    conserved_anchors_in_inference = FALSE,
    selection_uses_HI_HFR_or_stage = FALSE
  )
}

reference_samples_tbl <- do.call(rbind, reference_samples)
reference_profiles_tbl <- do.call(rbind, reference_profiles)
reference_tests_tbl <- do.call(rbind, reference_tests)

write_tsv(
  reference_samples_tbl,
  file.path(DIRS["tables"], "602_reference_distance_by_sample.tsv")
)
write_tsv(
  reference_profiles_tbl,
  file.path(DIRS["tables"], "602_reference_distance_by_profile.tsv")
)
write_tsv(
  reference_tests_tbl,
  file.path(DIRS["tables"], "602_reference_anchored_stage_tests.tsv")
)

reference_plot <- ggplot2::ggplot(
  reference_profiles_tbl,
  ggplot2::aes(restoration4, mean_distance_to_Conserved, color = locality)
) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.08, height = 0),
    size = 2.5
  ) +
  ggplot2::stat_summary(
    ggplot2::aes(group = restoration4),
    fun = mean, geom = "crossbar", width = 0.55, color = "black"
  ) +
  ggplot2::facet_wrap(~ ontology, scales = "free_y") +
  ggplot2::labs(
    x = "Restoration stage (Conserved anchors excluded)",
    y = "Profile mean Aitchison distance to local-depth Conserved reference",
    color = "Locality"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    legend.position = "bottom"
  )
ggplot2::ggsave(
  file.path(DIRS["figures"], "602_reference_anchored_functional_reassembly.png"),
  reference_plot, width = 10, height = 5.2, dpi = 400, bg = "white"
)
ggplot2::ggsave(
  file.path(DIRS["figures"], "602_reference_anchored_functional_reassembly.pdf"),
  reference_plot, width = 10, height = 5.2, bg = "white"
)

cat("[6/7] Saving figures and analysis objects\n")

for (ontology in unique(pcoa$ontology)) {
  plot_pcoa(
    pcoa[pcoa$ontology == ontology, ],
    ontology
  )
}

saveRDS(
  primary_data,
  file.path(
    DIRS["rds"],
    "602_primary_preprocessed_data.rds"
  ),
  compress = "xz"
)

saveRDS(
  list(
    metadata = md,
    permutations = perms,
    specifications = specs,
    permanova = permanova,
    betadisper = dispersion,
    primary = primary
  ),
  file.path(
    DIRS["rds"],
    "602_analysis_results.rds"
  ),
  compress = "xz"
)

cat("[7/7] Writing acceptance checks and provenance\n")

checks <- data.frame(
  check = c(
    "samples_51",
    "profiles_17",
    "all_profiles_complete",
    "all_locality_stage_cells_present_informative",
    "KO_IDs_match",
    "PFAM_IDs_match",
    "outcome_independent_filter",
    "unique_permutations_available",
    "all_permutations_valid",
    "primary_KO_result",
    "primary_PFAM_result"
  ),
  pass = c(
    nrow(md) == 51L,
    length(idx) == 17L,
    all(profile_ok),
    all(
      table(
        md$locality,
        md$restoration4
      ) > 0
    ),
    setequal(
      md$sample_id,
      rownames(raw$KEGG_ko)
    ),
    setequal(
      md$sample_id,
      rownames(raw$PFAM)
    ),
    TRUE,
    nrow(perms) > 0L,
    all(perm_audit$valid),
    sum(primary$ontology == "KEGG_ko") == 1L,
    sum(primary$ontology == "PFAM") == 1L
  ),
  critical = c(
    TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
  )
)

write_tsv(
  checks,
  file.path(
    DIRS["tables"],
    "602_acceptance_checks.tsv"
  )
)

if (!all(checks$pass[checks$critical])) {
  stop("Critical acceptance check failed")
}

write_tsv(
  data.frame(
    status = "PASS",
    samples = nrow(md),
    profiles = length(idx),
    complete_KO = ncol(raw$KEGG_ko),
    complete_PFAM = ncol(raw$PFAM),
    primary_KO = ncol(
      filtered$KEGG_ko$common
    ),
    primary_PFAM = ncol(
      filtered$PFAM$common
    ),
    permutations = nrow(perms),
    seed = SEED,
    primary_layer = "KEGG_ko",
    confirmation_layer = "PFAM",
    inference_unit = "complete sediment profile",
    permutation_scheme = paste(
      "profiles within locality;",
      "depth order preserved"
    ),
    transformation = "CLR",
    pseudocount = PC_PRIMARY,
    scope = paste(
      "functional gene potential;",
      "independent global analysis"
    )
  ),
  file.path(
    DIRS["tables"],
    "602_run_summary.tsv"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    DIRS["provenance"],
    "602_sessionInfo.txt"
  )
)

arg <- grep(
  "^--file=",
  commandArgs(FALSE),
  value = TRUE
)

if (length(arg)) {
  source_script <- sub(
    "^--file=",
    "",
    arg[1]
  )

  if (file.exists(source_script)) {
    file.copy(
      source_script,
      DIRS["provenance"],
      overwrite = TRUE
    )
  }
}

write_tsv(
  data.frame(
    input = names(input_files),
    path = unname(normalizePath(input_files, mustWork = TRUE)),
    md5 = unname(tools::md5sum(input_files))
  ),
  file.path(DIRS["provenance"], "602_input_files.tsv")
)

writeLines(
  OUT,
  file.path(
    RESULTS,
    "LATEST_602_validar_composicion_funcional_splitplot.txt"
  )
)

cat("\nPASS\n")
cat("Primary summary: tables/602_primary_summary.tsv\n")
cat("All PERMANOVA: tables/602_PERMANOVA_results.tsv\n")
cat("BETADISPER: tables/602_BETADISPER_results.tsv\n")
cat("Figures: figures/602_*_primary_Aitchison_PCoA.*\n")
cat("Finished UTC:", format(Sys.time(), tz = "UTC"), "\n")

sink(type = "message")
sink()
close(log_con)
