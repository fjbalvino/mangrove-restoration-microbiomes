#!/usr/bin/env Rscript
# ============================================================
# 03_03_diversity_composition.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Hill diversity, profile permutations, complete and selected taxonomic composition; Fig. 1b,c.
# Inputs (source expressions; complete list in docs/contracts/03_03_diversity_composition.json):
#   out <- trimws(readLines(path, n = 1L, warn = FALSE))
#   ps <- readRDS(PS_RDS)
#   external_meta <- read_csv(META_CSV, show_col_types = FALSE)
#   membership <- read_csv(MEMBERSHIP_CSV, show_col_types = FALSE) %>%
#   raw_metrics <- read_csv(METRICS_CSV, show_col_types = FALSE)
# Outputs (source expressions; complete list in contract):
#   writeLines(OUT_DIR, LATEST_ATTEMPT_FILE)
#   write_csv(parameter_gate, file.path(TABLE_DIR, "301_compuerta_parametros.csv"))
#   write_csv(input_gate, file.path(TABLE_DIR, "301_compuerta_inputs.csv"))
#   write_csv(
#   write_csv(profile_audit, file.path(TABLE_DIR, "301_auditoria_perfiles.csv"))
#   write_csv(design_gate, file.path(TABLE_DIR, "301_compuerta_diseno.csv"))
#   write_csv(lopo_gate, file.path(TABLE_DIR, "301_compuerta_LOPO.csv"))
#   ggsave(
#   write_csv(output_gate, file.path(TABLE_DIR, "301_compuerta_outputs.csv"))
#   write_csv(input_md5, file.path(TABLE_DIR, "301_input_md5.csv"))
# Algorithmic provenance:
# Hill diversity, Aitchison composition and profile-restricted permutation tests; dispersion.
#   McArdle & Anderson (2001), doi:10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2; Peres-Neto et al. (2006), doi:10.1890/0012-9658(2006)87[2614:VPOESD]2.0.CO;2.
# Source SHA-256: 6c0e4b9eaea616a23cc3aec68e57d4ac35877189fa35b10790c1469736aff78a
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


# 301_permanova_pcoa_restauracion_splitplot.R
#
# Primer bloque inferencial taxonomico del manuscrito.
#
# Orden obligatorio del manuscrito:
#   1) comunidad completa: diversidad Hill y ordenamiento/PERMANOVA;
#   2) comparacion de 1,031 candidatos crudos no unimodales frente a 84
#      candidatos residuales no unimodales;
#   3) sensibilidad comun: seleccionar las coordenadas de ambos universos
#      desde un unico CLR de la comunidad completa y conservar las 51 muestras;
#   4) sensibilidad de influencia dejando fuera un perfil completo cada vez.
#
# La etapa de restauracion y la localidad pertenecen al perfil. Las tres
# profundidades (5, 20 y 40 cm) se mantienen juntas y en el mismo orden durante
# las permutaciones de perfiles dentro de localidad. La comunidad completa usa
# las 51 muestras canonicas. Las subcomposiciones se conservan solo como
# sensibilidad y excluyen perfiles enteros si alguna profundidad tiene suma
# cero. La comparacion de coordenadas del CLR completo usa las mismas 51
# muestras, los mismos perfiles y las mismas permutaciones para 1,031 y 84.

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(permute)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Funciones generales
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

unique_permutation_count <- function(p) {
  nrow(unique(as.data.frame(p)))
}

# -----------------------------------------------------------------------------
# Permutaciones del diseno split-plot
# -----------------------------------------------------------------------------

profile_permutations_rows <- function(meta, n_perm, seed) {
  control <- permute::how(
    nperm = n_perm,
    blocks = factor(meta$locality),
    plots = permute::Plots(
      strata = factor(meta$profile_id),
      type = "free"
    ),
    within = permute::Within(type = "none")
  )
  set.seed(seed)
  raw <- permute::shuffleSet(nrow(meta), nset = n_perm, control = control)
  as.matrix(unique(as.data.frame(raw)))
}

profile_permutations_profiles <- function(meta_profile, n_perm, seed) {
  control <- permute::how(
    nperm = n_perm,
    blocks = factor(meta_profile$locality),
    within = permute::Within(type = "free")
  )
  set.seed(seed)
  raw <- permute::shuffleSet(
    nrow(meta_profile),
    nset = n_perm,
    control = control
  )
  as.matrix(unique(as.data.frame(raw)))
}

audit_row_permutations <- function(p, meta) {
  groups <- split(seq_len(nrow(meta)), meta$profile_id)
  valid <- vapply(
    seq_len(nrow(p)),
    function(i) {
      all(vapply(
        groups,
        function(idx) {
          source_rows <- p[i, idx]
          n_distinct(meta$profile_id[source_rows]) == 1L &&
            meta$locality[source_rows[[1L]]] == meta$locality[idx[[1L]]] &&
            identical(
              as.character(meta$depth_cm[source_rows]),
              as.character(meta$depth_cm[idx])
            )
        },
        logical(1)
      ))
    },
    logical(1)
  )
  tibble(
    permutation = seq_along(valid),
    valid_whole_profile_within_locality = valid
  )
}

audit_profile_permutations <- function(p, meta_profile) {
  valid <- vapply(
    seq_len(nrow(p)),
    function(i) {
      all(
        as.character(meta_profile$locality[p[i, ]]) ==
          as.character(meta_profile$locality)
      )
    },
    logical(1)
  )
  tibble(
    permutation = seq_along(valid),
    valid_profile_within_locality = valid
  )
}

# -----------------------------------------------------------------------------
# Diversidad alfa y pruebas por perfil
# -----------------------------------------------------------------------------

hill_numbers <- function(counts_samples_taxa) {
  totals <- rowSums(counts_samples_taxa)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stop("Hill requiere bibliotecas positivas y finitas", call. = FALSE)
  }
  proportions <- counts_samples_taxa / totals
  log_proportions <- log(proportions)
  log_proportions[!is.finite(log_proportions)] <- 0

  tibble(
    sample_id = rownames(counts_samples_taxa),
    library_size = as.numeric(totals),
    q0_richness = rowSums(counts_samples_taxa > 0),
    q1_exp_shannon = exp(-rowSums(proportions * log_proportions)),
    q2_inverse_simpson = 1 / rowSums(proportions^2)
  )
}

mean_rarefied_hill <- function(counts, depth, n_iterations, seed) {
  metric_names <- c(
    "q0_richness",
    "q1_exp_shannon",
    "q2_inverse_simpson"
  )
  accumulator <- matrix(
    0,
    nrow = nrow(counts),
    ncol = length(metric_names),
    dimnames = list(rownames(counts), metric_names)
  )

  set.seed(seed)
  for (i in seq_len(n_iterations)) {
    rarefied <- suppressWarnings(
      vegan::rrarefy(counts, sample = depth)
    )
    h <- hill_numbers(rarefied)
    accumulator <- accumulator + as.matrix(h[, metric_names, drop = FALSE])
  }

  out <- accumulator / n_iterations
  as_tibble(out, rownames = "sample_id") %>%
    mutate(
      rarefaction_depth = depth,
      rarefaction_iterations = n_iterations
    )
}

partial_f_value <- function(data, response, reduced_rhs, full_rhs) {
  reduced <- lm(
    as.formula(paste(response, "~", reduced_rhs)),
    data = data
  )
  full <- lm(
    as.formula(paste(response, "~", full_rhs)),
    data = data
  )
  rss_reduced <- sum(residuals(reduced)^2)
  rss_full <- sum(residuals(full)^2)
  df_num <- df.residual(reduced) - df.residual(full)
  df_den <- df.residual(full)
  f_value <- ((rss_reduced - rss_full) / df_num) / (rss_full / df_den)
  c(F = f_value, df_num = df_num, df_den = df_den)
}

profile_alpha_test <- function(
  data,
  response,
  permutations,
  analysis_role,
  adjust_library = FALSE,
  seed
) {
  reduced_rhs <- if (adjust_library) {
    "locality + log10_library_size"
  } else {
    "locality"
  }
  full_rhs <- paste(reduced_rhs, "+ restoration4")

  observed <- partial_f_value(data, response, reduced_rhs, full_rhs)
  reduced_fit <- lm(
    as.formula(paste(response, "~", reduced_rhs)),
    data = data
  )
  fitted_reduced <- fitted(reduced_fit)
  residual_reduced <- residuals(reduced_fit)

  permuted_f <- vapply(
    seq_len(nrow(permutations)),
    function(i) {
      permuted_data <- data
      permuted_data[[response]] <-
        fitted_reduced + residual_reduced[permutations[i, ]]
      partial_f_value(
        permuted_data,
        response,
        reduced_rhs,
        full_rhs
      )[["F"]]
    },
    numeric(1)
  )
  permuted_f <- permuted_f[is.finite(permuted_f)]

  tibble(
    response = response,
    analysis_role = analysis_role,
    n_profiles = nrow(data),
    observed_F = observed[["F"]],
    df_num = observed[["df_num"]],
    df_den = observed[["df_den"]],
    p_value = (1 + sum(permuted_f >= observed[["F"]])) /
      (1 + length(permuted_f)),
    reduced_model = paste(response, "~", reduced_rhs),
    full_model = paste(response, "~", full_rhs),
    permutation_method = "Freedman-Lane_residuals_within_locality",
    n_unique_permutations = unique_permutation_count(permutations),
    minimum_attainable_p = 1 /
      (unique_permutation_count(permutations) + 1),
    seed = seed
  )
}

# -----------------------------------------------------------------------------
# PERMANOVA, dispersion y PCoA
# -----------------------------------------------------------------------------

inference_role <- function(community, metric, model_role) {
  dplyr::case_when(
    model_role == "LOPO_profile_influence_sensitivity" ~
      "LOPO_profile_influence_sensitivity",
    model_role == "ECI_adjusted_sensitivity" ~
      "ECI_adjusted_sensitivity",
    community == "01_full_species_community_all51" &
      metric == "Aitchison_primary" ~
      "primary_full_community",
    community == "01_full_species_community_all51" &
      metric == "Bray_Hellinger_sensitivity" ~
      "metric_sensitivity_full_community",
    community == "02_full_species_community_matched_candidates" ~
      "population_control_full_community",
    community == "03_raw_nonunimodal_candidates_1031_available" ~
      "raw1031_available_population_sensitivity",
    community == "04_raw_nonunimodal_candidates_1031_matched" ~
      "candidate_universe_comparison_raw1031",
    community == "05_residual_nonunimodal_candidates_84_matched" ~
      "candidate_universe_comparison_residual84",
    community == "06_raw1031_selected_from_full_CLR_all51" ~
      "common_full_CLR_comparison_raw1031",
    community == "07_residual84_selected_from_full_CLR_all51" ~
      "common_full_CLR_comparison_residual84",
    TRUE ~ "other_sensitivity"
  )
}

inference_role_self_test <- inference_role(
  community = c(
    "01_full_species_community_all51",
    "01_full_species_community_all51",
    "02_full_species_community_matched_candidates",
    "03_raw_nonunimodal_candidates_1031_available",
    "04_raw_nonunimodal_candidates_1031_matched",
    "05_residual_nonunimodal_candidates_84_matched"
  ),
  metric = c(
    "Aitchison_primary",
    "Bray_Hellinger_sensitivity",
    "Aitchison_primary",
    "Aitchison_primary",
    "Aitchison_primary",
    "Aitchison_primary"
  ),
  model_role = c(
    "primary_adjusted_locality_depth",
    "primary_adjusted_locality_depth",
    "primary_adjusted_locality_depth",
    "primary_adjusted_locality_depth",
    "primary_adjusted_locality_depth",
    "primary_adjusted_locality_depth"
  )
)
if (!identical(
  inference_role_self_test,
  c(
    "primary_full_community",
    "metric_sensitivity_full_community",
    "population_control_full_community",
    "raw1031_available_population_sensitivity",
    "candidate_universe_comparison_raw1031",
    "candidate_universe_comparison_residual84"
  )
)) {
  stop("Fallo el auto-test vectorial de inference_role", call. = FALSE)
}

run_adonis <- function(
  distance_object,
  meta,
  permutations,
  community,
  community_order,
  metric,
  model_role,
  include_eci,
  seed
) {
  rhs <- if (include_eci) {
    "locality + depth_cm + ECI_PC1 + restoration4"
  } else {
    "locality + depth_cm + restoration4"
  }
  model_formula <- as.formula(
    paste("distance_object ~", rhs),
    env = environment()
  )
  fit <- vegan::adonis2(
    model_formula,
    data = meta,
    permutations = permutations,
    by = "margin"
  )
  unique_n <- unique_permutation_count(permutations)

  as.data.frame(fit) %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
    rename(
      df = Df,
      sum_squares = SumOfSqs,
      R2 = R2,
      pseudo_F = F,
      p_value_under_permutation_scheme = `Pr(>F)`
    ) %>%
    mutate(
      p_value = if_else(
        term == "restoration4",
        p_value_under_permutation_scheme,
        NA_real_
      ),
      p_value_valid = term == "restoration4",
      community = community,
      community_order = community_order,
      metric = metric,
      model_role = model_role,
      inference_role = inference_role(community, metric, model_role),
      analysis_population = unique(meta$analysis_population),
      n_samples_analysis = nrow(meta),
      n_profiles_analysis = n_distinct(meta$profile_id),
      model = paste("distance ~", rhs),
      permutation_scheme =
        "whole_profiles_within_locality_depth_order_preserved",
      n_requested_permutations = nrow(permutations),
      n_unique_permutations = unique_n,
      minimum_attainable_p = 1 / (unique_n + 1),
      seed = seed,
      interpretation = if_else(
        term == "restoration4",
        "restoration_stage_test",
        "nuisance_term_p_not_interpreted"
      )
    )
}

run_dispersion <- function(
  distance_object,
  meta,
  permutations,
  community,
  community_order,
  metric,
  seed
) {
  bd <- vegan::betadisper(
    distance_object,
    meta$restoration4,
    type = "median",
    bias.adjust = TRUE,
    add = "lingoes"
  )
  tab <- as.data.frame(
    vegan::permutest(bd, permutations = permutations)$tab
  ) %>%
    rownames_to_column("term") %>%
    as_tibble()
  names(tab) <- make.names(names(tab))
  group_row <- tab %>% slice(1L)
  unique_n <- unique_permutation_count(permutations)

  tibble(
    community = community,
    community_order = community_order,
    metric = metric,
    analysis_population = unique(meta$analysis_population),
    n_samples_analysis = nrow(meta),
    n_profiles_analysis = n_distinct(meta$profile_id),
    grouping_factor = "restoration4",
    df = group_row$Df,
    sum_squares = group_row$Sum.Sq,
    mean_squares = group_row$Mean.Sq,
    pseudo_F = group_row$F,
    p_value = group_row$Pr..F.,
    n_requested_permutations = nrow(permutations),
    n_unique_permutations = unique_n,
    minimum_attainable_p = 1 / (unique_n + 1),
    seed = seed,
    permutation_scheme =
      "whole_profiles_within_locality_depth_order_preserved",
    interpretation =
      "diagnostic_for_stage_location_dispersion_confounding"
  )
}

pcoa_scores <- function(
  distance_object,
  meta,
  community,
  community_order,
  metric
) {
  fit <- vegan::wcmdscale(
    distance_object,
    k = 2L,
    eig = TRUE,
    add = "lingoes"
  )
  points <- fit$points
  if (!identical(rownames(points), meta$sample_id)) {
    points <- points[match(meta$sample_id, rownames(points)), , drop = FALSE]
  }
  if (anyNA(points)) {
    stop("PCoA no pudo alinearse por sample_id", call. = FALSE)
  }
  positive_eigenvalues <- fit$eig[fit$eig > 0]
  variance <- if (length(positive_eigenvalues)) {
    fit$eig[seq_len(2L)] / sum(positive_eigenvalues)
  } else {
    c(NA_real_, NA_real_)
  }

  bind_cols(
    meta,
    tibble(PCoA1 = points[, 1L], PCoA2 = points[, 2L])
  ) %>%
    mutate(
      community = community,
      community_order = community_order,
      metric = metric,
      PCoA1_variance_positive = variance[[1L]],
      PCoA2_variance_positive = variance[[2L]],
      correction = "Lingoes"
    )
}

run_community <- function(
  counts,
  meta,
  permutations,
  community,
  community_order,
  pseudocount,
  seed
) {
  cat(
    "Community analysis start:", community,
    "| samples:", nrow(meta),
    "| profiles:", n_distinct(meta$profile_id),
    "\n"
  )
  flush(log_con)
  if (!identical(rownames(counts), meta$sample_id)) {
    stop(
      "Orden de counts y metadata inconsistente para ", community,
      call. = FALSE
    )
  }
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("Conteos no finitos o negativos para ", community, call. = FALSE)
  }
  if (any(rowSums(counts) <= 0)) {
    stop(
      "El universo ", community,
      " conserva muestras con suma cero tras aplicar su poblacion analitica",
      call. = FALSE
    )
  }
  aitchison <- stats::dist(clr_rows(counts, pseudocount))
  hellinger <- sqrt(counts / rowSums(counts))
  if (any(!is.finite(hellinger))) {
    stop("Hellinger no finito para ", community, call. = FALSE)
  }
  bray <- vegan::vegdist(hellinger, method = "bray")

  permanova <- bind_rows(
    run_adonis(
      aitchison,
      meta,
      permutations,
      community,
      community_order,
      "Aitchison_primary",
      "primary_adjusted_locality_depth",
      FALSE,
      seed
    ),
    run_adonis(
      bray,
      meta,
      permutations,
      community,
      community_order,
      "Bray_Hellinger_sensitivity",
      "primary_adjusted_locality_depth",
      FALSE,
      seed
    ),
    run_adonis(
      aitchison,
      meta,
      permutations,
      community,
      community_order,
      "Aitchison_primary",
      "ECI_adjusted_sensitivity",
      TRUE,
      seed
    ),
    run_adonis(
      bray,
      meta,
      permutations,
      community,
      community_order,
      "Bray_Hellinger_sensitivity",
      "ECI_adjusted_sensitivity",
      TRUE,
      seed
    )
  )

  dispersion <- bind_rows(
    run_dispersion(
      aitchison,
      meta,
      permutations,
      community,
      community_order,
      "Aitchison_primary",
      seed
    ),
    run_dispersion(
      bray,
      meta,
      permutations,
      community,
      community_order,
      "Bray_Hellinger_sensitivity",
      seed
    )
  )

  scores <- bind_rows(
    pcoa_scores(
      aitchison,
      meta,
      community,
      community_order,
      "Aitchison_primary"
    ),
    pcoa_scores(
      bray,
      meta,
      community,
      community_order,
      "Bray_Hellinger_sensitivity"
    )
  )

  cat("Community analysis complete:", community, "\n")
  flush(log_con)

  list(
    permanova = permanova,
    dispersion = dispersion,
    scores = scores,
    distance_finite = all(is.finite(as.vector(aitchison))) &&
      all(is.finite(as.vector(bray)))
  )
}

# Analisis de una distancia euclidiana ya construida. Se usa para seleccionar
# coordenadas de candidatos desde un unico CLR de la comunidad completa. Esta
# distancia no es la Aitchison de una subcomposicion recalculada: cuantifica el
# cambio de los candidatos con respecto a la media geometrica de la comunidad
# completa y conserva las 51 muestras, incluso si un subconjunto suma cero.
run_precomputed_distance <- function(
  distance_object,
  meta,
  permutations,
  community,
  community_order,
  metric,
  seed
) {
  if (!identical(attr(distance_object, "Labels"), meta$sample_id)) {
    stop("Distancia precomputada y metadata no estan alineadas", call. = FALSE)
  }
  if (any(!is.finite(as.vector(distance_object)))) {
    stop("Distancia precomputada no finita", call. = FALSE)
  }

  permanova <- bind_rows(
    run_adonis(
      distance_object,
      meta,
      permutations,
      community,
      community_order,
      metric,
      "primary_adjusted_locality_depth",
      FALSE,
      seed
    ),
    run_adonis(
      distance_object,
      meta,
      permutations,
      community,
      community_order,
      metric,
      "ECI_adjusted_sensitivity",
      TRUE,
      seed
    )
  )

  dispersion <- run_dispersion(
    distance_object,
    meta,
    permutations,
    community,
    community_order,
    metric,
    seed
  )

  scores <- pcoa_scores(
    distance_object,
    meta,
    community,
    community_order,
    metric
  )

  list(
    permanova = permanova,
    dispersion = dispersion,
    scores = scores,
    distance_finite = TRUE
  )
}

subset_distance <- function(distance_object, sample_ids) {
  labels <- attr(distance_object, "Labels")
  if (!all(sample_ids %in% labels)) {
    stop("La subdistancia LOPO contiene IDs desconocidos", call. = FALSE)
  }
  stats::as.dist(
    as.matrix(distance_object)[sample_ids, sample_ids, drop = FALSE]
  )
}

run_lopo_stage <- function(
  distance_object,
  meta,
  analysis,
  community_order,
  n_perm,
  seed
) {
  profiles <- unique(as.character(meta$profile_id))

  bind_rows(lapply(seq_along(profiles), function(i) {
    omitted <- profiles[[i]]
    cat(
      "LOPO:", analysis, "|", i, "/", length(profiles),
      "| omitted profile:", omitted, "\n"
    )
    flush(log_con)
    meta_i <- meta %>%
      filter(as.character(profile_id) != omitted) %>%
      droplevels() %>%
      arrange(locality, profile_id, depth_cm) %>%
      mutate(
        analysis_population = paste0(
          "LOPO_", n_distinct(profile_id), "_profiles"
        )
      )
    distance_i <- subset_distance(distance_object, meta_i$sample_id)
    permutations_i <- profile_permutations_rows(
      meta_i,
      n_perm,
      seed + i
    )
    audit_i <- audit_row_permutations(permutations_i, meta_i)
    if (!nrow(permutations_i) ||
        !all(audit_i$valid_whole_profile_within_locality)) {
      stop("Fallo la auditoria LOPO de perfil completo", call. = FALSE)
    }

    stage_i <- run_adonis(
      distance_i,
      meta_i,
      permutations_i,
      analysis,
      community_order,
      "LOPO_precomputed_distance",
      "LOPO_profile_influence_sensitivity",
      FALSE,
      seed + i
    ) %>%
      filter(term == "restoration4")

    omitted_meta <- meta %>%
      filter(as.character(profile_id) == omitted) %>%
      summarise(
        omitted_locality = first(as.character(locality)),
        omitted_stage = first(as.character(restoration4))
      )

    stage_i %>%
      transmute(
        analysis = analysis,
        omitted_profile = omitted,
        omitted_locality = omitted_meta$omitted_locality,
        omitted_stage = omitted_meta$omitted_stage,
        n_samples = nrow(meta_i),
        n_profiles = n_distinct(meta_i$profile_id),
        R2 = R2,
        pseudo_F = pseudo_F,
        p_value = p_value,
        n_unique_permutations = n_unique_permutations,
        minimum_attainable_p = minimum_attainable_p,
        all_permutations_preserve_whole_profiles =
          all(audit_i$valid_whole_profile_within_locality)
      )
  }))
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
RUN003 <- get_arg(
  "--run003",
  read_latest(OUT_ROOT, "003_integrar_ejes_ECI_HI_en_metadata")
)
RUN201 <- get_arg(
  "--run201",
  read_latest(OUT_ROOT, "201_detectar_especies_bimodales_bootstrap")
)
RUN202 <- get_arg(
  "--run202",
  read_latest(OUT_ROOT, "202_construir_matriz_CLR_y_eje_bimodal")
)

PS_RDS <- get_arg(
  "--ps_rds",
  file.path(RUN001, "rds", "001_phyloseq_canon_51.rds")
)
META_CSV <- get_arg(
  "--meta_csv",
  file.path(RUN003, "tables", "003_metadata_integrada_canon_51.csv")
)
MEMBERSHIP_CSV <- get_arg(
  "--membership_csv",
  file.path(RUN202, "tables", "202_membresia_universos_84_12_3.csv")
)
METRICS_CSV <- get_arg(
  "--metrics_csv",
  file.path(RUN201, "tables", "201_metricas_bimodalidad_todas_especies.csv")
)

N_PERM <- as_int(get_arg("--n_perm", "999"), "n_perm")
N_RAREFY <- as_int(get_arg("--n_rarefy", "100"), "n_rarefy")
SEED <- as_int(get_arg("--seed", "123"), "seed")
PSEUDOCOUNT <- as_num(get_arg("--pseudocount", "1"), "pseudocount")
EXPECTED_FULL <- as_int(get_arg("--expected_full", "12761"), "expected_full")
EXPECTED_RAW_CANDIDATES <- as_int(
  get_arg("--expected_raw_candidates", "1031"),
  "expected_raw_candidates"
)
EXPECTED_CANDIDATES <- as_int(
  get_arg("--expected_candidates", "84"),
  "expected_candidates"
)

SCRIPT_ID <- "301_permanova_pcoa_restauracion_splitplot"
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

log_con <- file(file.path(LOG_DIR, "301_log.txt"), "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0L) sink(type = "message")
  while (sink.number() > 0L) sink()
  close(log_con)
}, add = TRUE)

cat("===== 301 Full community first; candidate taxa second =====\n")
cat("Start UTC:", format(Sys.time(), tz = "UTC"), "\n")
cat("Phyloseq:", PS_RDS, "\n")
cat("Metadata:", META_CSV, "\n")
cat("RUN201:", RUN201, "\n")
cat("RUN202:", RUN202, "\n")
cat("Output:", OUT_DIR, "\n")

parameter_gate <- tibble(
  gate = c(
    "n_perm_at_least_999",
    "n_rarefy_at_least_100",
    "pseudocount_positive",
    "expected_universe_sizes_valid"
  ),
  pass = c(
    N_PERM >= 999L,
    N_RAREFY >= 100L,
    PSEUDOCOUNT > 0,
    EXPECTED_FULL >= EXPECTED_RAW_CANDIDATES &&
      EXPECTED_FULL >= EXPECTED_CANDIDATES
  ),
  detail = c(
    as.character(N_PERM),
    as.character(N_RAREFY),
    as.character(PSEUDOCOUNT),
    paste(
      EXPECTED_FULL,
      EXPECTED_RAW_CANDIDATES,
      EXPECTED_CANDIDATES,
      sep = "|"
    )
  )
)
write_csv(parameter_gate, file.path(TABLE_DIR, "301_compuerta_parametros.csv"))
if (!all(parameter_gate$pass)) {
  stop("Fallo la compuerta de parametros de 301", call. = FALSE)
}

PS_RDS <- require_file(PS_RDS)
META_CSV <- require_file(META_CSV)
MEMBERSHIP_CSV <- require_file(MEMBERSHIP_CSV)
METRICS_CSV <- require_file(METRICS_CSV)

# -----------------------------------------------------------------------------
# Carga y compuertas de inputs
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

external_meta <- read_csv(META_CSV, show_col_types = FALSE)
if (!"sample_id" %in% names(external_meta)) {
  stop("meta_csv requiere sample_id", call. = FALSE)
}
external_meta <- external_meta %>%
  mutate(sample_id = as.character(sample_id))

input_gate <- tibble(
  gate = c(
    "phyloseq_has_51_samples",
    "phyloseq_has_expected_full_taxa",
    "sample_data_existing_id_agrees_with_rownames",
    "otu_and_sample_data_ids_match",
    "external_metadata_has_51_unique_ids",
    "phyloseq_and_external_metadata_ids_match",
    "taxon_ids_unique",
    "otu_values_finite",
    "otu_values_nonnegative",
    "otu_values_integer_like",
    "all_samples_have_positive_library_size",
    "all_taxa_have_positive_total_count"
  ),
  pass = c(
    ncol(otu) == 51L && nrow(ps_meta) == 51L,
    nrow(otu) == EXPECTED_FULL,
    isTRUE(ps_meta_id_agreement),
    !anyDuplicated(colnames(otu)) &&
      !anyDuplicated(ps_meta$sample_id) &&
      setequal(colnames(otu), ps_meta$sample_id),
    nrow(external_meta) == 51L &&
      !anyNA(external_meta$sample_id) &&
      !anyDuplicated(external_meta$sample_id),
    setequal(ps_meta$sample_id, external_meta$sample_id),
    !is.null(rownames(otu)) && !anyDuplicated(rownames(otu)),
    all(is.finite(otu)),
    all(is.finite(otu)) && all(otu >= 0),
    all(is.finite(otu)) && all(abs(otu - round(otu)) < 1e-8),
    all(is.finite(otu)) && all(colSums(otu) > 0),
    all(is.finite(otu)) && all(rowSums(otu) > 0)
  ),
  detail = c(
    paste0("otu=", ncol(otu), "; sample_data=", nrow(ps_meta)),
    paste0("n_taxa=", nrow(otu)),
    if (ps_meta_has_sample_id) {
      paste0("column_present; exact_ordered_match=", ps_meta_id_agreement)
    } else {
      "column_absent; phyloseq rownames used"
    },
    paste0("otu_ids=", ncol(otu), "; metadata_ids=", nrow(ps_meta)),
    paste0(
      "n=", nrow(external_meta),
      "; duplicated=", anyDuplicated(external_meta$sample_id)
    ),
    paste0(
      "shared=",
      length(intersect(ps_meta$sample_id, external_meta$sample_id))
    ),
    paste0("n_taxa=", nrow(otu)),
    paste0("nonfinite=", sum(!is.finite(otu))),
    paste0("negative=", sum(otu < 0, na.rm = TRUE)),
    "raw count matrix required",
    paste0("empty_samples=", sum(colSums(otu) <= 0)),
    paste0("empty_taxa=", sum(rowSums(otu) <= 0))
  )
)
write_csv(input_gate, file.path(TABLE_DIR, "301_compuerta_inputs.csv"))
if (!all(input_gate$pass)) {
  stop("Fallo la compuerta de inputs de 301", call. = FALSE)
}

# -----------------------------------------------------------------------------
# Metadata canonica y estructura de perfiles
# -----------------------------------------------------------------------------

meta0 <- inner_join(
  ps_meta,
  external_meta,
  by = "sample_id",
  suffix = c(".ps", ".meta")
)

STAGES <- c(
  "Degraded",
  "Early restoration",
  "Intermediate restoration",
  "Advanced restoration",
  "Conserved"
)
DEPTHS <- c(5, 20, 40)

# Paleta original del manuscrito y de los scripts 011/051.
STAGE_COLORS <- c(
  "Degraded" = "#E41A1C",
  "Early restoration" = "#FDAE61",
  "Intermediate restoration" = "#FFF7A8",
  "Advanced restoration" = "#A6DDA0",
  "Conserved" = "#2B83BA"
)

STAGE_LABELS_SHORT <- c(
  "Degraded" = "D",
  "Early restoration" = "Early",
  "Intermediate restoration" = "Inter.",
  "Advanced restoration" = "Adv.",
  "Conserved" = "Cons."
)

stage_palette_contract <- tibble(
  manuscript_order = seq_along(STAGES),
  restoration4 = STAGES,
  short_label = unname(STAGE_LABELS_SHORT[STAGES]),
  hex_color = unname(STAGE_COLORS[STAGES]),
  source = "original manuscript scripts 011 and 051"
)
write_csv(
  stage_palette_contract,
  file.path(TABLE_DIR, "301_contrato_orden_y_colores_etapas.csv")
)

meta <- meta0 %>%
  transmute(
    sample_id = as.character(sample_id),
    locality = str_squish(
      coalesce_col(meta0, c("locality.meta", "locality", "locality.ps"))
    ),
    restoration4 = str_squish(
      coalesce_col(
        meta0,
        c(
          "restoration4.meta",
          "restoration4",
          "restoration4.ps",
          "collapsed_stage.meta",
          "collapsed_stage.ps"
        )
      )
    ),
    depth_cm = normalize_depth(
      coalesce_col(meta0, c("depth_cm.meta", "depth_cm", "depth_cm.ps"))
    ),
    lat_block = str_squish(
      coalesce_col(meta0, c("lat_block.meta", "lat_block", "lat_block.ps"))
    ),
    ECI_PC1 = suppressWarnings(
      as.numeric(
        coalesce_col(meta0, c("ECI_PC1.meta", "ECI_PC1", "ECI_PC1.ps"))
      )
    )
  ) %>%
  mutate(
    restoration4 = recode(restoration4, Preserved = "Conserved"),
    restoration4 = factor(restoration4, levels = STAGES, ordered = FALSE),
    locality = factor(locality),
    depth_cm = factor(depth_cm, levels = DEPTHS),
    profile_id = if_else(
      str_detect(lat_block, fixed(as.character(locality))),
      lat_block,
      paste(as.character(locality), lat_block, sep = "::")
    )
  ) %>%
  filter(sample_id %in% colnames(otu)) %>%
  arrange(locality, profile_id, depth_cm)

profile_audit <- meta %>%
  group_by(profile_id) %>%
  summarise(
    locality = first(as.character(locality)),
    restoration4 = first(as.character(restoration4)),
    n = n(),
    n_localities = n_distinct(locality),
    n_stages = n_distinct(restoration4),
    depths = paste(
      sort(as.numeric(as.character(depth_cm))),
      collapse = "|"
    ),
    valid = n == 3L &&
      n_localities == 1L &&
      n_stages == 1L &&
      setequal(as.numeric(as.character(depth_cm)), DEPTHS),
    .groups = "drop"
  )
write_csv(profile_audit, file.path(TABLE_DIR, "301_auditoria_perfiles.csv"))

X_primary <- model.matrix(
  ~ locality + depth_cm + restoration4,
  data = meta
)
X_eci <- model.matrix(
  ~ locality + depth_cm + ECI_PC1 + restoration4,
  data = meta
)

design_gate <- tibble(
  gate = c(
    "canonical_51_samples",
    "canonical_17_complete_profiles",
    "canonical_three_depths",
    "canonical_three_localities",
    "canonical_five_restoration_stages",
    "ECI_PC1_complete_and_variable",
    "primary_model_full_rank",
    "ECI_sensitivity_model_full_rank"
  ),
  pass = c(
    nrow(meta) == 51L && !anyDuplicated(meta$sample_id),
    nrow(profile_audit) == 17L && all(profile_audit$valid),
    setequal(as.numeric(as.character(meta$depth_cm)), DEPTHS),
    n_distinct(meta$locality) == 3L,
    setequal(as.character(meta$restoration4), STAGES),
    all(is.finite(meta$ECI_PC1)) && stats::sd(meta$ECI_PC1) > 0,
    qr(X_primary)$rank == ncol(X_primary),
    qr(X_eci)$rank == ncol(X_eci)
  ),
  detail = c(
    paste0("n=", nrow(meta)),
    paste0(
      "profiles=", nrow(profile_audit),
      "; valid=", sum(profile_audit$valid)
    ),
    paste(sort(unique(as.numeric(as.character(meta$depth_cm)))), collapse = ","),
    paste(levels(meta$locality), collapse = ","),
    paste(levels(meta$restoration4), collapse = " -> "),
    paste0(
      "finite=", sum(is.finite(meta$ECI_PC1)),
      "; sd=", signif(stats::sd(meta$ECI_PC1), 6)
    ),
    "distance ~ locality + depth_cm + restoration4",
    "distance ~ locality + depth_cm + ECI_PC1 + restoration4"
  )
)
write_csv(design_gate, file.path(TABLE_DIR, "301_compuerta_diseno.csv"))
if (!all(design_gate$pass)) {
  stop("Fallo la compuerta de diseno de 301", call. = FALSE)
}

otu <- otu[, meta$sample_id, drop = FALSE]

# -----------------------------------------------------------------------------
# Auditoria de linajes macroeucariotas
# -----------------------------------------------------------------------------

tax_object <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
if (is.null(tax_object)) {
  stop("El phyloseq requiere tax_table para auditar macroeucariotas", call. = FALSE)
}
tax_df_raw <- as.data.frame(tax_object, stringsAsFactors = FALSE)
class(tax_df_raw) <- "data.frame"
tax_df <- data.frame(
  Original_Taxon_ID = as.character(rownames(tax_df_raw)),
  tax_df_raw,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = NULL
)
tax_df <- tax_df[match(rownames(otu), tax_df$Original_Taxon_ID), , drop = FALSE]
if (anyNA(tax_df$Original_Taxon_ID)) {
  stop("Taxonomia incompleta para el universo phyloseq", call. = FALSE)
}

rank_columns <- setdiff(names(tax_df), "Original_Taxon_ID")
taxonomy_text <- apply(
  tax_df[, rank_columns, drop = FALSE],
  1L,
  function(z) paste(as.character(z), collapse = ";")
)
macro_eukaryote <- str_detect(
  taxonomy_text,
  regex(
    "(^|;|__)Homo sapiens($|;)|(^|;|__)(Metazoa|Animalia|Plantae|Viridiplantae)($|;)",
    ignore_case = TRUE
  )
)
macro_eukaryote[is.na(macro_eukaryote)] <- FALSE

tax_exclusion_audit <- tibble(
  Original_Taxon_ID = rownames(otu),
  taxonomy_text = taxonomy_text,
  excluded_macro_eukaryote = macro_eukaryote,
  exclusion_rule =
    "exclude Homo sapiens, Metazoa/Animalia and Plantae/Viridiplantae; retain Fungi"
)
write_csv(
  tax_exclusion_audit,
  file.path(TABLE_DIR, "301_auditoria_exclusion_macroeucariotas.csv")
)

full_taxa <- rownames(otu)[!macro_eukaryote]
if (length(full_taxa) < 2L) {
  stop("Muy pocos taxa despues de excluir macroeucariotas", call. = FALSE)
}
counts_full <- t(otu[full_taxa, , drop = FALSE])

# -----------------------------------------------------------------------------
# Universos crudo 1,031 y residual 84
# -----------------------------------------------------------------------------

membership <- read_csv(MEMBERSHIP_CSV, show_col_types = FALSE) %>%
  mutate(
    Original_Taxon_ID = as.character(Original_Taxon_ID),
    candidate84 = as_logical_strict(candidate84),
    in_phyloseq = as_logical_strict(in_phyloseq)
  )
raw_metrics <- read_csv(METRICS_CSV, show_col_types = FALSE)

required_raw_columns <- c("Original_Taxon_ID", "observed_joint_raw")
missing_raw_columns <- setdiff(required_raw_columns, names(raw_metrics))
if (length(missing_raw_columns)) {
  stop(
    "Faltan columnas en metricas 201: ",
    paste(missing_raw_columns, collapse = ", "),
    call. = FALSE
  )
}

raw_metrics <- raw_metrics %>%
  mutate(
    Original_Taxon_ID = as.character(Original_Taxon_ID),
    observed_joint_raw = as_logical_strict(observed_joint_raw)
  )

raw_candidate_ids <- raw_metrics %>%
  filter(observed_joint_raw %in% TRUE) %>%
  pull(Original_Taxon_ID) %>%
  unique()

candidate_ids <- membership %>%
  filter(candidate84 %in% TRUE) %>%
  pull(Original_Taxon_ID)

universe_gate <- tibble(
  gate = c(
    "raw_metrics_ids_unique",
    "raw1031_exact_count",
    "raw1031_all_map_to_phyloseq",
    "raw1031_all_retained_after_macro_eukaryote_exclusion",
    "membership_ids_unique",
    "candidate84_exact_count",
    "all_selected_taxa_in_phyloseq",
    "all_selected_taxa_retained_after_macro_eukaryote_exclusion"
  ),
  pass = c(
    !anyDuplicated(raw_metrics$Original_Taxon_ID),
    length(raw_candidate_ids) == EXPECTED_RAW_CANDIDATES,
    all(raw_candidate_ids %in% rownames(otu)),
    all(raw_candidate_ids %in% full_taxa),
    !anyDuplicated(membership$Original_Taxon_ID),
    length(candidate_ids) == EXPECTED_CANDIDATES,
    all(candidate_ids %in% rownames(otu)) &&
      all(membership$in_phyloseq %in% TRUE),
    all(candidate_ids %in% full_taxa)
  ),
  detail = c(
    paste0("n_metrics=", nrow(raw_metrics)),
    paste0("observed_raw=", length(raw_candidate_ids)),
    paste0(
      "mapped=",
      sum(raw_candidate_ids %in% rownames(otu)),
      "/",
      length(raw_candidate_ids)
    ),
    paste0(
      "retained=",
      sum(raw_candidate_ids %in% full_taxa),
      "/",
      length(raw_candidate_ids)
    ),
    paste0("n=", nrow(membership)),
    paste0("n=", length(candidate_ids)),
    paste0("mapped=", sum(candidate_ids %in% rownames(otu))),
    paste0("retained=", sum(candidate_ids %in% full_taxa))
  )
)
write_csv(
  universe_gate,
  file.path(TABLE_DIR, "301_compuerta_universos_taxonomicos.csv")
)
if (!all(universe_gate$pass)) {
  stop("Fallo la compuerta de universos taxonomicos de 301", call. = FALSE)
}

universe_overlap <- tibble(
  comparison = "raw1031_vs_residual84",
  n_intersection = length(intersect(raw_candidate_ids, candidate_ids)),
  n_union = length(union(raw_candidate_ids, candidate_ids))
) %>%
  mutate(jaccard = n_intersection / n_union)

write_csv(
  universe_overlap,
  file.path(TABLE_DIR, "301_solapamiento_universos_candidatos.csv")
)

counts_raw1031 <- counts_full[, raw_candidate_ids, drop = FALSE]
counts_candidate84 <- counts_full[, candidate_ids, drop = FALSE]

# ----------------------------------------------------------------------------
# Cobertura de los subconjuntos y poblacion analitica por perfiles completos
# ----------------------------------------------------------------------------

audit_subset_coverage <- function(counts, universe) {
  totals <- rowSums(counts)
  tibble(
    sample_id = rownames(counts),
    universe = universe,
    n_taxa_universe = ncol(counts),
    subset_library_size = as.numeric(totals),
    n_detected_taxa = as.integer(rowSums(counts > 0)),
    zero_sum = totals <= 0
  )
}

subset_coverage <- bind_rows(
  audit_subset_coverage(
    counts_full,
    "01_full_species_community"
  ),
  audit_subset_coverage(
    counts_raw1031,
    "02_raw_nonunimodal_candidates_1031"
  ),
  audit_subset_coverage(
    counts_candidate84,
    "02_residual_nonunimodal_candidates_84"
  )
) %>%
  left_join(
    meta %>%
      transmute(
        sample_id,
        profile_id,
        locality = as.character(locality),
        restoration4 = as.character(restoration4),
        depth_cm = as.numeric(as.character(depth_cm))
      ),
    by = "sample_id"
  ) %>%
  arrange(universe, locality, profile_id, depth_cm)

write_csv(
  subset_coverage,
  file.path(TABLE_DIR, "301_auditoria_cobertura_subconjuntos_por_muestra.csv")
)

zero_candidate_samples <- subset_coverage %>%
  filter(
    universe == "02_residual_nonunimodal_candidates_84",
    zero_sum
  ) %>%
  pull(sample_id)
zero_raw1031_samples <- subset_coverage %>%
  filter(
    universe == "02_raw_nonunimodal_candidates_1031",
    zero_sum
  ) %>%
  pull(sample_id)
zero_subset_samples <- zero_candidate_samples
excluded_subset_profiles <- meta %>%
  filter(sample_id %in% zero_subset_samples) %>%
  pull(profile_id) %>%
  unique()

zero_trigger <- subset_coverage %>%
  filter(
    universe == "02_residual_nonunimodal_candidates_84",
    zero_sum
  ) %>%
  group_by(sample_id) %>%
  summarise(
    zero_in_universes = paste(universe, collapse = "|"),
    .groups = "drop"
  )

excluded_profile_audit <- meta %>%
  filter(profile_id %in% excluded_subset_profiles) %>%
  transmute(
    profile_id,
    sample_id,
    locality = as.character(locality),
    restoration4 = as.character(restoration4),
    depth_cm = as.numeric(as.character(depth_cm))
  ) %>%
  left_join(zero_trigger, by = "sample_id") %>%
  mutate(
    triggers_profile_exclusion = !is.na(zero_in_universes),
    exclusion_reason =
      "complete_profile_removed_because_at_least_one_depth_has_zero_subset_sum"
  ) %>%
  arrange(locality, profile_id, depth_cm)

write_csv(
  excluded_profile_audit,
  file.path(TABLE_DIR, "301_auditoria_perfiles_excluidos_subconjuntos.csv")
)

meta_subset <- meta %>%
  filter(!profile_id %in% excluded_subset_profiles) %>%
  droplevels() %>%
  arrange(locality, profile_id, depth_cm)
counts_candidate84_subset <- counts_candidate84[
  meta_subset$sample_id,
  ,
  drop = FALSE
]

subset_profile_audit <- meta_subset %>%
  group_by(profile_id) %>%
  summarise(
    locality = first(as.character(locality)),
    restoration4 = first(as.character(restoration4)),
    n = n(),
    n_localities = n_distinct(locality),
    n_stages = n_distinct(restoration4),
    depths = paste(
      sort(as.numeric(as.character(depth_cm))),
      collapse = "|"
    ),
    valid = n == 3L &&
      n_localities == 1L &&
      n_stages == 1L &&
      setequal(as.numeric(as.character(depth_cm)), DEPTHS),
    .groups = "drop"
  )

write_csv(
  subset_profile_audit,
  file.path(TABLE_DIR, "301_auditoria_perfiles_subconjuntos.csv")
)

X_subset_primary <- model.matrix(
  ~ locality + depth_cm + restoration4,
  data = meta_subset
)
X_subset_eci <- model.matrix(
  ~ locality + depth_cm + ECI_PC1 + restoration4,
  data = meta_subset
)

subset_population_gate <- tibble(
  gate = c(
    "full_community_has_no_zero_samples",
    "candidate84_expected_three_zero_samples",
    "three_complete_profiles_excluded",
    "subset_population_has_42_samples",
    "subset_population_has_14_complete_profiles",
    "subset_population_retains_three_localities",
    "subset_population_retains_five_restoration_stages",
    "candidate84_positive_after_profile_exclusion",
    "subset_primary_model_full_rank",
    "subset_ECI_sensitivity_model_full_rank"
  ),
  pass = c(
    !any(subset_coverage$zero_sum[
      subset_coverage$universe == "01_full_species_community"
    ]),
    length(zero_candidate_samples) == 3L,
    length(excluded_subset_profiles) == 3L,
    nrow(meta_subset) == 42L,
    nrow(subset_profile_audit) == 14L && all(subset_profile_audit$valid),
    n_distinct(meta_subset$locality) == 3L,
    setequal(as.character(meta_subset$restoration4), STAGES),
    all(rowSums(counts_candidate84_subset) > 0),
    qr(X_subset_primary)$rank == ncol(X_subset_primary),
    qr(X_subset_eci)$rank == ncol(X_subset_eci)
  ),
  detail = c(
    paste0(
      "zero=",
      sum(subset_coverage$zero_sum[
        subset_coverage$universe == "01_full_species_community"
      ])
    ),
    paste(zero_candidate_samples, collapse = "|"),
    paste(excluded_subset_profiles, collapse = "|"),
    paste0("n=", nrow(meta_subset)),
    paste0(
      "profiles=", nrow(subset_profile_audit),
      "; valid=", sum(subset_profile_audit$valid)
    ),
    paste(levels(meta_subset$locality), collapse = "|"),
    paste(unique(as.character(meta_subset$restoration4)), collapse = "|"),
    paste0(
      "min_subset_library=",
      min(rowSums(counts_candidate84_subset))
    ),
    paste0(
      "rank=", qr(X_subset_primary)$rank,
      "/", ncol(X_subset_primary)
    ),
    paste0("rank=", qr(X_subset_eci)$rank, "/", ncol(X_subset_eci))
  )
)

write_csv(
  subset_population_gate,
  file.path(TABLE_DIR, "301_compuerta_poblacion_subconjuntos.csv")
)
if (!all(subset_population_gate$pass)) {
  stop(
    "Fallo la compuerta de poblacion analitica para los subconjuntos",
    call. = FALSE
  )
}

# Poblacion maxima positiva de 1,031 y poblacion comun para comparar 1,031 vs 84.
excluded_raw1031_profiles <- meta %>%
  filter(sample_id %in% zero_raw1031_samples) %>%
  pull(profile_id) %>%
  unique()

comparison_excluded_profiles <- union(
  excluded_subset_profiles,
  excluded_raw1031_profiles
)

meta_raw1031 <- meta %>%
  filter(!profile_id %in% excluded_raw1031_profiles) %>%
  droplevels() %>%
  arrange(locality, profile_id, depth_cm)

meta_comparison <- meta %>%
  filter(!profile_id %in% comparison_excluded_profiles) %>%
  droplevels() %>%
  arrange(locality, profile_id, depth_cm)

counts_raw1031_available <- counts_raw1031[
  meta_raw1031$sample_id,
  ,
  drop = FALSE
]
counts_full_comparison <- counts_full[
  meta_comparison$sample_id,
  ,
  drop = FALSE
]
counts_raw1031_comparison <- counts_raw1031[
  meta_comparison$sample_id,
  ,
  drop = FALSE
]
counts_candidate84_comparison <- counts_candidate84[
  meta_comparison$sample_id,
  ,
  drop = FALSE
]

raw1031_profile_audit <- meta_raw1031 %>%
  group_by(profile_id) %>%
  summarise(
    locality = first(as.character(locality)),
    restoration4 = first(as.character(restoration4)),
    n = n(),
    depths = paste(sort(as.numeric(as.character(depth_cm))), collapse = "|"),
    valid = n == 3L &&
      n_distinct(locality) == 1L &&
      n_distinct(restoration4) == 1L &&
      setequal(as.numeric(as.character(depth_cm)), DEPTHS),
    .groups = "drop"
  )

comparison_profile_audit <- meta_comparison %>%
  group_by(profile_id) %>%
  summarise(
    locality = first(as.character(locality)),
    restoration4 = first(as.character(restoration4)),
    n = n(),
    depths = paste(sort(as.numeric(as.character(depth_cm))), collapse = "|"),
    valid = n == 3L &&
      n_distinct(locality) == 1L &&
      n_distinct(restoration4) == 1L &&
      setequal(as.numeric(as.character(depth_cm)), DEPTHS),
    .groups = "drop"
  )

write_csv(
  raw1031_profile_audit,
  file.path(TABLE_DIR, "301_auditoria_perfiles_raw1031.csv")
)
write_csv(
  comparison_profile_audit,
  file.path(TABLE_DIR, "301_auditoria_perfiles_comparacion_1031_vs_84.csv")
)

raw1031_exclusion_audit <- meta %>%
  filter(profile_id %in% excluded_raw1031_profiles) %>%
  transmute(
    profile_id,
    sample_id,
    locality = as.character(locality),
    restoration4 = as.character(restoration4),
    depth_cm = as.numeric(as.character(depth_cm)),
    triggers_profile_exclusion = sample_id %in% zero_raw1031_samples,
    exclusion_reason =
      "complete_profile_removed_because_raw1031_has_zero_sum_at_one_or_more_depths"
  ) %>%
  arrange(locality, profile_id, depth_cm)

write_csv(
  raw1031_exclusion_audit,
  file.path(TABLE_DIR, "301_auditoria_perfiles_excluidos_raw1031.csv")
)

population_full_label <- "canonical_all_51"
population_raw1031_label <- paste0(
  "raw1031_positive_",
  nrow(meta_raw1031),
  "_samples_",
  n_distinct(meta_raw1031$profile_id),
  "_complete_profiles"
)
population_comparison_label <- paste0(
  "candidate_common_",
  nrow(meta_comparison),
  "_samples_",
  n_distinct(meta_comparison$profile_id),
  "_complete_profiles"
)

analysis_design_counts <- bind_rows(
  meta %>% mutate(analysis_population = population_full_label),
  meta_raw1031 %>%
    mutate(analysis_population = population_raw1031_label),
  meta_comparison %>%
    mutate(analysis_population = population_comparison_label)
) %>%
  group_by(analysis_population, locality, restoration4) %>%
  summarise(
    n_samples = n(),
    n_profiles = n_distinct(profile_id),
    .groups = "drop"
  ) %>%
  arrange(analysis_population, locality, restoration4)

write_csv(
  analysis_design_counts,
  file.path(TABLE_DIR, "301_diseno_analitico_por_localidad_etapa.csv")
)

X_raw1031_primary <- model.matrix(
  ~ locality + depth_cm + restoration4,
  data = meta_raw1031
)
X_raw1031_eci <- model.matrix(
  ~ locality + depth_cm + ECI_PC1 + restoration4,
  data = meta_raw1031
)
X_comparison_primary <- model.matrix(
  ~ locality + depth_cm + restoration4,
  data = meta_comparison
)
X_comparison_eci <- model.matrix(
  ~ locality + depth_cm + ECI_PC1 + restoration4,
  data = meta_comparison
)

candidate_comparison_gate <- tibble(
  gate = c(
    "raw1031_positive_after_profile_exclusion",
    "raw1031_profiles_complete",
    "raw1031_retains_three_localities",
    "raw1031_retains_five_restoration_stages",
    "raw1031_primary_model_full_rank",
    "raw1031_ECI_model_full_rank",
    "comparison_profiles_complete",
    "comparison_retains_three_localities",
    "comparison_retains_five_restoration_stages",
    "full_community_positive_in_comparison_population",
    "raw1031_positive_in_comparison_population",
    "candidate84_positive_in_comparison_population",
    "comparison_primary_model_full_rank",
    "comparison_ECI_model_full_rank"
  ),
  pass = c(
    all(rowSums(counts_raw1031_available) > 0),
    nrow(raw1031_profile_audit) > 0L && all(raw1031_profile_audit$valid),
    n_distinct(meta_raw1031$locality) == 3L,
    setequal(as.character(meta_raw1031$restoration4), STAGES),
    qr(X_raw1031_primary)$rank == ncol(X_raw1031_primary),
    qr(X_raw1031_eci)$rank == ncol(X_raw1031_eci),
    nrow(comparison_profile_audit) > 0L &&
      all(comparison_profile_audit$valid),
    n_distinct(meta_comparison$locality) == 3L,
    setequal(as.character(meta_comparison$restoration4), STAGES),
    all(rowSums(counts_full_comparison) > 0),
    all(rowSums(counts_raw1031_comparison) > 0),
    all(rowSums(counts_candidate84_comparison) > 0),
    qr(X_comparison_primary)$rank == ncol(X_comparison_primary),
    qr(X_comparison_eci)$rank == ncol(X_comparison_eci)
  ),
  detail = c(
    paste0("min_library=", min(rowSums(counts_raw1031_available))),
    paste0(
      "profiles=", nrow(raw1031_profile_audit),
      "; valid=", sum(raw1031_profile_audit$valid)
    ),
    paste(levels(meta_raw1031$locality), collapse = "|"),
    paste(unique(as.character(meta_raw1031$restoration4)), collapse = "|"),
    paste0(
      "rank=", qr(X_raw1031_primary)$rank,
      "/", ncol(X_raw1031_primary)
    ),
    paste0("rank=", qr(X_raw1031_eci)$rank, "/", ncol(X_raw1031_eci)),
    paste0(
      "profiles=", nrow(comparison_profile_audit),
      "; valid=", sum(comparison_profile_audit$valid)
    ),
    paste(levels(meta_comparison$locality), collapse = "|"),
    paste(unique(as.character(meta_comparison$restoration4)), collapse = "|"),
    paste0("min_library=", min(rowSums(counts_full_comparison))),
    paste0("min_library=", min(rowSums(counts_raw1031_comparison))),
    paste0("min_library=", min(rowSums(counts_candidate84_comparison))),
    paste0(
      "rank=", qr(X_comparison_primary)$rank,
      "/", ncol(X_comparison_primary)
    ),
    paste0(
      "rank=", qr(X_comparison_eci)$rank,
      "/", ncol(X_comparison_eci)
    )
  )
)

write_csv(
  candidate_comparison_gate,
  file.path(TABLE_DIR, "301_compuerta_comparacion_1031_vs_84.csv")
)
if (!all(candidate_comparison_gate$pass)) {
  stop(
    "Fallo la compuerta de poblacion comun para comparar 1031 frente a 84",
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# Permutaciones auditadas
# -----------------------------------------------------------------------------

permutations_rows <- profile_permutations_rows(meta, N_PERM, SEED)
row_permutation_audit <- audit_row_permutations(permutations_rows, meta)
write_csv(
  row_permutation_audit,
  file.path(TABLE_DIR, "301_auditoria_permutaciones_splitplot.csv")
)

permutations_rows_subset <- profile_permutations_rows(
  meta_subset,
  N_PERM,
  SEED + 40L
)
row_permutation_audit_subset <- audit_row_permutations(
  permutations_rows_subset,
  meta_subset
)
write_csv(
  row_permutation_audit_subset,
  file.path(
    TABLE_DIR,
    "301_auditoria_permutaciones_splitplot_subconjuntos.csv"
  )
)

meta_profile_key <- meta %>%
  distinct(profile_id, locality, restoration4) %>%
  arrange(locality, profile_id)
permutations_profiles <- profile_permutations_profiles(
  meta_profile_key,
  N_PERM,
  SEED + 20L
)
profile_permutation_audit <- audit_profile_permutations(
  permutations_profiles,
  meta_profile_key
)
write_csv(
  profile_permutation_audit,
  file.path(TABLE_DIR, "301_auditoria_permutaciones_perfiles.csv")
)

permutation_gate <- tibble(
  gate = c(
    "row_permutations_generated",
    "row_permutations_preserve_whole_profiles",
    "subset_row_permutations_generated",
    "subset_row_permutations_preserve_whole_profiles",
    "profile_permutations_generated",
    "profile_permutations_stay_within_locality"
  ),
  pass = c(
    nrow(permutations_rows) > 0L,
    all(row_permutation_audit$valid_whole_profile_within_locality),
    nrow(permutations_rows_subset) > 0L,
    all(
      row_permutation_audit_subset$valid_whole_profile_within_locality
    ),
    nrow(permutations_profiles) > 0L,
    all(profile_permutation_audit$valid_profile_within_locality)
  ),
  detail = c(
    paste0("unique=", nrow(permutations_rows)),
    paste0(
      "valid=",
      sum(row_permutation_audit$valid_whole_profile_within_locality),
      "/",
      nrow(row_permutation_audit)
    ),
    paste0("unique=", nrow(permutations_rows_subset)),
    paste0(
      "valid=",
      sum(
        row_permutation_audit_subset$valid_whole_profile_within_locality
      ),
      "/",
      nrow(row_permutation_audit_subset)
    ),
    paste0("unique=", nrow(permutations_profiles)),
    paste0(
      "valid=",
      sum(profile_permutation_audit$valid_profile_within_locality),
      "/",
      nrow(profile_permutation_audit)
    )
  )
)
write_csv(
  permutation_gate,
  file.path(TABLE_DIR, "301_compuerta_permutaciones.csv")
)
if (!all(permutation_gate$pass)) {
  stop("Fallo la compuerta de permutaciones de 301", call. = FALSE)
}

permutations_rows_raw1031 <- if (identical(
  meta_raw1031$sample_id,
  meta$sample_id
)) {
  permutations_rows
} else {
  profile_permutations_rows(meta_raw1031, N_PERM, SEED + 30L)
}

row_permutation_audit_raw1031 <- audit_row_permutations(
  permutations_rows_raw1031,
  meta_raw1031
)
write_csv(
  row_permutation_audit_raw1031,
  file.path(TABLE_DIR, "301_auditoria_permutaciones_splitplot_raw1031.csv")
)

permutations_rows_comparison <- if (identical(
  meta_comparison$sample_id,
  meta_subset$sample_id
)) {
  permutations_rows_subset
} else if (identical(meta_comparison$sample_id, meta$sample_id)) {
  permutations_rows
} else {
  profile_permutations_rows(meta_comparison, N_PERM, SEED + 50L)
}

row_permutation_audit_comparison <- audit_row_permutations(
  permutations_rows_comparison,
  meta_comparison
)
write_csv(
  row_permutation_audit_comparison,
  file.path(
    TABLE_DIR,
    "301_auditoria_permutaciones_splitplot_comparacion_1031_vs_84.csv"
  )
)

comparison_permutation_gate <- tibble(
  gate = c(
    "raw1031_permutations_generated",
    "raw1031_permutations_preserve_whole_profiles",
    "comparison_permutations_generated",
    "comparison_permutations_preserve_whole_profiles"
  ),
  pass = c(
    nrow(permutations_rows_raw1031) > 0L,
    all(
      row_permutation_audit_raw1031$valid_whole_profile_within_locality
    ),
    nrow(permutations_rows_comparison) > 0L,
    all(
      row_permutation_audit_comparison$valid_whole_profile_within_locality
    )
  ),
  detail = c(
    paste0("unique=", nrow(permutations_rows_raw1031)),
    paste0(
      "valid=",
      sum(
        row_permutation_audit_raw1031$valid_whole_profile_within_locality
      ),
      "/",
      nrow(row_permutation_audit_raw1031)
    ),
    paste0("unique=", nrow(permutations_rows_comparison)),
    paste0(
      "valid=",
      sum(
        row_permutation_audit_comparison$valid_whole_profile_within_locality
      ),
      "/",
      nrow(row_permutation_audit_comparison)
    )
  )
)

write_csv(
  comparison_permutation_gate,
  file.path(TABLE_DIR, "301_compuerta_permutaciones_comparacion_1031_vs_84.csv")
)
if (!all(comparison_permutation_gate$pass)) {
  stop(
    "Fallo la compuerta de permutaciones de la comparacion 1031 frente a 84",
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 1. Comunidad completa: diversidad Hill
# -----------------------------------------------------------------------------

cat("Stage 1/3: full-community Hill diversity and sensitivities\n")
flush(log_con)

alpha_sample_observed <- hill_numbers(counts_full) %>%
  left_join(meta, by = "sample_id")

if (!"library_size" %in% names(alpha_sample_observed) ||
    any(!is.finite(alpha_sample_observed$library_size)) ||
    any(alpha_sample_observed$library_size <= 0)) {
  stop(
    "La tabla Hill no conserva library_size positiva y finita",
    call. = FALSE
  )
}

rarefaction_depth <- min(alpha_sample_observed$library_size)
alpha_sample_rarefied <- mean_rarefied_hill(
  counts_full,
  rarefaction_depth,
  N_RAREFY,
  SEED + 10L
) %>%
  left_join(meta, by = "sample_id")

alpha_profile_observed <- alpha_sample_observed %>%
  group_by(profile_id, locality, restoration4) %>%
  summarise(
    q0_richness = mean(q0_richness),
    q1_exp_shannon = mean(q1_exp_shannon),
    q2_inverse_simpson = mean(q2_inverse_simpson),
    log10_library_size = mean(log10(library_size)),
    n_depths = n_distinct(depth_cm),
    .groups = "drop"
  ) %>%
  arrange(locality, profile_id)

alpha_profile_rarefied <- alpha_sample_rarefied %>%
  group_by(profile_id, locality, restoration4) %>%
  summarise(
    q0_richness = mean(q0_richness),
    q1_exp_shannon = mean(q1_exp_shannon),
    q2_inverse_simpson = mean(q2_inverse_simpson),
    rarefaction_depth = first(rarefaction_depth),
    rarefaction_iterations = first(rarefaction_iterations),
    n_depths = n_distinct(depth_cm),
    .groups = "drop"
  ) %>%
  arrange(locality, profile_id)

if (!identical(alpha_profile_observed$profile_id, meta_profile_key$profile_id)) {
  stop("Orden de perfiles inconsistente en diversidad alfa", call. = FALSE)
}
if (!identical(alpha_profile_rarefied$profile_id, meta_profile_key$profile_id)) {
  stop("Orden de perfiles rarefaccion inconsistente", call. = FALSE)
}

hill_metrics <- c(
  "q0_richness",
  "q1_exp_shannon",
  "q2_inverse_simpson"
)

alpha_tests_primary <- bind_rows(lapply(
  hill_metrics,
  function(metric) {
    profile_alpha_test(
      alpha_profile_observed,
      metric,
      permutations_profiles,
      "primary_observed_Hill_profile_mean",
      FALSE,
      SEED + 20L
    )
  }
))

alpha_tests_rarefied <- bind_rows(lapply(
  hill_metrics,
  function(metric) {
    profile_alpha_test(
      alpha_profile_rarefied,
      metric,
      permutations_profiles,
      "rarefied_depth_sensitivity",
      FALSE,
      SEED + 20L
    )
  }
))

alpha_tests_library_adjusted <- bind_rows(lapply(
  hill_metrics,
  function(metric) {
    profile_alpha_test(
      alpha_profile_observed,
      metric,
      permutations_profiles,
      "library_size_adjusted_sensitivity",
      TRUE,
      SEED + 20L
    )
  }
))

alpha_tests <- bind_rows(
  alpha_tests_primary,
  alpha_tests_rarefied,
  alpha_tests_library_adjusted
) %>%
  group_by(analysis_role) %>%
  mutate(q_value_BH_within_role = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    manuscript_order = case_when(
      analysis_role == "primary_observed_Hill_profile_mean" ~ 1L,
      analysis_role == "rarefied_depth_sensitivity" ~ 2L,
      TRUE ~ 3L
    )
  ) %>%
  arrange(manuscript_order, match(response, hill_metrics))

write_csv(
  alpha_sample_observed,
  file.path(TABLE_DIR, "301_diversidad_Hill_completa_por_muestra.csv")
)
write_csv(
  alpha_sample_rarefied,
  file.path(TABLE_DIR, "301_diversidad_Hill_rarefaccion_por_muestra.csv")
)
write_csv(
  alpha_profile_observed,
  file.path(TABLE_DIR, "301_diversidad_Hill_completa_por_perfil.csv")
)
write_csv(
  alpha_profile_rarefied,
  file.path(TABLE_DIR, "301_diversidad_Hill_rarefaccion_por_perfil.csv")
)
write_csv(
  alpha_tests,
  file.path(TABLE_DIR, "301_pruebas_Hill_comunidad_completa.csv")
)

rarefaction_audit <- tibble(
  parameter = c(
    "minimum_library_size",
    "median_library_size",
    "maximum_library_size",
    "rarefaction_depth",
    "rarefaction_iterations",
    "minimum_over_median"
  ),
  value = c(
    min(alpha_sample_observed$library_size),
    median(alpha_sample_observed$library_size),
    max(alpha_sample_observed$library_size),
    rarefaction_depth,
    N_RAREFY,
    rarefaction_depth / median(alpha_sample_observed$library_size)
  )
)
write_csv(
  rarefaction_audit,
  file.path(TABLE_DIR, "301_auditoria_profundidad_secuenciacion.csv")
)

cat("Stage 1/3 complete: full-community Hill diversity\n")
flush(log_con)

# -----------------------------------------------------------------------------
# 2. Ordenamiento y PERMANOVA: comunidad completa y sensibilidades candidatas
# -----------------------------------------------------------------------------

cat("Stage 2/3: community composition and candidate sensitivities\n")
flush(log_con)

community_specs <- list(
  list(
    community = "01_full_species_community_all51",
    order = 1L,
    counts = counts_full,
    metadata = meta %>%
      mutate(analysis_population = population_full_label),
    permutations = permutations_rows
  ),
  list(
    community = "02_full_species_community_matched_candidates",
    order = 2L,
    counts = counts_full_comparison,
    metadata = meta_comparison %>%
      mutate(analysis_population = population_comparison_label),
    permutations = permutations_rows_comparison
  ),
  list(
    community = "03_raw_nonunimodal_candidates_1031_available",
    order = 3L,
    counts = counts_raw1031_available,
    metadata = meta_raw1031 %>%
      mutate(analysis_population = population_raw1031_label),
    permutations = permutations_rows_raw1031
  ),
  list(
    community = "04_raw_nonunimodal_candidates_1031_matched",
    order = 4L,
    counts = counts_raw1031_comparison,
    metadata = meta_comparison %>%
      mutate(analysis_population = population_comparison_label),
    permutations = permutations_rows_comparison
  ),
  list(
    community = "05_residual_nonunimodal_candidates_84_matched",
    order = 5L,
    counts = counts_candidate84_comparison,
    metadata = meta_comparison %>%
      mutate(analysis_population = population_comparison_label),
    permutations = permutations_rows_comparison
  )
)

count_community_results <- lapply(
  community_specs,
  function(spec) {
    run_community(
      counts = spec$counts,
      meta = spec$metadata,
      permutations = spec$permutations,
      community = spec$community,
      community_order = spec$order,
      pseudocount = PSEUDOCOUNT,
      seed = SEED
    )
  }
)

# Un solo CLR completo para los tres objetos comparables. Al seleccionar
# columnas despues del centrado, 1,031 y 84 se prueban en las mismas 51
# muestras, sin redefinir el cierre composicional ni excluir perfiles.
full_clr_all51 <- clr_rows(counts_full, PSEUDOCOUNT)
distance_full_aitchison_all51 <- stats::dist(full_clr_all51)
distance_raw1031_full_clr_all51 <- stats::dist(
  full_clr_all51[, raw_candidate_ids, drop = FALSE]
)
distance_residual84_full_clr_all51 <- stats::dist(
  full_clr_all51[, candidate_ids, drop = FALSE]
)

common_clr_gate <- tibble(
  gate = c(
    "full_CLR_is_51_by_12761",
    "raw1031_coordinates_selected_after_full_CLR",
    "residual84_coordinates_selected_after_full_CLR",
    "common_CLR_distances_are_finite",
    "common_CLR_distances_share_all51_sample_order"
  ),
  pass = c(
    nrow(full_clr_all51) == 51L && ncol(full_clr_all51) == length(full_taxa),
    all(raw_candidate_ids %in% colnames(full_clr_all51)),
    all(candidate_ids %in% colnames(full_clr_all51)),
    all(is.finite(as.vector(distance_raw1031_full_clr_all51))) &&
      all(is.finite(as.vector(distance_residual84_full_clr_all51))),
    identical(attr(distance_raw1031_full_clr_all51, "Labels"), meta$sample_id) &&
      identical(
        attr(distance_residual84_full_clr_all51, "Labels"),
        meta$sample_id
      )
  ),
  detail = c(
    paste(dim(full_clr_all51), collapse = "x"),
    paste0("n_coordinates=", length(raw_candidate_ids)),
    paste0("n_coordinates=", length(candidate_ids)),
    "finite Euclidean distances on selected full-CLR coordinates",
    "same 51 samples in canonical order"
  )
)
write_csv(
  common_clr_gate,
  file.path(TABLE_DIR, "301_compuerta_CLR_comun_coordenadas_candidatas.csv")
)
if (!all(common_clr_gate$pass)) {
  stop("Fallo la compuerta del CLR comun", call. = FALSE)
}

meta_all51_analysis <- meta %>%
  mutate(analysis_population = population_full_label)

common_clr_results <- list(
  run_precomputed_distance(
    distance_raw1031_full_clr_all51,
    meta_all51_analysis,
    permutations_rows,
    "06_raw1031_selected_from_full_CLR_all51",
    6L,
    "Euclidean_selected_full_CLR_sensitivity",
    SEED + 60L
  ),
  run_precomputed_distance(
    distance_residual84_full_clr_all51,
    meta_all51_analysis,
    permutations_rows,
    "07_residual84_selected_from_full_CLR_all51",
    7L,
    "Euclidean_selected_full_CLR_sensitivity",
    SEED + 70L
  )
)

community_results <- c(count_community_results, common_clr_results)

permanova <- bind_rows(lapply(community_results, `[[`, "permanova")) %>%
  arrange(community_order, model_role, metric, term)
stage_rows <- permanova$term == "restoration4"
permanova$p_BH_across_all_stage_tests <- NA_real_
permanova$p_BH_across_all_stage_tests[stage_rows] <- p.adjust(
  permanova$p_value[stage_rows],
  method = "BH"
)

# Comparacion dirigida entre los dos universos candidatos sobre coordenadas
# seleccionadas de un unico CLR completo. Ambos usan exactamente las mismas 51
# muestras, los mismos 17 perfiles y las mismas permutaciones.
target_candidate_rows <- with(
  permanova,
  term == "restoration4" &
    model_role == "primary_adjusted_locality_depth" &
    metric == "Euclidean_selected_full_CLR_sensitivity" &
    community %in% c(
      "06_raw1031_selected_from_full_CLR_all51",
      "07_residual84_selected_from_full_CLR_all51"
    )
)
permanova$p_BH_candidate_universe_primary_common_CLR <- NA_real_
permanova$p_BH_candidate_universe_primary_common_CLR[target_candidate_rows] <-
  p.adjust(permanova$p_value[target_candidate_rows], method = "BH")

dispersion <- bind_rows(lapply(community_results, `[[`, "dispersion")) %>%
  arrange(community_order, metric) %>%
  mutate(p_BH_across_dispersion_tests = p.adjust(p_value, method = "BH"))

scores <- bind_rows(lapply(community_results, `[[`, "scores")) %>%
  arrange(community_order, metric, locality, profile_id, depth_cm)

stage_results <- permanova %>%
  filter(term == "restoration4") %>%
  select(
    community_order,
    community,
    metric,
    model_role,
    inference_role,
    analysis_population,
    n_samples_analysis,
    n_profiles_analysis,
    df,
    sum_squares,
    R2,
    pseudo_F,
    p_value,
    p_BH_across_all_stage_tests,
    p_BH_candidate_universe_primary_common_CLR,
    n_unique_permutations,
    minimum_attainable_p,
    model,
    permutation_scheme
  ) %>%
  arrange(community_order, model_role, metric)

write_csv(
  permanova,
  file.path(TABLE_DIR, "301_PERMANOVA_todos_terminos.csv")
)
write_csv(
  stage_results,
  file.path(TABLE_DIR, "301_resultados_restauracion_orden_manuscrito.csv")
)
write_csv(
  stage_results %>% filter(community_order == 1L),
  file.path(TABLE_DIR, "301_01_resultado_comunidad_completa_51.csv")
)
write_csv(
  stage_results %>% filter(community_order == 2L),
  file.path(TABLE_DIR, "301_02_control_comunidad_completa_poblacion_comun.csv")
)
write_csv(
  stage_results %>% filter(community_order == 3L),
  file.path(TABLE_DIR, "301_03_resultado_candidatos1031_poblacion_disponible.csv")
)
write_csv(
  stage_results %>% filter(community_order == 4L),
  file.path(TABLE_DIR, "301_04_resultado_candidatos1031_poblacion_comun.csv")
)
write_csv(
  stage_results %>% filter(community_order == 5L),
  file.path(TABLE_DIR, "301_05_resultado_candidatos84_poblacion_comun.csv")
)
write_csv(
  stage_results %>% filter(community_order == 6L),
  file.path(TABLE_DIR, "301_06_resultado_1031_coordenadas_CLR_completo_51.csv")
)
write_csv(
  stage_results %>% filter(community_order == 7L),
  file.path(TABLE_DIR, "301_07_resultado_84_coordenadas_CLR_completo_51.csv")
)
write_csv(
  stage_results %>%
    filter(
      community_order %in% c(6L, 7L),
      metric == "Euclidean_selected_full_CLR_sensitivity",
      model_role == "primary_adjusted_locality_depth"
    ),
  file.path(TABLE_DIR, "301_comparacion_dirigida_1031_vs_84_CLR_comun_51.csv")
)
write_csv(
  dispersion,
  file.path(TABLE_DIR, "301_dispersion_restauracion.csv")
)
write_csv(
  scores,
  file.path(TABLE_DIR, "301_PCoA_source_data.csv")
)

cat("Stage 2/3 complete: community composition\n")
flush(log_con)

# -----------------------------------------------------------------------------
# 3. Sensibilidad de influencia: dejar fuera un perfil completo cada vez
# -----------------------------------------------------------------------------

cat("Stage 3/3: leave-one-complete-profile-out influence sensitivity\n")
flush(log_con)

lopo_results <- bind_rows(
  run_lopo_stage(
    distance_full_aitchison_all51,
    meta_all51_analysis,
    "01_full_species_community_all51",
    1L,
    N_PERM,
    SEED + 1000L
  ),
  run_lopo_stage(
    distance_raw1031_full_clr_all51,
    meta_all51_analysis,
    "06_raw1031_selected_from_full_CLR_all51",
    6L,
    N_PERM,
    SEED + 2000L
  ),
  run_lopo_stage(
    distance_residual84_full_clr_all51,
    meta_all51_analysis,
    "07_residual84_selected_from_full_CLR_all51",
    7L,
    N_PERM,
    SEED + 3000L
  )
)

lopo_baselines <- stage_results %>%
  filter(
    model_role == "primary_adjusted_locality_depth",
    (community == "01_full_species_community_all51" &
      metric == "Aitchison_primary") |
      (community %in% c(
        "06_raw1031_selected_from_full_CLR_all51",
        "07_residual84_selected_from_full_CLR_all51"
      ) & metric == "Euclidean_selected_full_CLR_sensitivity")
  ) %>%
  transmute(
    analysis = community,
    baseline_R2 = R2,
    baseline_p = p_value
  )

lopo_results <- lopo_results %>%
  left_join(lopo_baselines, by = "analysis") %>%
  mutate(
    delta_R2_from_all17 = R2 - baseline_R2,
    nominal_alpha_class_changed = (p_value < 0.05) != (baseline_p < 0.05)
  ) %>%
  arrange(analysis, omitted_locality, omitted_profile)

lopo_summary <- lopo_results %>%
  group_by(analysis) %>%
  summarise(
    baseline_R2 = first(baseline_R2),
    baseline_p = first(baseline_p),
    n_profile_omissions = n(),
    R2_min = min(R2),
    R2_median = median(R2),
    R2_max = max(R2),
    max_abs_delta_R2 = max(abs(delta_R2_from_all17)),
    p_min = min(p_value),
    p_median = median(p_value),
    p_max = max(p_value),
    n_nominal_p_below_0_05 = sum(p_value < 0.05),
    n_nominal_alpha_class_changes = sum(nominal_alpha_class_changed),
    all_permutations_preserve_whole_profiles =
      all(all_permutations_preserve_whole_profiles),
    .groups = "drop"
  )

write_csv(
  lopo_results,
  file.path(TABLE_DIR, "301_LOPO_resultados_por_perfil_omitido.csv")
)
write_csv(
  lopo_summary,
  file.path(TABLE_DIR, "301_LOPO_resumen_influencia.csv")
)

lopo_gate <- tibble(
  gate = c(
    "LOPO_has_three_analyses_times_17_profiles",
    "LOPO_omits_whole_profiles",
    "LOPO_retains_16_complete_profiles_each",
    "LOPO_statistics_are_finite",
    "LOPO_permutations_preserve_whole_profiles"
  ),
  pass = c(
    nrow(lopo_results) == 3L * 17L,
    all(lopo_results$n_samples == 48L),
    all(lopo_results$n_profiles == 16L),
    all(is.finite(lopo_results$R2)) &&
      all(is.finite(lopo_results$pseudo_F)) &&
      all(is.finite(lopo_results$p_value)),
    all(lopo_results$all_permutations_preserve_whole_profiles)
  ),
  detail = c(
    paste0("n=", nrow(lopo_results)),
    "three depths removed together",
    paste(unique(lopo_results$n_profiles), collapse = "|"),
    "R2, pseudo-F and P finite",
    "audited within locality with depth order preserved"
  )
)
write_csv(lopo_gate, file.path(TABLE_DIR, "301_compuerta_LOPO.csv"))
if (!all(lopo_gate$pass)) {
  stop("Fallo la compuerta LOPO", call. = FALSE)
}

cat("Stage 3/3 complete: LOPO influence sensitivity\n")
flush(log_con)

# -----------------------------------------------------------------------------
# Figuras fuente del manuscrito: 301 produce solo Fig. 1B y Fig. 1C.
# Los universos de 1,031 y 84 taxa se exportan como sensibilidades
# suplementarias; no se insertan automaticamente en la Figura 1.
# -----------------------------------------------------------------------------

theme_manuscript <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.title = element_text(color = "black", face = "plain"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = grid::unit(1.4, "mm"),
      panel.grid.major.y = element_line(
        color = "#E7E7E7",
        linewidth = 0.25
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(
        color = "black",
        face = "bold",
        size = rel(1.02),
        margin = margin(b = 2)
      ),
      plot.title = element_text(face = "bold", hjust = 0, size = rel(1.02)),
      plot.title.position = "plot",
      plot.subtitle = element_text(color = "#333333"),
      plot.caption = element_text(
        color = "#4D4D4D",
        size = rel(0.78),
        hjust = 0
      ),
      plot.tag = element_text(face = "bold", size = rel(1.35)),
      legend.key = element_blank(),
      legend.key.height = grid::unit(4.2, "mm"),
      legend.key.width = grid::unit(4.2, "mm"),
      legend.background = element_blank(),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = rel(0.92)),
      plot.margin = margin(8, 10, 8, 10)
    )
}

format_p_figure <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# Medias marginales ajustadas e IC95 calculados con el mismo modelo primario
# usado para cada numero de Hill: respuesta ~ localidad + etapa.
profile_hill_stage_ci <- function(data, response) {
  analysis_data <- data %>%
    mutate(
      locality = factor(locality),
      restoration4 = factor(restoration4, levels = STAGES)
    )
  fit <- lm(
    reformulate(c("locality", "restoration4"), response = response),
    data = analysis_data
  )
  if (anyNA(coef(fit))) {
    stop("Modelo Hill con coeficientes no estimables: ", response, call. = FALSE)
  }

  prediction_grid <- expand.grid(
    locality = levels(analysis_data$locality),
    restoration4 = STAGES,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      locality = factor(locality, levels = levels(analysis_data$locality)),
      restoration4 = factor(restoration4, levels = STAGES)
    )
  design <- model.matrix(
    delete.response(terms(fit)),
    prediction_grid,
    contrasts.arg = fit$contrasts
  )
  critical <- qt(0.975, df = df.residual(fit))

  bind_rows(lapply(STAGES, function(stage) {
    idx <- prediction_grid$restoration4 == stage
    contrast <- colMeans(design[idx, , drop = FALSE])
    estimate <- as.numeric(contrast %*% coef(fit))
    standard_error <- sqrt(
      as.numeric(t(contrast) %*% vcov(fit) %*% contrast)
    )
    tibble(
      response = response,
      restoration4 = factor(stage, levels = STAGES),
      estimate = estimate,
      standard_error = standard_error,
      conf_low = estimate - critical * standard_error,
      conf_high = estimate + critical * standard_error,
      confidence_level = 0.95,
      model = paste0(response, " ~ locality + restoration4"),
      estimand = "equal-weight marginal mean across localities"
    )
  }))
}

hill_titles <- c(
  "q0_richness" = "q0  Richness",
  "q1_exp_shannon" = "q1  exp(Shannon)",
  "q2_inverse_simpson" = "q2  Inverse Simpson"
)

hill_ci <- bind_rows(lapply(
  hill_metrics,
  function(metric) profile_hill_stage_ci(alpha_profile_observed, metric)
)) %>%
  mutate(
    hill_label = factor(
      unname(hill_titles[response]),
      levels = unname(hill_titles[hill_metrics])
    )
  )

write_csv(
  hill_ci,
  file.path(TABLE_DIR, "301_Hill_estimaciones_ajustadas_IC95_por_etapa.csv")
)

alpha_long <- alpha_profile_observed %>%
  pivot_longer(
    all_of(hill_metrics),
    names_to = "response",
    values_to = "value"
  ) %>%
  mutate(
    hill_label = factor(
      unname(hill_titles[response]),
      levels = unname(hill_titles[hill_metrics])
    )
  )

alpha_annotations <- alpha_tests %>%
  filter(analysis_role == "primary_observed_Hill_profile_mean") %>%
  transmute(
    response,
    hill_label = factor(
      unname(hill_titles[response]),
      levels = unname(hill_titles[hill_metrics])
    ),
    label = paste0(
      "Stage: F = ", sprintf("%.2f", observed_F),
      "; P = ", vapply(p_value, format_p_figure, character(1)),
      "; BH q = ",
      vapply(q_value_BH_within_role, format_p_figure, character(1))
    )
  )

profile_n <- alpha_profile_observed %>%
  count(restoration4, name = "n_profiles") %>%
  complete(
    restoration4 = factor(STAGES, levels = STAGES),
    fill = list(n_profiles = 0L)
  ) %>%
  arrange(restoration4)
stage_axis_labels <- setNames(
  paste0(
    unname(STAGE_LABELS_SHORT[STAGES]),
    "\n(n=", profile_n$n_profiles, ")"
  ),
  STAGES
)

p_alpha <- ggplot(alpha_long, aes(restoration4, value)) +
  geom_point(
    aes(fill = restoration4),
    shape = 21,
    color = "black",
    stroke = 0.3,
    size = 2.05,
    alpha = 0.84,
    position = position_jitter(width = 0.105, height = 0),
    show.legend = FALSE
  ) +
  geom_linerange(
    data = hill_ci,
    aes(
      x = restoration4,
      ymin = conf_low,
      ymax = conf_high
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.55
  ) +
  geom_point(
    data = hill_ci,
    aes(x = restoration4, y = estimate, fill = restoration4),
    inherit.aes = FALSE,
    shape = 23,
    color = "black",
    stroke = 0.55,
    size = 2.8,
    show.legend = FALSE
  ) +
  geom_text(
    data = alpha_annotations,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.02,
    vjust = 1.3,
    size = 2.05,
    color = "#2B2B2B"
  ) +
  facet_wrap(~ hill_label, ncol = 1L, scales = "free_y") +
  scale_fill_manual(values = STAGE_COLORS, limits = STAGES, drop = FALSE) +
  scale_x_discrete(labels = stage_axis_labels, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.07, 0.22))) +
  labs(
    tag = "B",
    x = "Restoration stage",
    y = "Profile-level Hill number"
  ) +
  theme_manuscript(base_size = 8.4) +
  theme(
    panel.spacing.y = grid::unit(2.2, "mm"),
    axis.text.x = element_text(size = 7.1),
    plot.tag.position = c(0.01, 0.99)
  )

plot_pcoa <- function(
  score_table,
  selected_community,
  selected_metric,
  panel_tag = NULL,
  title = NULL,
  show_population_caption = TRUE
) {
  plot_data <- score_table %>%
    filter(
      community == selected_community,
      metric == selected_metric
    ) %>%
    mutate(
      restoration4 = factor(restoration4, levels = STAGES),
      locality = factor(locality),
      depth_cm = factor(depth_cm, levels = as.character(DEPTHS))
    )
  if (!nrow(plot_data)) {
    stop("No hay scores para: ", selected_community, call. = FALSE)
  }

  x_variance <- unique(plot_data$PCoA1_variance_positive)
  y_variance <- unique(plot_data$PCoA2_variance_positive)
  x_label <- if (length(x_variance) == 1L && is.finite(x_variance)) {
    paste0("PCoA1 (", sprintf("%.1f", 100 * x_variance), "%)")
  } else {
    "PCoA1"
  }
  y_label <- if (length(y_variance) == 1L && is.finite(y_variance)) {
    paste0("PCoA2 (", sprintf("%.1f", 100 * y_variance), "%)")
  } else {
    "PCoA2"
  }
  population_label <- unique(plot_data$analysis_population)
  if (length(population_label) != 1L) {
    stop("La PCoA mezcla poblaciones analiticas", call. = FALSE)
  }

  stage_test <- stage_results %>%
    filter(
      community == selected_community,
      metric == selected_metric,
      model_role == "primary_adjusted_locality_depth"
    )
  dispersion_test <- dispersion %>%
    filter(
      community == selected_community,
      metric == selected_metric
    )
  if (nrow(stage_test) != 1L || nrow(dispersion_test) != 1L) {
    stop("Resultado inferencial ambiguo para la PCoA", call. = FALSE)
  }
  statistics_label <- paste0(
    "Stage PERMANOVA\nR\u00b2 = ", sprintf("%.3f", stage_test$R2),
    "; pseudo-F = ", sprintf("%.2f", stage_test$pseudo_F),
    "; P = ", format_p_figure(stage_test$p_value),
    "\nStage dispersion P = ",
    format_p_figure(dispersion_test$p_value)
  )
  caption_label <- if (show_population_caption) {
    paste0(
      n_distinct(plot_data$sample_id), " samples; ",
      n_distinct(plot_data$profile_id), " complete vertical profiles. ",
      "Ellipses show 68% normal-data regions and are descriptive."
    )
  } else {
    NULL
  }

  ggplot(plot_data, aes(PCoA1, PCoA2)) +
    stat_ellipse(
      aes(group = restoration4, fill = restoration4),
      geom = "polygon",
      type = "norm",
      level = 0.68,
      alpha = 0.075,
      color = NA,
      show.legend = FALSE
    ) +
    stat_ellipse(
      aes(group = restoration4, color = restoration4),
      type = "norm",
      level = 0.68,
      linewidth = 0.55,
      alpha = 0.9,
      show.legend = FALSE
    ) +
    geom_hline(yintercept = 0, color = "#D6D6D6", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "#D6D6D6", linewidth = 0.3) +
    geom_point(
      aes(
        fill = restoration4,
        shape = locality,
        size = depth_cm
      ),
      color = "black",
      stroke = 0.38,
      alpha = 0.92
    ) +
    annotate(
      "label",
      x = -Inf,
      y = -Inf,
      label = statistics_label,
      hjust = -0.03,
      vjust = -0.16,
      size = 2.05,
      label.size = 0.18,
      label.padding = grid::unit(1.0, "mm"),
      color = "#222222",
      fill = "white"
    ) +
    scale_fill_manual(
      values = STAGE_COLORS,
      limits = STAGES,
      drop = FALSE
    ) +
    scale_color_manual(
      values = STAGE_COLORS,
      limits = STAGES,
      drop = FALSE
    ) +
    scale_shape_manual(values = c(21, 22, 24)) +
    scale_size_manual(
      values = c("5" = 2.0, "20" = 2.7, "40" = 3.45),
      drop = FALSE
    ) +
    coord_equal(clip = "off") +
    labs(
      tag = panel_tag,
      title = title,
      x = x_label,
      y = y_label,
      fill = "Restoration stage",
      shape = "Locality",
      size = "Depth (cm)",
      caption = caption_label
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(shape = 21, size = 2.7, color = "black")
      ),
      shape = guide_legend(
        order = 2,
        override.aes = list(size = 2.7, fill = "white")
      ),
      size = guide_legend(
        order = 3,
        override.aes = list(shape = 21, fill = "#BDBDBD")
      )
    ) +
    theme_manuscript(base_size = 8.4) +
    theme(
      panel.grid.major.x = element_line(
        color = "#ECECEC",
        linewidth = 0.25
      ),
      panel.grid.major.y = element_line(
        color = "#ECECEC",
        linewidth = 0.25
      ),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.35, "mm"),
      legend.margin = margin(0, 0, 0, 2),
      plot.tag.position = c(0.01, 0.99)
    )
}

save_manuscript_figure <- function(plot, stem, width_mm, height_mm) {
  ggsave(
    file.path(PLOT_DIR, paste0(stem, ".png")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    file.path(PLOT_DIR, paste0(stem, ".pdf")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white"
  )
  ggsave(
    file.path(PLOT_DIR, paste0(stem, ".tiff")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
}

# Fuentes de Figura 1, respetando el orden A-E fijado para el manuscrito.
p_fig1b_hill <- p_alpha
p_fig1c_full_aitchison <- plot_pcoa(
  scores,
  "01_full_species_community_all51",
  "Aitchison_primary",
  panel_tag = "C",
  show_population_caption = FALSE
)

save_manuscript_figure(
  p_fig1b_hill,
  "301_Fig1B_Hill_profile_level",
  width_mm = 100,
  height_mm = 132
)
save_manuscript_figure(
  p_fig1c_full_aitchison,
  "301_Fig1C_PCoA_Aitchison_full_community",
  width_mm = 155,
  height_mm = 110
)

# Figuras suplementarias: sensibilidades de metrica, poblacion y universo.
p_full_bray <- plot_pcoa(
  scores,
  "01_full_species_community_all51",
  "Bray_Hellinger_sensitivity",
  show_population_caption = FALSE
)
p_full_matched <- plot_pcoa(
  scores,
  "02_full_species_community_matched_candidates",
  "Aitchison_primary",
  show_population_caption = FALSE
)
p_raw1031_available <- plot_pcoa(
  scores,
  "03_raw_nonunimodal_candidates_1031_available",
  "Aitchison_primary",
  show_population_caption = FALSE
)
p_raw1031_matched <- plot_pcoa(
  scores,
  "04_raw_nonunimodal_candidates_1031_matched",
  "Aitchison_primary",
  show_population_caption = FALSE
)
p_candidate84_matched <- plot_pcoa(
  scores,
  "05_residual_nonunimodal_candidates_84_matched",
  "Aitchison_primary",
  show_population_caption = FALSE
)
p_raw1031_common_clr <- plot_pcoa(
  scores,
  "06_raw1031_selected_from_full_CLR_all51",
  "Euclidean_selected_full_CLR_sensitivity",
  show_population_caption = FALSE
)
p_candidate84_common_clr <- plot_pcoa(
  scores,
  "07_residual84_selected_from_full_CLR_all51",
  "Euclidean_selected_full_CLR_sensitivity",
  show_population_caption = FALSE
)

save_manuscript_figure(
  p_full_bray,
  "301_FigureS_metric_sensitivity_full_Bray",
  150,
  112
)
save_manuscript_figure(
  p_full_matched,
  "301_FigureS_population_control_full_matched",
  150,
  112
)
save_manuscript_figure(
  p_raw1031_available,
  "301_FigureS_raw1031_available_profiles",
  150,
  112
)
save_manuscript_figure(
  p_raw1031_matched,
  "301_FigureS_raw1031_matched_profiles",
  150,
  112
)
save_manuscript_figure(
  p_candidate84_matched,
  "301_FigureS_residual84_matched_profiles",
  150,
  112
)
save_manuscript_figure(
  p_raw1031_common_clr,
  "301_FigureS_raw1031_selected_full_CLR_all51",
  150,
  112
)
save_manuscript_figure(
  p_candidate84_common_clr,
  "301_FigureS_residual84_selected_full_CLR_all51",
  150,
  112
)

figure_contract <- tibble(
  manuscript_panel = c("Figure 1B", "Figure 1C"),
  content = c(
    "Hill q0/q1/q2 at complete-profile level with adjusted 95% CI",
    "Full-community Aitchison PCoA with stage, locality and depth visible"
  ),
  file_stem = c(
    "301_Fig1B_Hill_profile_level",
    "301_Fig1C_PCoA_Aitchison_full_community"
  ),
  inferential_annotation = c(
    "global stage pseudo-F, permutation P and within-role BH q",
    "stage PERMANOVA R2, pseudo-F, P and stage-dispersion P"
  )
)
write_csv(
  figure_contract,
  file.path(TABLE_DIR, "301_contrato_paneles_Figura1.csv")
)

# -----------------------------------------------------------------------------
# Compuerta final, procedencia y cierre transaccional
# -----------------------------------------------------------------------------

full_primary_row <- stage_results %>%
  filter(
    community == "01_full_species_community_all51",
    metric == "Aitchison_primary",
    model_role == "primary_adjusted_locality_depth"
  )

candidate_common_clr_primary <- stage_results %>%
  filter(
    community %in% c(
      "06_raw1031_selected_from_full_CLR_all51",
      "07_residual84_selected_from_full_CLR_all51"
    ),
    metric == "Euclidean_selected_full_CLR_sensitivity",
    model_role == "primary_adjusted_locality_depth"
  )

expected_pcoa_rows <- sum(vapply(
  community_specs,
  function(spec) 2L * nrow(spec$metadata),
  integer(1)
)) + 2L * nrow(meta)

output_gate <- tibble(
  gate = c(
    "all_distance_objects_finite",
    "expected_stage_test_rows",
    "exactly_one_primary_full_community_test",
    "all_stage_p_values_finite_and_bounded",
    "dispersion_rows_complete",
    "PCoA_rows_complete",
    "alpha_primary_tests_complete",
    "alpha_rarefied_tests_complete",
    "alpha_library_adjusted_tests_complete",
    "no_topN_filter_full_community",
    "candidate_common_CLR_has_two_primary_tests",
    "candidate_common_CLR_uses_identical_all51_population",
    "candidate_common_CLR_uses_identical_permutation_count",
    "LOPO_gate_complete",
    "figure1_source_panels_exported_png_pdf_tiff"
  ),
  pass = c(
    all(vapply(
      community_results,
      function(x) x$distance_finite,
      logical(1)
    )),
    nrow(stage_results) == 24L,
    nrow(full_primary_row) == 1L,
    all(is.finite(stage_results$p_value)) &&
      all(stage_results$p_value >= 0 & stage_results$p_value <= 1),
    nrow(dispersion) == 12L && all(is.finite(dispersion$p_value)),
    nrow(scores) == expected_pcoa_rows &&
      all(is.finite(scores$PCoA1)) &&
      all(is.finite(scores$PCoA2)),
    sum(alpha_tests$analysis_role ==
      "primary_observed_Hill_profile_mean") == 3L,
    sum(alpha_tests$analysis_role ==
      "rarefied_depth_sensitivity") == 3L,
    sum(alpha_tests$analysis_role ==
      "library_size_adjusted_sensitivity") == 3L,
    ncol(counts_full) == length(full_taxa),
    nrow(candidate_common_clr_primary) == 2L,
    n_distinct(candidate_common_clr_primary$analysis_population) == 1L &&
      all(candidate_common_clr_primary$n_samples_analysis == 51L),
    n_distinct(candidate_common_clr_primary$n_unique_permutations) == 1L,
    all(lopo_gate$pass),
    all(file.exists(file.path(
      PLOT_DIR,
      c(
        "301_Fig1B_Hill_profile_level.png",
        "301_Fig1B_Hill_profile_level.pdf",
        "301_Fig1B_Hill_profile_level.tiff",
        "301_Fig1C_PCoA_Aitchison_full_community.png",
        "301_Fig1C_PCoA_Aitchison_full_community.pdf",
        "301_Fig1C_PCoA_Aitchison_full_community.tiff"
      )
    )))
  ),
  detail = c(
    paste0(
      "valid=",
      sum(vapply(
        community_results,
        function(x) x$distance_finite,
        logical(1)
      )),
      "/7"
    ),
    paste0("n=", nrow(stage_results)),
    paste0("n=", nrow(full_primary_row)),
    "stage P values constrained to [0,1]",
    paste0("n=", nrow(dispersion)),
    paste0("n=", nrow(scores), "; expected=", expected_pcoa_rows),
    paste0(
      "n=",
      sum(alpha_tests$analysis_role ==
        "primary_observed_Hill_profile_mean")
    ),
    paste0(
      "n=",
      sum(alpha_tests$analysis_role ==
        "rarefied_depth_sensitivity")
    ),
    paste0(
      "n=",
      sum(alpha_tests$analysis_role ==
        "library_size_adjusted_sensitivity")
    ),
    paste0("all_retained_non_macro_eukaryote_taxa=", ncol(counts_full)),
    paste0("n=", nrow(candidate_common_clr_primary)),
    paste0(
      "populations=",
      paste(
        unique(candidate_common_clr_primary$analysis_population),
        collapse = "|"
      )
    ),
    paste0(
      "unique_permutation_counts=",
      paste(
        unique(candidate_common_clr_primary$n_unique_permutations),
        collapse = "|"
      )
    ),
    paste0("LOPO_gates_passed=", sum(lopo_gate$pass), "/", nrow(lopo_gate)),
    "Figure 1B and 1C source panels present in PNG, PDF and TIFF"
  )
)
write_csv(output_gate, file.path(TABLE_DIR, "301_compuerta_outputs.csv"))
if (!all(output_gate$pass)) {
  stop("Fallo la compuerta final de outputs de 301", call. = FALSE)
}

input_paths <- c(
  PS_RDS,
  META_CSV,
  MEMBERSHIP_CSV,
  METRICS_CSV
)
input_md5 <- tibble(
  input = c(
    "ps_rds",
    "meta_csv",
    "membership_csv",
    "raw_metrics_csv"
  ),
  path = input_paths,
  md5 = unname(as.character(tools::md5sum(input_paths)))
)
write_csv(input_md5, file.path(TABLE_DIR, "301_input_md5.csv"))

run_info_values <- c(
  script = SCRIPT_ID,
  run001 = RUN001,
  run003 = RUN003,
  run201 = RUN201,
  run202 = RUN202,
  ps_rds = PS_RDS,
  meta_csv = META_CSV,
  raw_metrics_csv = METRICS_CSV,
  seed = SEED,
  n_perm_requested = N_PERM,
  n_profile_permutations_unique_all51 = nrow(permutations_rows),
  minimum_attainable_stage_p_all51 = 1 / (nrow(permutations_rows) + 1),
  n_profile_permutations_unique_raw1031 =
    nrow(permutations_rows_raw1031),
  minimum_attainable_stage_p_raw1031 =
    1 / (nrow(permutations_rows_raw1031) + 1),
  n_profile_permutations_unique_candidate_common =
    nrow(permutations_rows_comparison),
  minimum_attainable_stage_p_candidate_common =
    1 / (nrow(permutations_rows_comparison) + 1),
  n_rarefy = N_RAREFY,
  rarefaction_depth = rarefaction_depth,
  pseudocount = PSEUDOCOUNT,
  n_samples_all51 = nrow(meta),
  n_profiles_all51 = n_distinct(meta$profile_id),
  n_samples_raw1031_available = nrow(meta_raw1031),
  n_profiles_raw1031_available = n_distinct(meta_raw1031$profile_id),
  n_samples_candidate_common = nrow(meta_comparison),
  n_profiles_candidate_common = n_distinct(meta_comparison$profile_id),
  n_excluded_profiles_raw1031 = length(excluded_raw1031_profiles),
  excluded_profiles_raw1031 = paste(
    excluded_raw1031_profiles,
    collapse = "|"
  ),
  n_excluded_profiles_candidate_common =
    length(comparison_excluded_profiles),
  excluded_profiles_candidate_common = paste(
    comparison_excluded_profiles,
    collapse = "|"
  ),
  zero_samples_raw1031 = paste(zero_raw1031_samples, collapse = "|"),
  zero_samples_residual84 = paste(
    zero_subset_samples,
    collapse = "|"
  ),
  n_full_taxa_input = nrow(otu),
  n_macro_eukaryote_taxa_excluded = sum(macro_eukaryote),
  n_full_microbial_taxa_analyzed = ncol(counts_full),
  n_raw_candidates1031 = length(raw_candidate_ids),
  n_residual_candidates84 = length(candidate_ids),
  primary_metric = "Aitchison_CLR_pseudocount_1",
  common_CLR_candidate_sensitivity =
    "Euclidean_distance_on_coordinates_selected_after_full_community_CLR",
  metric_sensitivity = "Bray_Curtis_after_Hellinger",
  primary_model = "distance ~ locality + depth_cm + restoration4",
  environmental_sensitivity_model =
    "distance ~ locality + depth_cm + ECI_PC1 + restoration4",
  primary_inferential_unit = "complete_vertical_profile",
  figure1_sources = "B=Hill_profile_level; C=full_Aitchison_PCoA",
  candidate_comparison =
    "raw1031_vs_residual84_on_identical_all51_full_CLR_reference",
  LOPO_sensitivity =
    "leave_one_complete_vertical_profile_out_regenerate_profile_permutations",
  manuscript_order =
    "full community first; candidate common-CLR and LOPO sensitivities second",
  latest_attempt_file = LATEST_ATTEMPT_FILE,
  latest_file = LATEST_FILE
)
run_info <- enframe(
  run_info_values,
  name = "parameter",
  value = "value"
) %>%
  mutate(value = as.character(value))
write_csv(run_info, file.path(TABLE_DIR, "301_parametros_y_resumen.csv"))

report <- c(
  "301 taxonomic community analysis",
  "=================================",
  "",
  paste0("Phyloseq: ", PS_RDS),
  paste0("Metadata: ", META_CSV),
  paste0("RUN201: ", RUN201),
  paste0("RUN202: ", RUN202),
  paste0("Output: ", OUT_DIR),
  "",
  "Manuscript order:",
  paste0(
    "  1. full microbial community (primary): ",
    ncol(counts_full),
    " taxa; ", nrow(meta), " samples / ",
    n_distinct(meta$profile_id), " profiles"
  ),
  paste0(
    "  2. raw non-unimodal candidates: ",
    length(raw_candidate_ids),
    " taxa"
  ),
  paste0(
    "  3. residual non-unimodal candidates: ",
    length(candidate_ids),
    " taxa"
  ),
  "",
  "Candidate-population contract:",
  paste0(
    "  - raw1031 maximum positive population: ",
    nrow(meta_raw1031), " samples / ",
    n_distinct(meta_raw1031$profile_id), " complete profiles"
  ),
  paste0(
    "  - matched 1031-vs-84 population: ",
    nrow(meta_comparison), " samples / ",
    n_distinct(meta_comparison$profile_id), " complete profiles"
  ),
  paste0(
    "  - profiles excluded from matched comparison: ",
    paste(comparison_excluded_profiles, collapse = ", ")
  ),
  "  - the full community is rerun on the matched population as a population control",
  "  - matched 1031 and 84 analyses use identical rows and permutations",
  "  - the common-CLR comparison selects 1031 and 84 coordinates after one full-community CLR",
  "  - the common-CLR comparison retains all 51 samples and all 17 profiles",
  "",
  "Full-community alpha diversity:",
  "  - q0 richness, q1 exponential Shannon and q2 inverse Simpson",
  "  - primary tests use profile means across the three depths",
  "  - Figure 1B shows all 17 profiles plus locality-adjusted means and 95% CI",
  paste0(
    "  - rarefaction sensitivity: ",
    N_RAREFY,
    " iterations at ",
    rarefaction_depth,
    " reads"
  ),
  "  - library-size-adjusted sensitivity uses Freedman-Lane permutations",
  "",
  "Community composition:",
  "  - primary metric: Aitchison distance after CLR with pseudocount 1",
  "  - sensitivity metric: Bray-Curtis after Hellinger transformation",
  "  - primary model adjusts locality and depth",
  "  - environmental sensitivity additionally adjusts ECI_PC1",
  "  - complete profiles are permuted within locality",
  "  - all three depths move together and retain their 5/20/40 order",
  "  - permutations are regenerated for each analysis population",
  "  - LOPO omits all three depths of one profile and regenerates permutations",
  "  - Figure 1C is only the full-community Aitchison PCoA",
  "  - Figure 1C maps stage to color, locality to shape and depth to point size",
  "  - Figure 1C reports R2, pseudo-F, stage P and stage-dispersion P",
  "",
  "Figure-output contract:",
  "  - 301 exports source panels Figure 1B and Figure 1C separately",
  "  - 301 does not assemble Figure 1A-E",
  "  - 1,031- and 84-taxon PCoAs remain supplementary sensitivities",
  "  - no 12-taxon analysis or figure is produced",
  "  - colors follow the original 011/051 manuscript palette",
  "  - outputs are PNG 600 dpi, vector PDF and LZW TIFF 600 dpi",
  "",
  "Interpretation contract:",
  "  - the full-community result is interpreted first",
  "  - raw1031 and residual84 are compared on identical all-51 full-CLR coordinates",
  "  - matched subcomposition analyses remain population/closure sensitivities",
  "  - LOPO is an influence diagnostic, not 17 independent confirmatory tests",
  "  - betadisper is diagnostic, not evidence of centroid separation",
  "  - no result is forced to agree with the preliminary manuscript"
)
writeLines(report, file.path(OUT_DIR, "301_report.txt"))

writeLines(
  capture.output(sessionInfo()),
  file.path(LOG_DIR, "301_sessionInfo.txt")
)

writeLines(OUT_DIR, LATEST_FILE)

cat("===== DONE 301 =====\n")
cat("Samples:", nrow(meta), "\n")
cat("Complete profiles:", n_distinct(meta$profile_id), "\n")
cat("Raw1031 available samples:", nrow(meta_raw1031), "\n")
cat(
  "Raw1031 available complete profiles:",
  n_distinct(meta_raw1031$profile_id),
  "\n"
)
cat("Candidate-common samples:", nrow(meta_comparison), "\n")
cat(
  "Candidate-common complete profiles:",
  n_distinct(meta_comparison$profile_id),
  "\n"
)
cat(
  "Excluded candidate-common profiles:",
  paste(comparison_excluded_profiles, collapse = "|"),
  "\n"
)
cat("Full microbial taxa:", ncol(counts_full), "\n")
cat("Raw candidate1031 taxa:", ncol(counts_raw1031), "\n")
cat("Residual candidate84 taxa:", ncol(counts_candidate84), "\n")
cat("Unique whole-profile permutations:", nrow(permutations_rows), "\n")
cat(
  "Unique raw1031 whole-profile permutations:",
  nrow(permutations_rows_raw1031),
  "\n"
)
cat(
  "Unique candidate-common whole-profile permutations:",
  nrow(permutations_rows_comparison),
  "\n"
)
cat("Rarefaction depth:", rarefaction_depth, "\n")
cat("Output:", OUT_DIR, "\n")
cat("Latest:", LATEST_FILE, "\n")
cat("Finished UTC:", format(Sys.time(), tz = "UTC"), "\n")
