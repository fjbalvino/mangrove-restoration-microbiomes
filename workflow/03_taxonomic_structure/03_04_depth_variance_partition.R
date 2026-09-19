#!/usr/bin/env Rscript
# ============================================================
# 03_04_depth_variance_partition.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Fig. 1d used historical ECI from phyloseq. Its numeric vector differs from 003b; manuscript reconciliation required.
# Inputs (source expressions; complete list in docs/contracts/03_04_depth_variance_partition.json):
#   trimws(readLines(f, n = 1L, warn = FALSE))
#   ps <- readRDS(PS_RDS)
#   master <- read_csv(MASTER_CSV, show_col_types = FALSE)
#   selected <- read_csv(BIMODAL_CSV, show_col_types = FALSE)
# Outputs (source expressions; complete list in contract):
#   writeLines(OUT_DIR, file.path(OUT_ROOT, paste0("LATEST_", SCRIPT_ID, ".txt")))
#   write_csv(
#   write_csv(eci_audit, file.path(TABLE_DIR, "302_auditoria_ECI_PC1.csv"))
#   write_csv(profile_audit, file.path(TABLE_DIR, "302_auditoria_perfiles.csv"))
#   write_csv(results, file.path(TABLE_DIR, "302_resultados_restauracion_por_profundidad.csv"))
#   write_csv(scores_all, file.path(TABLE_DIR, "302_CAP_scores_source_data.csv"))
#   write_csv(dispersion_primary, file.path(TABLE_DIR, "302_dispersion_por_profundidad.csv"))
#   write_csv(full_tests, file.path(TABLE_DIR, "302_auditoria_completa_anova_CAP.csv"))
#   write_csv(varpart_all, file.path(TABLE_DIR, "302_particion_varianza_por_profundidad.csv"))
#   ggsave(
# Algorithmic provenance:
# Depth-specific partial dbRDA and adjusted unique variance fractions; separate test schemes.
#   McArdle & Anderson (2001), doi:10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2; Peres-Neto et al. (2006), doi:10.1890/0012-9658(2006)87[2614:VPOESD]2.0.CO;2.
# Source SHA-256: becd1c666a66930a8e9c68b9deedd67c141fb40971da201635fc07d69a8115c1
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


# 302_CAP_parcial_bimodal_restauracion_profundidad.R
# dbRDA/PERMANOVA parcial por profundidad para dos poblaciones predefinidas:
#   1) comunidad completa: 51 muestras, 17 perfiles, 12,761 especies;
#   2) comunidad bimodal residual: 42 muestras, 14 perfiles, 84 especies.
#
# Especificacion primaria reproducida de
# 53_alt_depth_varpart_horizontal_all_taxa.R:
#   comunidad completa sin cap de taxa y CLR con pseudoconteo 1;
#   dbRDA: restoration4 + ECI_PC1 + Condition(locality);
#   P y pseudo-F secuenciales de anova.cca(..., by = "terms"), con
#   restoration4 antes de ECI_PC1 y permutaciones no restringidas;
#   efecto de restauracion como fraccion pura de R2 ajustado de varpart;
#   valores negativos de varpart conservados en tablas y truncados solo para
#   reproducir la barra normalizada del script 53.
# La PERMANOVA marginal con permutaciones restringidas por localidad se
# conserva como sensibilidad y nunca sustituye al P primario de la replica 53.
# A profundidad fija hay una observacion por perfil. Los tres perfiles excluidos
# del universo bimodal se eliminan completos porque al menos una de sus
# profundidades tiene suma cero en la subcomposicion de 84 especies.

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(permute)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (!length(i)) return(default)
  if (i == length(args)) stop("Argumento sin valor: ", flag, call. = FALSE)
  args[[i + 1L]]
}
as_int <- function(x, default) { z <- suppressWarnings(as.integer(x)); if (is.na(z)) default else z }
as_num <- function(x, default) { z <- suppressWarnings(as.numeric(x)); if (is.na(z)) default else z }
as_flag <- function(x) {
  tolower(str_squish(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}
stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S", tz = "UTC")
read_latest <- function(root, id) {
  f <- file.path(root, paste0("LATEST_", id, ".txt"))
  if (!file.exists(f)) stop("No existe ", f)
  trimws(readLines(f, n = 1L, warn = FALSE))
}
find_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) stop("No se encontro input: ", paste(paths, collapse = "; "))
  hit[[1L]]
}
coalesce_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) return(rep(NA_character_, nrow(df)))
  out <- as.character(df[[hit[[1L]]]])
  if (length(hit) > 1L) for (nm in hit[-1L]) out <- dplyr::coalesce(out, as.character(df[[nm]]))
  out
}
normalize_depth <- function(x) {
  suppressWarnings(as.numeric(str_extract(as.character(x), "[0-9]+(?:\\.[0-9]+)?")))
}
clr_rows <- function(m, pc) {
  z <- m + pc
  log_prop <- log(z / rowSums(z))
  sweep(log_prop, 1, rowMeans(log_prop), "-")
}
find_result_column <- function(df, candidates, context) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) {
    stop(
      "No se encontro una columna compatible para ", context,
      ". Columnas disponibles: ", paste(names(df), collapse = " | "),
      call. = FALSE
    )
  }
  hit[[1L]]
}

OUT_ROOT <- get_arg("--out_root", "/home/fjbalvino/Tipping_points/resultados_finales")
RUN001 <- get_arg("--run001", read_latest(OUT_ROOT, "001_auditar_estructura_metadata"))
PS_RDS <- get_arg("--ps_rds", file.path(RUN001, "rds", "001_phyloseq_canon_51.rds"))
MASTER_CSV <- get_arg("--master_csv", file.path(RUN001, "tables", "001_metadata_canon_51.csv"))
RUN202 <- get_arg("--run202", read_latest(OUT_ROOT, "202_construir_matriz_CLR_y_eje_bimodal"))
BIMODAL_CSV <- get_arg(
  "--bimodal_csv",
  find_file(c(
    file.path(RUN202, "tables", "202_auditoria_mapeo_taxa.csv"),
    file.path(RUN202, "202_auditoria_mapeo_taxa.csv")
  ))
)
N_PERM <- as_int(get_arg("--n_perm", "999"), 999L)
SEED <- as_int(get_arg("--seed", "123"), 123L)
PSEUDOCOUNT <- as_num(get_arg("--pseudocount", "1"), 1)
EXPECTED_FULL_SAMPLES <- as_int(get_arg("--expected_full_samples", "51"), 51L)
EXPECTED_FULL_PROFILES <- as_int(get_arg("--expected_full_profiles", "17"), 17L)
EXPECTED_FULL_TAXA <- as_int(get_arg("--expected_full_taxa", "12761"), 12761L)
EXPECTED_BIMODAL_SAMPLES <- as_int(get_arg("--expected_bimodal_samples", "42"), 42L)
EXPECTED_BIMODAL_PROFILES <- as_int(get_arg("--expected_bimodal_profiles", "14"), 14L)
EXPECTED_BIMODAL_TAXA <- as_int(get_arg("--expected_bimodal_taxa", "84"), 84L)
if (!is.finite(PSEUDOCOUNT) || abs(PSEUDOCOUNT - 1) > .Machine$double.eps^0.5) {
  stop("La replica del script 53 exige --pseudocount 1", call. = FALSE)
}
if (is.na(N_PERM) || N_PERM != 999L) {
  stop("La replica del script 53 exige --n_perm 999", call. = FALSE)
}
if (is.na(SEED) || SEED != 123L) {
  stop("La replica del script 53 exige --seed 123", call. = FALSE)
}

SCRIPT_ID <- "302_CAP_parcial_bimodal_restauracion_profundidad"
OUT_DIR <- file.path(OUT_ROOT, paste0(SCRIPT_ID, "_", stamp()))
TABLE_DIR <- file.path(OUT_DIR, "tables")
PLOT_DIR <- file.path(OUT_DIR, "plots")
LOG_DIR <- file.path(OUT_DIR, "logs")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
writeLines(OUT_DIR, file.path(OUT_ROOT, paste0("LATEST_", SCRIPT_ID, ".txt")))

log_con <- file(file.path(LOG_DIR, "302_log.txt"), "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0L) sink(type = "message")
  while (sink.number() > 0L) sink()
  close(log_con)
}, add = TRUE)
msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"), "|", ..., "\n")
  flush.console()
}

STAGES <- c(
  "Degraded", "Early restoration", "Intermediate restoration",
  "Advanced restoration", "Conserved"
)
DEPTHS <- c(5, 20, 40)
DEPTH_ANALYSIS_ORDER <- c(40, 20, 5)
STAGE_COLORS <- c(
  "Degraded" = "#E41A1C",
  "Early restoration" = "#FDAE61",
  "Intermediate restoration" = "#FFF7A8",
  "Advanced restoration" = "#A6DDA0",
  "Conserved" = "#2B83BA"
)
COMMUNITY_COLORS <- c(
  "full_species_community" = "#4D4D4D",
  "residual_bimodal_84" = "#7B3294"
)
COMMUNITY_SHAPES <- c(
  "full_species_community" = 21,
  "residual_bimodal_84" = 24
)
write_csv(
  tibble(
    restoration4 = STAGES,
    hex_color = unname(STAGE_COLORS[STAGES]),
    source = "same restoration-stage palette used in Figure 1C / script 301"
  ),
  file.path(TABLE_DIR, "302_contrato_colores_etapas.csv")
)
write_csv(
  tibble(
    community = names(COMMUNITY_COLORS),
    hex_color = unname(COMMUNITY_COLORS),
    point_shape = unname(COMMUNITY_SHAPES),
    figure_role = c(
      "complete community: 17 profiles and 12,761 species",
      "candidate bimodal community: 14 profiles and 84 species"
    )
  ),
  file.path(TABLE_DIR, "302_contrato_colores_universos.csv")
)
write_csv(
  tibble(
    component = c(
      "taxon universe", "CLR pseudocount", "dbRDA formula",
      "primary inferential test",
      "permutation restriction", "restoration effect size",
      "multiple testing in original 53", "dispersion test"
    ),
    script_53_alt_specification = c(
      "all taxa; no top-3000 cap",
      "1",
      "restoration + ECI_PC1 + Condition(locality)",
      "anova.cca capscale by=terms; restoration entered first",
      "unrestricted",
      "pure adjusted R2 from varpart [a]",
      "none; raw P < 0.05 used for asterisk",
      "betadisper permutest"
    ),
    corrected_302_specification = c(
      paste0(EXPECTED_FULL_TAXA, " species in complete-community branch"),
      as.character(PSEUDOCOUNT),
      "restoration4 + ECI_PC1 + Condition(locality)",
      "anova.cca capscale by=terms; restoration4 entered first",
      "unrestricted for primary 53 replication",
      "pure adjusted R2 from varpart [a]",
      "raw P retained as primary; BH also exported as sensitivity",
      "betadisper type=median (default in original 53); permutest"
    )
  ),
  file.path(TABLE_DIR, "302_contrato_reproduccion_script_53_all_taxa.csv")
)
set.seed(SEED)

make_permutations <- function(md) {
  control <- permute::how(
    nperm = N_PERM,
    within = permute::Within(type = "free")
  )
  as.matrix(permute::shuffleSet(nrow(md), nset = N_PERM, control = control))
}
unique_perm_n <- function(p) nrow(unique(as.data.frame(p)))

run_depth_model <- function(counts, meta, depth, community) {
  if (!identical(rownames(counts), meta$sample_id)) {
    stop("El orden de counts y metadata no coincide para ", community)
  }
  idx <- which(meta$depth_cm == depth)
  md <- droplevels(meta[idx, , drop = FALSE])
  m <- counts[idx, , drop = FALSE]
  if (nrow(md) != n_distinct(md$profile_id)) {
    stop("Mas de una muestra por perfil a ", depth, " cm")
  }
  X <- model.matrix(~ locality + restoration4 + ECI_PC1, data = md)
  if (qr(X)$rank < ncol(X)) {
    stop(
      "Modelo locality + restoration4 + ECI_PC1 sin rango completo a ",
      depth, " cm"
    )
  }
  if (any(!is.finite(md$ECI_PC1)) || stats::sd(md$ECI_PC1) <= 0) {
    stop("ECI_PC1 ausente o sin variacion a ", depth, " cm")
  }
  if (any(rowSums(m) <= 0)) {
    stop(
      "Hay muestras de suma cero en ", community, " a ", depth,
      " cm; no se permite fabricar una composicion con pseudoconteos"
    )
  }
  clr <- clr_rows(m, PSEUDOCOUNT)
  d <- stats::dist(clr, method = "euclidean")

  cap <- vegan::capscale(
    clr ~ restoration4 + ECI_PC1 + Condition(locality),
    data = md,
    distance = "euclidean",
    add = FALSE
  )
  # Prueba primaria: replica del script 53_alt. restoration4 es el primer
  # termino y su P es secuencial, no marginal respecto de ECI_PC1. No se
  # reinicia aqui la semilla: el bucle externo reproduce el flujo 40 -> 20 -> 5
  # del script original, con set.seed(123) una vez por comunidad.
  p <- make_permutations(md)
  u <- unique_perm_n(p)
  terms <- as.data.frame(anova.cca(cap, by = "terms", permutations = p)) %>%
    rownames_to_column("term") %>% as_tibble()
  names(terms) <- make.names(names(terms))

  X_stage <- model.matrix(~ restoration4, data = md)[, -1, drop = FALSE]
  X_eci <- as.matrix(md$ECI_PC1)
  rownames(X_eci) <- md$sample_id
  colnames(X_eci) <- "ECI_PC1"
  X_loc <- model.matrix(~ locality, data = md)[, -1, drop = FALSE]
  vp <- vegan::varpart(clr, X_stage, X_eci, X_loc)
  ind <- as.data.frame(vp$part$indfract)
  ind$component <- rownames(ind)
  rownames(ind) <- NULL
  adj_col <- find_result_column(
    ind,
    c("Adj.R.square", "Adj.R.squared", "Adjusted.R.square"),
    paste0("R2 ajustado de varpart a ", depth, " cm")
  )
  fraction_value <- function(letter) {
    pattern <- paste0("^\\[", letter, "\\]")
    z <- ind[
      stringr::str_detect(trimws(ind$component), pattern),
      adj_col,
      drop = TRUE
    ]
    if (length(z) != 1L) return(NA_real_)
    as.numeric(z)
  }
  stage_pure_adj_r2 <- fraction_value("a")
  eci_pure_adj_r2 <- fraction_value("b")
  locality_pure_adj_r2 <- fraction_value("c")

  # El 53 llamaba betadisper() sin declarar type: el valor por defecto es
  # median. Este permutest es la segunda y ultima operacion aleatoria del
  # estrato en el flujo original.
  bd <- vegan::betadisper(d, md$restoration4, type = "median")
  bd_permutation <- as.data.frame(
    vegan::permutest(bd, permutations = N_PERM)$tab
  ) %>%
    rownames_to_column("term") %>%
    as_tibble()
  names(bd_permutation) <- make.names(names(bd_permutation))
  rng_after_original53 <- get(".Random.seed", envir = .GlobalEnv)

  # Auditorias y sensibilidades posteriores: se restaurara el estado RNG del
  # flujo original antes de pasar a la siguiente profundidad.
  set.seed(SEED)
  global <- as.data.frame(anova.cca(cap, permutations = N_PERM)) %>%
    rownames_to_column("term") %>% as_tibble()
  set.seed(SEED)
  axes <- as.data.frame(anova.cca(cap, by = "axis", permutations = N_PERM)) %>%
    rownames_to_column("term") %>% as_tibble()
  names(global) <- make.names(names(global))
  names(axes) <- make.names(names(axes))

  set.seed(SEED)
  adonis_terms <- vegan::adonis2(
    d ~ restoration4 + ECI_PC1,
    data = md,
    permutations = N_PERM,
    by = "terms",
    strata = md$locality
  ) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    as_tibble()
  set.seed(SEED)
  adonis_margin <- vegan::adonis2(
    d ~ restoration4 + ECI_PC1,
    data = md,
    permutations = N_PERM,
    by = "margin",
    strata = md$locality
  ) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    as_tibble()
  names(adonis_terms) <- make.names(names(adonis_terms))
  names(adonis_margin) <- make.names(names(adonis_margin))

  cap_stage_row <- terms %>% filter(term == "restoration4")
  cap_eci_row <- terms %>% filter(term == "ECI_PC1")
  if (nrow(cap_stage_row) != 1L || nrow(cap_eci_row) != 1L) {
    stop("No se obtuvieron restoration4 y ECI_PC1 en anova.cca a ", depth, " cm")
  }
  cap_ss_col <- find_result_column(
    cap_stage_row,
    c("Variance", "SumOfSqs", "ChiSquare", "Sum.Sq", "SumOfSquares"),
    paste0("varianza secuencial de restoration4 a ", depth, " cm")
  )
  cap_f_col <- find_result_column(
    cap_stage_row,
    c("F", "pseudo.F", "F.Model"),
    paste0("pseudo-F anova.cca de restoration4 a ", depth, " cm")
  )
  cap_p_col <- find_result_column(
    cap_stage_row,
    c("Pr..F.", "Pr...F..", "P", "p.value"),
    paste0("P anova.cca de restoration4 a ", depth, " cm")
  )
  cap_eci_f_col <- find_result_column(
    cap_eci_row,
    c("F", "pseudo.F", "F.Model"),
    paste0("pseudo-F anova.cca de ECI_PC1 a ", depth, " cm")
  )
  cap_eci_p_col <- find_result_column(
    cap_eci_row,
    c("Pr..F.", "Pr...F..", "P", "p.value"),
    paste0("P anova.cca de ECI_PC1 a ", depth, " cm")
  )

  adonis_stage_row <- adonis_margin %>% filter(term == "restoration4")
  adonis_eci_row <- adonis_margin %>% filter(term == "ECI_PC1")
  if (nrow(adonis_stage_row) != 1L || nrow(adonis_eci_row) != 1L) {
    stop("No se obtuvieron las filas marginales de sensibilidad a ", depth, " cm")
  }
  adonis_r2_col <- find_result_column(
    adonis_stage_row,
    c("R2", "R.squared", "Rsquare"),
    paste0("R2 marginal de sensibilidad a ", depth, " cm")
  )
  adonis_f_col <- find_result_column(
    adonis_stage_row,
    c("F", "pseudo.F", "F.Model"),
    paste0("pseudo-F marginal de sensibilidad a ", depth, " cm")
  )
  adonis_p_col <- find_result_column(
    adonis_stage_row,
    c("Pr..F.", "Pr...F..", "P", "p.value"),
    paste0("P marginal de sensibilidad a ", depth, " cm")
  )

  bd_anova <- as.data.frame(anova(bd)) %>%
    rownames_to_column("term") %>%
    as_tibble()
  names(bd_anova) <- make.names(names(bd_anova))
  assign(".Random.seed", rng_after_original53, envir = .GlobalEnv)

  ss_stage <- as.numeric(cap_stage_row[[cap_ss_col]][[1L]])
  adj <- tryCatch(vegan::RsquareAdj(cap), error = function(e) list(r.squared = NA_real_, adj.r.squared = NA_real_))

  result <- tibble(
    community = community, depth_cm = depth,
    model = "Aitchison ~ restoration4 + ECI_PC1 + Condition(locality)",
    primary_test = paste0(
      "original53 anova.cca by=terms: restoration4 first; ",
      "unrestricted permutations"
    ),
    n_samples_analysis = nrow(md),
    n_profiles_analysis = n_distinct(md$profile_id),
    n_taxa_analysis = ncol(m),
    term = "restoration4", df = cap_stage_row$Df[[1L]],
    sum_squares = ss_stage,
    sum_squares_source_column = cap_ss_col,
    partial_R2 = stage_pure_adj_r2,
    effect_definition = "pure_adjusted_R2_from_varpart_[a]_restoration_given_ECI_and_locality",
    model_R2 = unname(adj$r.squared),
    model_adjusted_R2 = unname(adj$adj.r.squared),
    pseudo_F = as.numeric(cap_stage_row[[cap_f_col]][[1L]]),
    p_value = as.numeric(cap_stage_row[[cap_p_col]][[1L]]),
    eci_pure_adjusted_R2 = eci_pure_adj_r2,
    locality_pure_adjusted_R2 = locality_pure_adj_r2,
    eci_pseudo_F = as.numeric(cap_eci_row[[cap_eci_f_col]][[1L]]),
    eci_p_value = as.numeric(cap_eci_row[[cap_eci_p_col]][[1L]]),
    sensitivity_adonis_margin_R2 = as.numeric(adonis_stage_row[[adonis_r2_col]][[1L]]),
    sensitivity_adonis_margin_F = as.numeric(adonis_stage_row[[adonis_f_col]][[1L]]),
    sensitivity_adonis_margin_p = as.numeric(adonis_stage_row[[adonis_p_col]][[1L]]),
    n_requested_permutations = N_PERM,
    n_permutation_rows_audited = nrow(p),
    n_unique_permutations = u,
    minimum_attainable_p = 1 / (N_PERM + 1), seed = SEED,
    permutation_scheme = "unrestricted_original53_seed_stream_40_20_5"
  )

  site <- as.data.frame(scores(cap, display = "sites", choices = 1:2, scaling = 1)) %>%
    rownames_to_column("sample_id") %>% as_tibble()
  names(site)[2:3] <- c("CAP1", "CAP2")
  site <- site %>% left_join(md, by = "sample_id") %>%
    mutate(community = community, depth_cm_model = depth)

  standardize_dispersion <- function(tab, test_name) {
    f_name <- find_result_column(
      tab,
      c("F", "F.value", "pseudo.F", "F.Model"),
      paste0("F de dispersion (", test_name, ") a ", depth, " cm")
    )
    p_name <- find_result_column(
      tab,
      c("Pr..F.", "Pr...F..", "P", "p.value"),
      paste0("P de dispersion (", test_name, ") a ", depth, " cm")
    )
    tibble(
      community = community,
      depth_cm = depth,
      test = test_name,
      term = tab$term,
      df = tab$Df,
      pseudo_F = as.numeric(tab[[f_name]]),
      p_value = as.numeric(tab[[p_name]]),
      n_requested_permutations = if_else(
        test_name == "permutest_original53", N_PERM, NA_integer_
      ),
      n_unique_permutations = if_else(
        test_name == "permutest_original53", u, NA_integer_
      ),
      minimum_attainable_p = if_else(
        test_name == "permutest_original53", 1 / (N_PERM + 1), NA_real_
      ),
      seed = if_else(test_name == "permutest_original53", SEED, NA_integer_)
    )
  }
  dispersion <- bind_rows(
    standardize_dispersion(bd_permutation, "permutest_original53"),
    standardize_dispersion(bd_anova, "anova_sensitivity")
  )

  full_tables <- bind_rows(
    global %>% mutate(analysis = "capscale", test = "global"),
    terms %>% mutate(analysis = "capscale", test = "terms_sequential"),
    axes %>% mutate(analysis = "capscale", test = "axis"),
    adonis_terms %>% mutate(analysis = "adonis2", test = "terms_sequential"),
    adonis_margin %>% mutate(analysis = "adonis2", test = "margin_blocked_sensitivity")
  ) %>% mutate(community = community, depth_cm = depth)

  varpart_table <- as_tibble(ind) %>%
    mutate(
      community = community,
      depth_cm = depth,
      predictor_set = case_when(
        str_detect(trimws(component), "^\\[a\\]") ~ "restoration_pure",
        str_detect(trimws(component), "^\\[b\\]") ~ "ECI_PC1_pure",
        str_detect(trimws(component), "^\\[c\\]") ~ "locality_pure",
        str_detect(trimws(component), "^\\[h\\]") ~ "residuals",
        TRUE ~ "shared_fraction"
      )
    )

  list(
    result = result,
    scores = site,
    dispersion = dispersion,
    full = full_tables,
    varpart = varpart_table
  )
}

for (f in c(PS_RDS, MASTER_CSV, BIMODAL_CSV)) {
  if (!file.exists(f)) stop("Falta input: ", f)
}
input_paths <- c(PS_RDS, MASTER_CSV, BIMODAL_CSV)
input_info <- file.info(input_paths)
write_csv(
  tibble(
    input_role = c("phyloseq", "canonical_metadata", "bimodal_taxon_map"),
    path = input_paths,
    size_bytes = as.numeric(input_info$size),
    modified_utc = format(input_info$mtime, tz = "UTC", usetz = TRUE),
    md5 = unname(tools::md5sum(input_paths))
  ),
  file.path(TABLE_DIR, "302_input_provenance.csv")
)
ps <- readRDS(PS_RDS)
otu <- as(phyloseq::otu_table(ps), "matrix")
if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)
storage.mode(otu) <- "numeric"
otu[!is.finite(otu)] <- 0

master <- read_csv(MASTER_CSV, show_col_types = FALSE)
ps_meta_raw <- as.data.frame(
  phyloseq::sample_data(ps),
  stringsAsFactors = FALSE
)
# phyloseq::sample_data conserva una subclase que no es compatible con
# vec_slice() de dplyr >= 1.1; aqui se fuerza un data.frame base genuino.
class(ps_meta_raw) <- "data.frame"
ps_sample_ids <- rownames(ps_meta_raw)
if ("sample_id" %in% names(ps_meta_raw)) {
  existing_sample_ids <- as.character(ps_meta_raw$sample_id)
  informative <- !is.na(existing_sample_ids) & nzchar(existing_sample_ids)
  if (any(existing_sample_ids[informative] != ps_sample_ids[informative])) {
    stop(
      "La columna sample_id de sample_data no coincide con sample_names(ps)",
      call. = FALSE
    )
  }
  ps_meta_raw$sample_id <- NULL
}
ps_meta <- ps_meta_raw %>%
  rownames_to_column("sample_id")
msg("Inputs loaded; phyloseq samples:", nrow(ps_meta))
if (!"sample_id" %in% names(master)) stop("master_csv requiere sample_id")
master <- master %>%
  mutate(sample_id = as.character(sample_id))
if (anyDuplicated(master$sample_id)) {
  stop("master_csv contiene sample_id duplicados", call. = FALSE)
}
if (nrow(master) != EXPECTED_FULL_SAMPLES) {
  stop(
    "master_csv debe definir exactamente ", EXPECTED_FULL_SAMPLES,
    " muestras canonicas; encontradas: ", nrow(master),
    call. = FALSE
  )
}

missing_master_in_phyloseq <- setdiff(master$sample_id, ps_meta$sample_id)
extra_phyloseq_samples <- setdiff(ps_meta$sample_id, master$sample_id)
sample_alignment_audit <- tibble(
  sample_id = union(ps_meta$sample_id, master$sample_id)
) %>%
  mutate(
    in_source_phyloseq = sample_id %in% ps_meta$sample_id,
    in_canonical_master = sample_id %in% master$sample_id,
    analysis_action = case_when(
      in_source_phyloseq & in_canonical_master ~ "retain_canonical_sample",
      in_source_phyloseq & !in_canonical_master ~
        "exclude_source_sample_not_in_canonical_master",
      !in_source_phyloseq & in_canonical_master ~
        "ERROR_canonical_sample_missing_from_phyloseq",
      TRUE ~ "ERROR_unclassified"
    )
  ) %>%
  arrange(analysis_action, sample_id)
write_csv(
  sample_alignment_audit,
  file.path(TABLE_DIR, "302_auditoria_alineacion_muestras_fuente.csv")
)
if (length(missing_master_in_phyloseq) > 0L) {
  stop(
    "Faltan en phyloseq muestras canonicas del master: ",
    paste(missing_master_in_phyloseq, collapse = " | "),
    call. = FALSE
  )
}
msg(
  "Canonical alignment:", nrow(master), "samples retained;",
  length(extra_phyloseq_samples), "source-only samples excluded"
)

ps_meta_canonical <- ps_meta %>%
  filter(sample_id %in% master$sample_id)
meta0 <- inner_join(
  ps_meta_canonical,
  master,
  by = "sample_id",
  suffix = c(".ps", ".master")
)
meta0$.source_order_master <- match(meta0$sample_id, master$sample_id)
meta_pre <- meta0 %>%
  transmute(
    sample_id = as.character(sample_id),
    source_order_master = as.integer(.source_order_master),
    locality = str_squish(coalesce_col(meta0, c("locality.master", "locality", "locality.ps"))),
    restoration4 = str_squish(coalesce_col(meta0, c("restoration4.master", "restoration4", "restoration4.ps"))),
    depth_cm = normalize_depth(coalesce_col(meta0, c("depth_cm.master", "depth_cm", "depth_cm.ps"))),
    lat_block = str_squish(coalesce_col(meta0, c("lat_block.master", "lat_block", "lat_block.ps"))),
    ECI_PC1 = suppressWarnings(as.numeric(coalesce_col(
      meta0,
      c("ECI_PC1.master", "ECI_PC1", "ECI_PC1.ps", "eci_pc1.master", "eci_pc1", "eci_pc1.ps")
    )))
  ) %>%
  filter(sample_id %in% colnames(otu)) %>%
  mutate(
    restoration4 = recode(restoration4, Preserved = "Conserved"),
    restoration4 = factor(restoration4, levels = STAGES, ordered = FALSE),
    profile_id = if_else(
      str_detect(lat_block, fixed(locality)), lat_block,
      paste(locality, lat_block, sep = "::")
    )
  )
eci_audit <- meta_pre %>%
  transmute(
    sample_id,
    depth_cm,
    locality,
    ECI_PC1,
    finite_ECI_PC1 = is.finite(ECI_PC1)
  )
write_csv(eci_audit, file.path(TABLE_DIR, "302_auditoria_ECI_PC1.csv"))
if (nrow(meta_pre) != EXPECTED_FULL_SAMPLES || any(!is.finite(meta_pre$ECI_PC1))) {
  stop(
    "ECI_PC1 debe estar disponible y ser finito para las 51 muestras canonicas",
    call. = FALSE
  )
}

meta <- meta_pre %>%
  filter(
    !is.na(restoration4), !is.na(depth_cm), !is.na(locality),
    !is.na(profile_id)
  ) %>%
  arrange(source_order_master) %>%
  mutate(locality = factor(locality), depth_cm = as.numeric(depth_cm))
if (anyDuplicated(meta$sample_id)) stop("sample_id duplicado")

profile_audit <- meta %>%
  group_by(profile_id) %>%
  summarise(
    n = n(), n_locality = n_distinct(locality), n_stage = n_distinct(restoration4),
    depths = paste(sort(depth_cm), collapse = "|"),
    valid = n == 3L && n_locality == 1L && n_stage == 1L && setequal(depth_cm, DEPTHS),
    .groups = "drop"
  )
write_csv(profile_audit, file.path(TABLE_DIR, "302_auditoria_perfiles.csv"))
if (nrow(meta) != EXPECTED_FULL_SAMPLES ||
    n_distinct(meta$sample_id) != EXPECTED_FULL_SAMPLES ||
    nrow(profile_audit) != EXPECTED_FULL_PROFILES ||
    !all(profile_audit$valid)) {
  stop("302 requiere exactamente 51 muestras en 17 perfiles completos 5/20/40", call. = FALSE)
}

otu <- otu[, meta$sample_id, drop = FALSE]
counts_full <- t(otu)
counts_full <- counts_full[, colSums(counts_full) > 0, drop = FALSE]
if (ncol(counts_full) != EXPECTED_FULL_TAXA) {
  stop(
    "La comunidad completa debe contener exactamente ", EXPECTED_FULL_TAXA,
    " taxa; encontrados: ", ncol(counts_full)
  )
}
selected <- read_csv(BIMODAL_CSV, show_col_types = FALSE)
if (!all(c("Original_Taxon_ID", "candidate84", "in_phyloseq") %in% names(selected))) {
  stop(
    "BIMODAL_CSV debe ser 202_auditoria_mapeo_taxa.csv y contener ",
    "Original_Taxon_ID, candidate84 e in_phyloseq"
  )
}
selected <- selected %>%
  mutate(
    candidate84 = as_flag(candidate84),
    in_phyloseq = as_flag(in_phyloseq)
  )
bimodal_taxa <- selected %>%
  filter(candidate84 %in% TRUE, in_phyloseq %in% TRUE) %>%
  pull(Original_Taxon_ID) %>%
  as.character() %>%
  unique() %>%
  intersect(colnames(counts_full))
if (length(bimodal_taxa) != EXPECTED_BIMODAL_TAXA) {
  stop(
    "El universo bimodal residual debe contener exactamente ",
    EXPECTED_BIMODAL_TAXA, " taxa; encontrados: ", length(bimodal_taxa)
  )
}

counts_bimodal_all <- counts_full[, bimodal_taxa, drop = FALSE]
zero_bimodal_samples <- rownames(counts_bimodal_all)[
  rowSums(counts_bimodal_all) <= 0
]
excluded_bimodal_profiles <- meta %>%
  filter(sample_id %in% zero_bimodal_samples) %>%
  distinct(profile_id) %>%
  pull(profile_id) %>%
  as.character()

excluded_bimodal_audit <- meta %>%
  filter(profile_id %in% excluded_bimodal_profiles) %>%
  mutate(
    candidate84_sum = rowSums(counts_bimodal_all[sample_id, , drop = FALSE]),
    zero_sum_trigger = candidate84_sum <= 0,
    exclusion_reason = paste0(
      "complete_profile_removed_because_at_least_one_depth_has_zero_sum_",
      "in_residual_bimodal_84"
    )
  ) %>%
  select(
    profile_id, sample_id, locality, restoration4, depth_cm,
    candidate84_sum, zero_sum_trigger, exclusion_reason
  ) %>%
  arrange(locality, profile_id, depth_cm)
write_csv(
  excluded_bimodal_audit,
  file.path(TABLE_DIR, "302_perfiles_excluidos_universo_bimodal_84.csv")
)

meta_bimodal <- meta %>%
  filter(!profile_id %in% excluded_bimodal_profiles) %>%
  arrange(source_order_master)
counts_bimodal <- counts_bimodal_all[meta_bimodal$sample_id, , drop = FALSE]

bimodal_profile_audit <- meta_bimodal %>%
  group_by(profile_id) %>%
  summarise(
    n = n(),
    n_locality = n_distinct(locality),
    n_stage = n_distinct(restoration4),
    depths = paste(sort(depth_cm), collapse = "|"),
    valid = n == 3L && n_locality == 1L && n_stage == 1L &&
      setequal(depth_cm, DEPTHS),
    .groups = "drop"
  )
write_csv(
  bimodal_profile_audit,
  file.path(TABLE_DIR, "302_auditoria_perfiles_bimodales_14.csv")
)

exclusion_locality_audit <- meta %>%
  distinct(profile_id, locality) %>%
  filter(profile_id %in% excluded_bimodal_profiles) %>%
  count(locality, name = "n_excluded_profiles") %>%
  tidyr::complete(
    locality = levels(meta$locality),
    fill = list(n_excluded_profiles = 0L)
  )
write_csv(
  exclusion_locality_audit,
  file.path(TABLE_DIR, "302_exclusiones_bimodales_por_localidad.csv")
)

if (length(excluded_bimodal_profiles) !=
      EXPECTED_FULL_PROFILES - EXPECTED_BIMODAL_PROFILES ||
    nrow(meta_bimodal) != EXPECTED_BIMODAL_SAMPLES ||
    n_distinct(meta_bimodal$profile_id) != EXPECTED_BIMODAL_PROFILES ||
    !all(bimodal_profile_audit$valid) ||
    any(rowSums(counts_bimodal) <= 0) ||
    nrow(exclusion_locality_audit) != 3L ||
    any(exclusion_locality_audit$n_excluded_profiles != 1L)) {
  stop(
    "La poblacion bimodal debe contener 42 muestras/14 perfiles, ",
    "sin sumas cero y con un perfil excluido por localidad",
    call. = FALSE
  )
}

all_runs <- list()
k <- 1L
analysis_sets <- list(
  full_species_community = list(counts = counts_full, meta = meta),
  residual_bimodal_84 = list(counts = counts_bimodal, meta = meta_bimodal)
)
for (community in names(analysis_sets)) {
  counts <- analysis_sets[[community]]$counts
  meta_analysis <- analysis_sets[[community]]$meta
  set.seed(SEED)
  for (depth in DEPTH_ANALYSIS_ORDER) {
    msg(
      "Starting depth-specific dbRDA:", community,
      "at", paste0(depth, " cm"),
      "with", nrow(meta_analysis), "samples"
    )
    all_runs[[k]] <- run_depth_model(
      counts,
      meta_analysis,
      depth,
      community
    )
    msg(
      "Completed depth-specific dbRDA:", community,
      "at", paste0(depth, " cm")
    )
    k <- k + 1L
  }
}

results <- bind_rows(lapply(all_runs, `[[`, "result")) %>%
  group_by(community) %>%
  mutate(
    supported_original53_raw = p_value < 0.05,
    q_BH_across_depths = p.adjust(p_value, method = "BH"),
    supported_after_BH = q_BH_across_depths <= 0.05
  ) %>%
  ungroup() %>%
  mutate(
    community_label = case_when(
      community == "full_species_community" ~
        paste0(
          "Complete community\n",
          EXPECTED_FULL_PROFILES, " profiles; ", n_taxa_analysis, " species"
        ),
      community == "residual_bimodal_84" ~
        paste0(
          "Candidate bimodal community\n",
          EXPECTED_BIMODAL_PROFILES, " profiles; ", n_taxa_analysis, " species"
        ),
      TRUE ~ community
    ),
    result_label = sprintf(
      "Pure adj. R² = %.3f\nF = %.2f; P53 = %.3f; q = %.3f",
      partial_R2,
      pseudo_F,
      p_value,
      q_BH_across_depths
    )
  )
scores_all <- bind_rows(lapply(all_runs, `[[`, "scores"))
dispersion_all <- bind_rows(lapply(all_runs, `[[`, "dispersion"))
dispersion_primary <- dispersion_all %>% filter(test == "permutest_original53")
dispersion_sensitivity <- dispersion_all %>%
  filter(test == "anova_sensitivity")
full_tests <- bind_rows(lapply(all_runs, `[[`, "full"))
varpart_all <- bind_rows(lapply(all_runs, `[[`, "varpart"))
write_csv(results, file.path(TABLE_DIR, "302_resultados_restauracion_por_profundidad.csv"))
write_csv(scores_all, file.path(TABLE_DIR, "302_CAP_scores_source_data.csv"))
write_csv(dispersion_primary, file.path(TABLE_DIR, "302_dispersion_por_profundidad.csv"))
write_csv(
  dispersion_sensitivity,
  file.path(TABLE_DIR, "302_dispersion_anova_sensibilidad.csv")
)
write_csv(
  dispersion_all,
  file.path(TABLE_DIR, "302_dispersion_auditoria_completa.csv")
)
write_csv(full_tests, file.path(TABLE_DIR, "302_auditoria_completa_anova_CAP.csv"))
write_csv(varpart_all, file.path(TABLE_DIR, "302_particion_varianza_por_profundidad.csv"))

vertical_summary <- results %>%
  arrange(community, depth_cm) %>%
  group_by(community) %>%
  summarise(
    n_profiles = first(n_profiles_analysis),
    n_taxa = first(n_taxa_analysis),
    supported_depths_original53_raw = if_else(
      any(supported_original53_raw),
      paste0(depth_cm[supported_original53_raw], " cm", collapse = " | "),
      "none"
    ),
    supported_depths_BH_sensitivity = if_else(
      any(supported_after_BH),
      paste0(depth_cm[supported_after_BH], " cm", collapse = " | "),
      "none"
    ),
    interpretation_guardrail = case_when(
      any(supported_original53_raw) && !any(supported_after_BH) ~
        paste0(
          "The original-53 raw test detected at least one depth, but none survived ",
          "BH across depths; report both results and do not call the raw result robust."
        ),
      any(supported_original53_raw) && any(supported_after_BH) ~
        paste0(
          "Report the original-53 raw support and the BH-supported depths separately; ",
          "do not infer continuous vertical persistence unless all depths are supported."
        ),
      TRUE ~ paste0(
        "No sampled depth was supported by the original-53 raw test or its BH ",
        "sensitivity; do not claim depth-specific restoration structure."
      )
    ),
    .groups = "drop"
  )
write_csv(
  vertical_summary,
  file.path(TABLE_DIR, "302_resumen_heterogeneidad_vertical.csv")
)

plot_results <- results %>%
  mutate(
    community = factor(
      community,
      levels = c("full_species_community", "residual_bimodal_84")
    ),
    depth_factor = factor(depth_cm, levels = DEPTHS)
  )

pd <- position_dodge(width = 0.28)
p_effect <- ggplot(
  plot_results,
  aes(depth_factor, partial_R2, color = community, shape = community)
) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "#BDBDBD") +
  geom_point(
    data = filter(plot_results, !supported_original53_raw),
    fill = "white",
    size = 3.4,
    stroke = 0.9,
    position = pd
  ) +
  geom_point(
    data = filter(plot_results, supported_original53_raw),
    aes(fill = community),
    size = 3.4,
    stroke = 0.9,
    position = pd
  ) +
  geom_text(
    aes(label = result_label),
    position = pd,
    vjust = -0.75,
    lineheight = 0.95,
    size = 2.75,
    color = "#202020",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = COMMUNITY_COLORS,
    breaks = names(COMMUNITY_COLORS),
    labels = c(
      "Complete community (17 profiles)",
      "Candidate bimodal community (14 profiles)"
    ),
    drop = FALSE
  ) +
  scale_fill_manual(
    values = COMMUNITY_COLORS,
    breaks = names(COMMUNITY_COLORS),
    labels = c(
      "Complete community (17 profiles)",
      "Candidate bimodal community (14 profiles)"
    ),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = COMMUNITY_SHAPES,
    breaks = names(COMMUNITY_SHAPES),
    labels = c(
      "Complete community (17 profiles)",
      "Candidate bimodal community (14 profiles)"
    ),
    drop = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.35))) +
  labs(
    x = "Sediment depth",
    y = expression("Pure adjusted " * R^2 * " for restoration stage"),
    color = NULL,
    fill = NULL,
    shape = NULL,
    caption = paste0(
      "Filled symbols: raw P < 0.05 in the original-53 sequential anova.cca test; ",
      "BH q values are shown as a multiple-testing sensitivity. Permutations are ",
      "unrestricted, restoration4 enters before ECI_PC1, and locality is conditioned. ",
      "Effects are raw pure adjusted fractions from varpart."
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10, base_family = "sans") +
  theme(
    axis.title = element_text(color = "#202020"),
    axis.text = element_text(color = "#202020"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(2, 2, 2, 2),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.key.width = grid::unit(5, "mm"),
    plot.caption = element_text(size = 7.2, color = "#4D4D4D", hjust = 0),
    plot.margin = margin(9, 14, 8, 9)
  ) +
  guides(
    color = guide_legend(override.aes = list(fill = unname(COMMUNITY_COLORS))),
    fill = "none",
    shape = "none"
  )
ggsave(
  file.path(PLOT_DIR, "302_efecto_restauracion_por_profundidad.png"),
  p_effect, width = 7.6, height = 4.7, dpi = 600, bg = "white",
  limitsize = FALSE
)
ggsave(
  file.path(PLOT_DIR, "302_efecto_restauracion_por_profundidad.pdf"),
  p_effect, width = 7.6, height = 4.7, device = grDevices::cairo_pdf,
  bg = "white", limitsize = FALSE
)
ggsave(
  file.path(PLOT_DIR, "302_efecto_restauracion_por_profundidad.tiff"),
  p_effect, width = 7.6, height = 4.7, dpi = 600,
  compression = "lzw", bg = "white", limitsize = FALSE
)

# Replica visual del 53_alt: solo comunidad completa, todos los taxa, fracciones
# puras negativas truncadas exclusivamente para el grafico y normalizadas al
# total de Restoration + Environment + Locality. La tabla conserva ambos valores.
rep53_bd <- dispersion_primary %>%
  filter(community == "full_species_community", term == "Groups") %>%
  transmute(depth_cm, p_betadisper_original53 = p_value)
rep53_summary <- results %>%
  filter(community == "full_species_community") %>%
  transmute(
    depth_cm,
    n = n_samples_analysis,
    p_stage_original53 = p_value,
    q_BH_sensitivity = q_BH_across_depths,
    Restoration = partial_R2,
    Environment = eci_pure_adjusted_R2,
    Locality = locality_pure_adjusted_R2
  ) %>%
  left_join(rep53_bd, by = "depth_cm") %>%
  arrange(match(depth_cm, DEPTH_ANALYSIS_ORDER))
write_csv(
  rep53_summary,
  file.path(TABLE_DIR, "302_replicacion_53_all_taxa_summary.csv")
)
rep53_long <- rep53_summary %>%
  tidyr::pivot_longer(
    cols = c(Restoration, Environment, Locality),
    names_to = "component",
    values_to = "adjusted_R2_raw"
  ) %>%
  group_by(depth_cm) %>%
  mutate(
    adjusted_R2_plot = pmax(0, adjusted_R2_raw),
    denominator_plot = sum(adjusted_R2_plot),
    share_explained_plot = if_else(
      denominator_plot > 0,
      adjusted_R2_plot / denominator_plot,
      NA_real_
    ),
    depth_label = paste0(
      depth_cm, " cm",
      if_else(p_stage_original53 < 0.05, "*", ""),
      if_else(p_betadisper_original53 < 0.05, "\u2020", "")
    )
  ) %>%
  ungroup() %>%
  arrange(depth_cm, component)
if (any(!is.finite(rep53_long$share_explained_plot))) {
  stop("No se pudo normalizar la replica grafica del script 53", call. = FALSE)
}
rep53_depth_labels <- rep53_long %>%
  distinct(depth_cm, depth_label) %>%
  arrange(desc(depth_cm))
rep53_long <- rep53_long %>%
  mutate(
    depth_label = factor(
      depth_label,
      levels = rep53_depth_labels$depth_label
    ),
    component = factor(
      component,
      levels = c("Restoration", "Environment", "Locality")
    )
  )
write_csv(
  rep53_long,
  file.path(TABLE_DIR, "302_replicacion_53_all_taxa_source_data.csv")
)

p_rep53 <- ggplot(
  rep53_long,
  aes(x = share_explained_plot, y = depth_label, fill = component)
) +
  geom_col(
    width = 0.78,
    color = "white",
    linewidth = 0.25,
    position = position_stack(reverse = TRUE)
  ) +
  scale_x_continuous(
    labels = function(z) paste0(round(z * 100), "%"),
    breaks = seq(0, 1, by = 0.25),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    expand = expansion(mult = c(0.06, 0.06))
  ) +
  scale_fill_manual(
    values = c(
      "Environment" = "#E41A1C",
      "Locality" = "#377EB8",
      "Restoration" = "#4DAF4A"
    ),
    breaks = c("Restoration", "Environment", "Locality")
  ) +
  labs(
    x = "Relative share of pure adjusted variance",
    y = NULL,
    fill = NULL
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.key.size = grid::unit(0.42, "cm"),
    legend.spacing.x = grid::unit(0.14, "cm"),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.text.y = element_text(color = "#303030"),
    panel.grid.major.x = element_line(color = "#E6E6E6", linewidth = 0.35),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(4, 5, 2, 4)
  )
ggsave(
  file.path(PLOT_DIR, "302_replicacion_53_varpart_all_taxa.png"),
  p_rep53, width = 7.2, height = 3.0, dpi = 600, bg = "white"
)
ggsave(
  file.path(PLOT_DIR, "302_replicacion_53_varpart_all_taxa.pdf"),
  p_rep53, width = 7.2, height = 3.0, device = grDevices::cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(PLOT_DIR, "302_replicacion_53_varpart_all_taxa.tiff"),
  p_rep53, width = 7.2, height = 3.0, dpi = 600,
  compression = "lzw", bg = "white"
)

p_cap <- ggplot(
  scores_all %>%
    mutate(
      community_label = case_when(
        community == "full_species_community" ~
          "Complete community (17 profiles)",
        community == "residual_bimodal_84" ~
          "Candidate bimodal community (14 profiles)",
        TRUE ~ community
      ),
      community_label = factor(
        community_label,
        levels = c(
          "Complete community (17 profiles)",
          "Candidate bimodal community (14 profiles)"
        )
      )
    ),
  aes(CAP1, CAP2, fill = restoration4, shape = locality)
) +
  geom_point(
    size = 2.6,
    alpha = 0.90,
    color = "#303030",
    stroke = 0.55
  ) +
  facet_grid(community_label ~ depth_cm_model, scales = "free") +
  scale_fill_manual(values = STAGE_COLORS, limits = STAGES, drop = FALSE) +
  scale_shape_manual(values = c(21, 22, 24)) +
  labs(
    x = "CAP1", y = "CAP2",
    fill = "Restoration stage", shape = "Locality"
  ) +
  theme_bw(base_size = 9, base_family = "sans") +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(2, 2, 2, 2),
    strip.text.y = element_text(size = 8),
    plot.margin = margin(8, 12, 8, 8)
  ) +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 1)
  )
ggsave(
  file.path(PLOT_DIR, "302_CAP_por_profundidad_source_plot.png"),
  p_cap, width = 10.5, height = 7.2, dpi = 600, bg = "white",
  limitsize = FALSE
)
ggsave(
  file.path(PLOT_DIR, "302_CAP_por_profundidad_source_plot.pdf"),
  p_cap, width = 10.5, height = 7.2, device = grDevices::cairo_pdf,
  bg = "white", limitsize = FALSE
)
ggsave(
  file.path(PLOT_DIR, "302_CAP_por_profundidad_source_plot.tiff"),
  p_cap, width = 10.5, height = 7.2, dpi = 600,
  compression = "lzw", bg = "white", limitsize = FALSE
)

precision_audit <- results %>%
  select(community, depth_cm, df, sum_squares, partial_R2, pseudo_F, p_value,
         n_requested_permutations, n_unique_permutations,
         minimum_attainable_p, seed, permutation_scheme) %>%
  mutate(
    rounded_F_2digits = round(pseudo_F, 2),
    duplicated_rounded_F = duplicated(paste(community, rounded_F_2digits)) |
      duplicated(paste(community, rounded_F_2digits), fromLast = TRUE),
    note = if_else(
      duplicated_rounded_F,
      "same rounded F can coexist with different full-precision F and permutation p",
      "unique after rounding within community"
    )
  )
write_csv(precision_audit, file.path(TABLE_DIR, "302_auditoria_precision_pseudoF_completa.csv"))

write_csv(
  tibble(
    parameter = c(
      "script", "ps_rds", "master_csv", "run202", "seed_base",
      "n_perm_requested", "n_source_phyloseq_samples",
      "n_source_only_samples_excluded", "n_full_samples", "n_full_profiles",
      "n_full_taxa", "n_bimodal_samples", "n_bimodal_profiles",
      "n_bimodal_taxa", "n_excluded_bimodal_profiles",
      "n_zero_sum_bimodal_samples", "primary_metric", "pseudocount",
      "primary_inferential_test", "primary_effect_size",
      "capscale_add_correction", "inference_unit",
      "multiple_testing_family", "global_repeated_depth_CAP"
    ),
    value = as.character(c(
      SCRIPT_ID, PS_RDS, MASTER_CSV, RUN202, SEED, N_PERM,
      nrow(ps_meta), length(extra_phyloseq_samples),
      nrow(meta), n_distinct(meta$profile_id), ncol(counts_full),
      nrow(meta_bimodal), n_distinct(meta_bimodal$profile_id),
      ncol(counts_bimodal), length(excluded_bimodal_profiles),
      length(zero_bimodal_samples), "Aitchison_CLR", PSEUDOCOUNT,
      paste0(
        "original53_anova.cca_by_terms_restoration4_first_",
        "unrestricted_permutations"
      ),
      "pure_adjusted_R2_from_varpart_[a]",
      "FALSE",
      "profile_at_each_fixed_depth",
      paste0(
        "original53_raw_P_primary;_BH_across_5_20_40_cm_",
        "within_each_community_as_sensitivity"
      ),
      "not_used_for_p_values"
    ))
  ),
  file.path(TABLE_DIR, "302_parametros_y_resumen.csv")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(LOG_DIR, "302_sessionInfo.txt")
)

final_gate <- tibble(
  gate = c(
    "all_canonical_master_samples_present_in_phyloseq",
    "source_only_samples_excluded_before_analysis",
    "full_population_51_samples_17_profiles",
    "full_species_universe_exactly_12761_taxa",
    "bimodal_population_42_samples_14_profiles",
    "residual_bimodal_universe_exactly_84_taxa",
    "three_complete_profiles_excluded",
    "three_zero_sum_samples_trigger_exclusion",
    "one_bimodal_profile_excluded_per_locality",
    "bimodal_population_has_no_zero_sum_samples",
    "ECI_PC1_complete_and_variable",
    "six_depth_specific_models_complete",
    "original53_all_taxa_summary_has_three_depths",
    "all_model_p_values_finite_and_bounded",
    "all_varpart_restoration_fractions_finite",
    "all_depthwise_BH_values_finite_and_bounded",
    "all_original53_permutation_dispersion_tests_complete",
    "all_anova_dispersion_sensitivities_complete",
    "vertical_interpretation_guardrail_written",
    "main_panel_exported_png_pdf_tiff",
    "original53_all_taxa_varpart_exported_png_pdf_tiff",
    "source_CAP_exported_png_pdf_tiff"
  ),
  pass = c(
    length(missing_master_in_phyloseq) == 0L,
    !any(meta$sample_id %in% extra_phyloseq_samples),
    nrow(meta) == EXPECTED_FULL_SAMPLES &&
      n_distinct(meta$profile_id) == EXPECTED_FULL_PROFILES,
    ncol(counts_full) == EXPECTED_FULL_TAXA,
    nrow(meta_bimodal) == EXPECTED_BIMODAL_SAMPLES &&
      n_distinct(meta_bimodal$profile_id) == EXPECTED_BIMODAL_PROFILES,
    ncol(counts_bimodal) == EXPECTED_BIMODAL_TAXA,
    length(excluded_bimodal_profiles) == 3L,
    length(zero_bimodal_samples) == 3L,
    nrow(exclusion_locality_audit) == 3L &&
      all(exclusion_locality_audit$n_excluded_profiles == 1L),
    all(rowSums(counts_bimodal) > 0),
    all(is.finite(meta$ECI_PC1)) &&
      all(vapply(split(meta$ECI_PC1, meta$depth_cm), stats::sd, numeric(1)) > 0),
    nrow(results) == 6L,
    nrow(rep53_summary) == 3L &&
      setequal(rep53_summary$depth_cm, DEPTHS) &&
      all(rep53_summary$n == EXPECTED_FULL_PROFILES),
    all(is.finite(results$p_value) & results$p_value >= 0 & results$p_value <= 1),
    nrow(varpart_all) > 0L && all(is.finite(results$partial_R2)),
    all(
      is.finite(results$q_BH_across_depths) &
        results$q_BH_across_depths >= 0 &
        results$q_BH_across_depths <= 1
    ),
    nrow(filter(
      dispersion_all,
      test == "permutest_original53", term == "Groups"
    )) == 6L,
    nrow(filter(
      dispersion_all,
      test == "anova_sensitivity", term == "Groups"
    )) == 6L,
    nrow(vertical_summary) == 2L &&
      all(nchar(vertical_summary$interpretation_guardrail) > 0L),
    all(file.exists(file.path(
      PLOT_DIR,
      paste0("302_efecto_restauracion_por_profundidad.", c("png", "pdf", "tiff"))
    ))),
    all(file.exists(file.path(
      PLOT_DIR,
      paste0("302_replicacion_53_varpart_all_taxa.", c("png", "pdf", "tiff"))
    ))),
    all(file.exists(file.path(
      PLOT_DIR,
      paste0("302_CAP_por_profundidad_source_plot.", c("png", "pdf", "tiff"))
    )))
  ),
  detail = c(
    paste0("missing=", length(missing_master_in_phyloseq)),
    paste(extra_phyloseq_samples, collapse = "|"),
    paste0("samples=", nrow(meta), "; profiles=", n_distinct(meta$profile_id)),
    paste0("taxa=", ncol(counts_full)),
    paste0(
      "samples=", nrow(meta_bimodal),
      "; profiles=", n_distinct(meta_bimodal$profile_id)
    ),
    paste0("taxa=", ncol(counts_bimodal)),
    paste(excluded_bimodal_profiles, collapse = "|"),
    paste(zero_bimodal_samples, collapse = "|"),
    paste0(
      as.character(exclusion_locality_audit$locality), "=",
      exclusion_locality_audit$n_excluded_profiles,
      collapse = "|"
    ),
    paste0("minimum_row_sum=", min(rowSums(counts_bimodal))),
    paste0(
      "finite=", sum(is.finite(meta$ECI_PC1)),
      "; depth_sd=",
      paste(
        signif(vapply(split(meta$ECI_PC1, meta$depth_cm), stats::sd, numeric(1)), 6),
        collapse = "|"
      )
    ),
    paste0("models=", nrow(results)),
    paste0(
      "depths=", paste(rep53_summary$depth_cm, collapse = "|"),
      "; samples_per_depth=", paste(rep53_summary$n, collapse = "|")
    ),
    "P constrained to [0,1]",
    paste0("finite fractions=", sum(is.finite(results$partial_R2))),
    "BH q constrained to [0,1]",
    paste0(
      "permutest original53 group rows=",
      nrow(filter(
        dispersion_all,
        test == "permutest_original53", term == "Groups"
      ))
    ),
    paste0(
      "anova sensitivity group rows=",
      nrow(filter(
        dispersion_all,
        test == "anova_sensitivity", term == "Groups"
      ))
    ),
    paste(vertical_summary$supported_depths_original53_raw, collapse = " || "),
    "PNG/PDF/TIFF",
    "PNG/PDF/TIFF",
    "PNG/PDF/TIFF"
  )
)
write_csv(final_gate, file.path(TABLE_DIR, "302_compuertas_finales.csv"))
if (!all(final_gate$pass)) {
  stop("Fallo una o mas compuertas finales de 302", call. = FALSE)
}

cat("Depth-specific partial dbRDA completed.\n")
cat("Complete community:", n_distinct(meta$profile_id), "profiles;", ncol(counts_full), "taxa.\n")
cat("Candidate bimodal community:", n_distinct(meta_bimodal$profile_id), "profiles;", ncol(counts_bimodal), "taxa.\n")
cat("Excluded bimodal profiles:", paste(excluded_bimodal_profiles, collapse = " | "), "\n")
cat("Finished UTC:", format(Sys.time(), tz = "UTC"), "\n")
