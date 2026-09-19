#!/usr/bin/env Rscript
# ============================================================
# 04_04_plot_global_TITAN.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Renders stored TITAN results; does not repair an obsolete predictor.
# Inputs (source expressions; complete list in docs/contracts/04_04_plot_global_TITAN.json):
#   DEFAULT401 <- if (file.exists(LATEST401)) trimws(readLines(LATEST401, warn = FALSE, n = 1L)) else NA_character_
#   x <- readr::read_csv(ANN_ALL_FILE, show_col_types = FALSE)
#   sumz <- readr::read_csv(SUMZ_FILE, show_col_types = FALSE)
#   support <- readr::read_csv(support_file, show_col_types = FALSE) %>%
# Outputs (source expressions; complete list in contract):
#   ggsave(
#   readr::write_csv(
#   writeLines(
# Algorithmic provenance:
# Plot saved TITAN2 community thresholds and indicators; no model fitting.
#   Baker & King (2010), doi:10.1111/j.2041-210X.2009.00007.x.
# Source SHA-256: 5d0f5acf3cb7ab10c3f5f814c758915dbf58132147efec69b296f976d3a529c4
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


# ==============================================================================
# 404_TITAN2_nature_style_summary_plots.R
#
# Objetivo:
#   Generar figuras listas para manuscrito a partir del análisis TITAN2 global:
#
#   1) Forest plot editorial de umbrales comunitarios filtrados sum(z)
#   2) Ridge/density plot de puntos de cambio taxonómicos por gradiente
#   3) Dotplot compacto de taxa indicadores confiables anotados
#   4) Dotplots suplementarios por gradiente
#
# Ajustes editoriales:
#   - Sin títulos, subtítulos ni captions dentro de la figura.
#   - Tipografía base 9 pt.
#   - Exportación PNG/PDF/SVG; raster a 600 dpi.
#   - Dimensiones en cm, compatibles con A4.
#   - Leyendas compactas, sin duplicación.
#   - Anchos mejorados para evitar recorte de etiquetas.
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(forcats)
  library(grid)
})

has_ggridges <- requireNamespace("ggridges", quietly = TRUE)
has_svglite  <- requireNamespace("svglite", quietly = TRUE)

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

to_num <- function(x) suppressWarnings(as.numeric(x))

nice_gradient_label <- function(x) {
  dplyr::recode(
    x,
    "MHI_local" = "Mangrove Health Index (MHI)",
    "ECI_PC1" = "Environmental condition PC1 (ECI)",
    "water_inundation_PC1" = "Water / inundation",
    "vegetation_landscape_PC1" = "Vegetation / landscape",
    "physicochemical_PC1" = "Physicochemical",
    "nutrients_redox_PC1" = "Nutrients / redox",
    "moisture_stress_PC1" = "Moisture (NDMI)",
    .default = x
  )
}

direction_label <- function(x) {
  dplyr::case_when(
    x %in% c("z-", "1", "-", "negative", "decrease", "decreasing") ~ "z−",
    x %in% c("z+", "2", "+", "positive", "increase", "increasing") ~ "z+",
    TRUE ~ as.character(x)
  )
}

short_tax_label <- function(x, width = 34) {
  x <- as.character(x)
  x <- gsub("\\s+", " ", x)
  stringr::str_wrap(x, width = width)
}

calc_single_tax_plot_height_cm <- function(n_taxa) {
  h <- 4.8 + 0.62 * n_taxa
  h <- max(h, 9.5)
  h <- min(h, 24.5)
  h
}

calc_single_tax_plot_width_cm <- function(max_chars) {
  w <- 12.5 + 0.045 * max_chars
  w <- max(w, 14.0)
  w <- min(w, 18.0)
  w
}

save_plot <- function(p, file_base, width_cm, height_cm) {
  ggsave(
    filename = paste0(file_base, ".png"),
    plot = p,
    width = width_cm,
    height = height_cm,
    units = "cm",
    dpi = FIG_DPI,
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    filename = paste0(file_base, ".pdf"),
    plot = p,
    width = width_cm,
    height = height_cm,
    units = "cm",
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )

  if (has_svglite) {
    ggsave(
      filename = paste0(file_base, ".svg"),
      plot = p,
      width = width_cm,
      height = height_cm,
      units = "cm",
      device = svglite::svglite,
      bg = "white",
      limitsize = FALSE
    )
  }
}

theme_publication <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),

      axis.title = element_text(color = "black", size = base_size + 1),
      axis.text = element_text(color = "black", size = base_size),

      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.30),
      axis.ticks.length = unit(1.6, "mm"),

      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),

      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black", size = base_size),

      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size),
      legend.key = element_blank(),
      legend.key.size = unit(4.0, "mm"),
      legend.spacing.x = unit(2.0, "mm"),
      legend.spacing.y = unit(0.5, "mm"),
      legend.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm"),
      legend.box.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm"),

      plot.margin = margin(t = 3, r = 6, b = 5, l = 6, unit = "mm")
    )
}

# ==============================================================================
# Parameters
# ==============================================================================

RESULTS_DIR <- "/home/fjbalvino/Tipping_points/resultados_finales"
LATEST401 <- file.path(RESULTS_DIR, "LATEST_401_TITAN2_global_taxa_bimodales_gradientes_ambientales.txt")
DEFAULT401 <- if (file.exists(LATEST401)) trimws(readLines(LATEST401, warn = FALSE, n = 1L)) else NA_character_
OUT401 <- get_arg("--out401", Sys.getenv("OUT401", unset = DEFAULT401))
RELIABILITY <- get_arg("--reliability", "95")
TOP_N_TAXA <- get_integer_arg("--top_n_taxa", 25L)
MIN_ABS_Z <- get_numeric_arg("--min_abs_z", 0)
DROP_MOISTURE <- get_integer_arg("--drop_moisture", 0L)

if (is.na(OUT401) || !nzchar(OUT401) || !dir.exists(OUT401)) {
  stop("Usa --out401 /ruta/al/output/401", call. = FALSE)
}

DIR_TABLES <- file.path(OUT401, "tables")
DIR_PLOTS  <- file.path(OUT401, "plots_nature")

safe_mkdir(DIR_PLOTS)

SUMZ_FILE    <- file.path(DIR_TABLES, "402_titan_sumz_thresholds_all_gradients.csv")
ANN_ALL_FILE <- file.path(DIR_TABLES, "403_titan_indval_all_gradients_annotated.csv")
ANN90_FILE   <- file.path(DIR_TABLES, "403_titan_indval_reliable90_annotated.csv")
ANN95_FILE   <- file.path(DIR_TABLES, "403_titan_indval_reliable95_annotated.csv")

stop_if_missing_file(SUMZ_FILE, "402 sumz thresholds")
stop_if_missing_file(ANN_ALL_FILE, "403 annotated all gradients")

FIG_DPI <- 600
FONT_BASE <- 9

FIG_FOREST_W_CM <- 18.0
FIG_FOREST_H_CM <- 9.5

FIG_RIDGE_W_CM  <- 18.0
FIG_RIDGE_H_CM  <- 9.8

FIG_TAXA_W_CM   <- 18.0
FIG_TAXA_H_CM   <- 25.0

choose_indicator_file <- function() {
  if (!RELIABILITY %in% c("90", "95")) {
    stop("--reliability debe ser 90 o 95.")
  }
  cutoff <- as.numeric(RELIABILITY) / 100
  x <- readr::read_csv(ANN_ALL_FILE, show_col_types = FALSE)
  stopifnot(all(c("purity", "reliability") %in% names(x)))
  x <- x %>%
    mutate(purity = to_num(purity), reliability = to_num(reliability)) %>%
    filter(is.finite(purity), is.finite(reliability),
           purity >= cutoff, reliability >= cutoff)
  if (!nrow(x)) stop("No hay indicadores para el filtro solicitado.")
  list(file = ANN_ALL_FILE, mode = RELIABILITY, data = x)
}

chosen <- choose_indicator_file()
ANN_FILE <- chosen$file
RELIABILITY_USED <- chosen$mode

msg("OUT401: ", OUT401)
msg("Using sumz file: ", SUMZ_FILE)
msg("Using indicator file: ", ANN_FILE)
msg("Reliability mode used: ", RELIABILITY_USED)
msg("DROP_MOISTURE: ", DROP_MOISTURE)
msg("TOP_N_TAXA: ", TOP_N_TAXA)
msg("MIN_ABS_Z: ", MIN_ABS_Z)

gradient_order <- c(
  "Mangrove Health Index (MHI)",
  "Nutrients / redox",
  "Physicochemical",
  "Vegetation / landscape",
  "Water / inundation"
)

if (DROP_MOISTURE == 0) {
  gradient_order <- c(
    "Mangrove Health Index (MHI)",
    "Moisture (NDMI)",
    "Nutrients / redox",
    "Physicochemical",
    "Vegetation / landscape",
    "Water / inundation"
  )
}

pal_dir <- c(
  "z−" = "#D55E00",
  "z+" = "#009E73"
)

shape_dir <- c(
  "z−" = 21,
  "z+" = 24
)

# ==============================================================================
# 1. Forest plot
# ==============================================================================

sumz <- readr::read_csv(SUMZ_FILE, show_col_types = FALSE)

# Los umbrales filtrados de esta corrida usan pureza/reliability 0.95.
support_file <- file.path(
  DIR_TABLES, "402_titan_indval_all_gradients.csv"
)
support <- readr::read_csv(support_file, show_col_types = FALSE) %>%
  mutate(
    direction = case_when(
      direction %in% c("z-", "z−") ~ "z−",
      direction == "z+" ~ "z+",
      TRUE ~ NA_character_
    ),
    passes95 = is.finite(purity) & is.finite(reliability) &
      purity >= 0.95 & reliability >= 0.95
  ) %>%
  group_by(gradient, direction) %>%
  summarise(n_indicators = sum(passes95), .groups = "drop")

forest_all <- sumz %>%
  filter(component %in% c("fsumz-", "fsumz+")) %>%
  { if (DROP_MOISTURE == 1) filter(., gradient != "moisture_stress_PC1") else . } %>%
  mutate(
    gradient_label = nice_gradient_label(gradient),
    gradient_label = factor(gradient_label, levels = rev(gradient_order)),
    direction = ifelse(grepl("\\+", component), "z+", "z−"),
    across(c(cp, q05, q10, q50, q90, q95), to_num)
  ) %>%
  left_join(support, by = c("gradient", "direction"))

if (anyNA(forest_all$gradient_label) ||
    anyNA(forest_all$n_indicators)) {
  stop("Falta etiqueta o conteo de indicadores en el forest plot.")
}

forest_all <- forest_all %>%
  mutate(
    y_base = as.numeric(gradient_label),
    y_plot = y_base + ifelse(direction == "z−", 0.14, -0.14),
    direction = factor(direction, levels = c("z−", "z+"))
  )

axis_labels <- forest_all %>%
  group_by(gradient_label, y_base) %>%
  summarise(
    n_minus = sum(n_indicators[direction == "z−"]),
    n_plus = sum(n_indicators[direction == "z+"]),
    .groups = "drop"
  ) %>%
  arrange(y_base) %>%
  mutate(
    label = paste0(
      as.character(gradient_label),
      "  [", n_minus, " / ", n_plus, "]"
    )
  )

sumz_f <- forest_all %>%
  filter(is.finite(cp))

if (any(!is.finite(sumz_f$q05)) ||
    any(!is.finite(sumz_f$q95))) {
  stop("Intervalos incompletos en el forest plot.")
}

readr::write_csv(
  forest_all,
  file.path(DIR_TABLES, "404_threshold_forest_data_used.csv")
)

writeLines(
  c(
    "Labels: [number of z-minus / number of z-plus indicators].",
    "Indicators require purity and reliability >= 0.95.",
    "Points: observed filtered community change points.",
    "Light intervals: bootstrap percentiles 5-95 (central 90%).",
    "Dark intervals: bootstrap percentiles 10-90 (central 80%).",
    "Zero is an axis reference, not a significance threshold.",
    "Gradient axes retain their own scales; widths are not directly comparable."
  ),
  file.path(DIR_TABLES, "404_threshold_forest_legend.txt")
)

p_forest <- ggplot(
  sumz_f,
  aes(
    y = y_plot, x = cp,
    color = direction, fill = direction, shape = direction
  )
) +
  geom_errorbar(
    aes(xmin = q05, xmax = q95),
    orientation = "y",
    width = 0.12, linewidth = 0.65, alpha = 0.30,
    show.legend = FALSE
  ) +
  geom_errorbar(
    aes(xmin = q10, xmax = q90),
    orientation = "y",
    width = 0.07, linewidth = 1.00, alpha = 0.90,
    show.legend = FALSE
  ) +
  geom_point(size = 3.2, stroke = 0.75, color = "black") +
  scale_color_manual(values = pal_dir, drop = FALSE) +
  scale_fill_manual(values = pal_dir, drop = FALSE) +
  scale_shape_manual(values = shape_dir, drop = FALSE) +
  scale_y_continuous(
    breaks = axis_labels$y_base,
    labels = axis_labels$label,
    limits = c(0.5, nrow(axis_labels) + 0.5),
    expand = expansion(mult = 0)
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  labs(
    x = "Environmental gradient value (original axis scale)",
    y = NULL, fill = "Response"
  ) +
  guides(
    color = "none",
    shape = "none",
    fill = guide_legend(
      nrow = 1,
      override.aes = list(
        shape = unname(shape_dir),
        color = "black", size = 3, alpha = 1
      )
    )
  ) +
  theme_publication(base_size = FONT_BASE) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "bottom"
  )

save_plot(
  p_forest,
  file.path(DIR_PLOTS, "404_threshold_forest_nature"),
  width_cm = FIG_FOREST_W_CM,
  height_cm = FIG_FOREST_H_CM
)

# ==============================================================================
# 2. Indicator taxa table
# ==============================================================================

ann <- chosen$data

needed <- c("gradient", "taxon", "tax_label", "direction")
missing_needed <- setdiff(needed, names(ann))

if (length(missing_needed) > 0) {
  stop("Indicator table missing columns: ", paste(missing_needed, collapse = ", "), call. = FALSE)
}

for (cc in c("obs_cp", "cp_50", "cp_10", "cp_90", "cp_05", "cp_95", "z_score", "taxonomy_full")) {
  if (!cc %in% names(ann)) ann[[cc]] <- NA_real_
}

ann <- ann %>%
  { if (DROP_MOISTURE == 1) filter(., gradient != "moisture_stress_PC1") else . } %>%
  mutate(
    gradient_label = nice_gradient_label(gradient),
    gradient_label = factor(gradient_label, levels = rev(gradient_order)),
    direction = direction_label(direction),
    direction = factor(direction, levels = c("z−", "z+")),

    obs_cp = to_num(obs_cp),
    cp_50 = to_num(cp_50),
    cp_10 = to_num(cp_10),
    cp_90 = to_num(cp_90),
    cp_05 = to_num(cp_05),
    cp_95 = to_num(cp_95),
    z_score = to_num(z_score),

    cp_plot = dplyr::coalesce(obs_cp, cp_50),

    tax_label = ifelse(is.na(tax_label) | tax_label == "", taxon, tax_label),
    taxonomy_full = ifelse(is.na(taxonomy_full) | taxonomy_full == "", tax_label, taxonomy_full),

    abs_z = abs(z_score),
    abs_z = ifelse(is.finite(abs_z), abs_z, 1)
  ) %>%
  filter(
    !is.na(gradient_label),
    !is.na(direction),
    is.finite(cp_plot),
    abs_z >= MIN_ABS_Z
  )

readr::write_csv(
  ann,
  file.path(DIR_TABLES, paste0("404_indicator_taxa_used_reliable", RELIABILITY_USED, ".csv"))
)

msg("Indicator taxa available for plotting: ", nrow(ann))

# ==============================================================================
# 3. Ridge / density plot
# ==============================================================================

if (nrow(ann) > 1) {
  if (has_ggridges) {
    p_ridge <- ggplot(
      ann,
      aes(
        x = cp_plot,
        y = gradient_label,
        fill = direction
      )
    ) +
      ggridges::geom_density_ridges(
        alpha = 0.72,
        scale = 1.05,
        rel_min_height = 0.01,
        color = "white",
        linewidth = 0.35,
        show.legend = TRUE
      ) +
      geom_vline(
        xintercept = 0,
        linewidth = 0.35,
        linetype = "dashed",
        color = "grey45"
      ) +
      scale_fill_manual(values = pal_dir, drop = FALSE) +
      scale_x_continuous(
        breaks = pretty_breaks(n = 6),
        expand = expansion(mult = c(0.04, 0.08))
      ) +
      labs(
        x = "Taxon-specific change point",
        y = NULL,
        fill = "Response"
      ) +
      guides(
        fill = guide_legend(
          title = "Response",
          nrow = 1,
          byrow = TRUE,
          override.aes = list(alpha = 0.9)
        )
      ) +
      theme_publication(base_size = FONT_BASE) +
      theme(
        axis.text.y = element_text(size = 9),
        legend.position = "bottom"
      )
  } else {
    p_ridge <- ggplot(
      ann,
      aes(x = cp_plot, fill = direction)
    ) +
      geom_density(alpha = 0.55, linewidth = 0.4) +
      facet_wrap(~ gradient_label, scales = "free_y", ncol = 2) +
      geom_vline(
        xintercept = 0,
        linewidth = 0.35,
        linetype = "dashed",
        color = "grey45"
      ) +
      scale_fill_manual(values = pal_dir, drop = FALSE) +
      labs(
        x = "Taxon-specific change point",
        y = "Density",
        fill = "Response"
      ) +
      guides(
        fill = guide_legend(
          title = "Response",
          nrow = 1,
          byrow = TRUE
        )
      ) +
      theme_publication(base_size = FONT_BASE)
  }

  save_plot(
    p_ridge,
    file.path(DIR_PLOTS, paste0("404_indicator_cp_ridge_nature_reliable", RELIABILITY_USED)),
    width_cm = FIG_RIDGE_W_CM,
    height_cm = FIG_RIDGE_H_CM
  )
} else {
  msg("Skipping ridge: not enough indicator taxa.")
}

# ==============================================================================
# 4. Top annotated indicator taxa: combined plot
# ==============================================================================

# PANEL_B_STRATIFIED_INTERVALS_V1
# Seleccion reproducible por gradiente y direccion; desempate por taxon ID.

stopifnot(
  !anyDuplicated(ann[, c("gradient", "taxon")]),
  all(is.finite(ann$z_score)),
  all(is.finite(ann$obs_cp)),
  all(is.finite(ann$cp_05)),
  all(is.finite(ann$cp_10)),
  all(is.finite(ann$cp_90)),
  all(is.finite(ann$cp_95)),
  all(ann$cp_05 <= ann$cp_10),
  all(ann$cp_10 <= ann$cp_90),
  all(ann$cp_90 <= ann$cp_95)
)

ann_top <- ann %>%
  mutate(
    abs_z = abs(z_score),
    cp_plot = obs_cp,
    gradient_label = factor(
      nice_gradient_label(gradient), levels = gradient_order
    )
  ) %>%
  group_by(gradient, direction) %>%
  arrange(desc(abs_z), taxon, .by_group = TRUE) %>%
  slice_head(n = 3L) %>%
  ungroup() %>%
  arrange(gradient_label, direction, cp_plot, taxon) %>%
  mutate(row_id = paste(gradient, taxon, sep = "::"))

stopifnot(nrow(ann_top) > 0, !anyDuplicated(ann_top$row_id))

# IDs internos distintos evitan fusionar taxa con etiquetas coincidentes.
ann_top <- ann_top %>%
  group_by(gradient, tax_label) %>%
  mutate(
    display_label = ifelse(
      n() > 1L, paste0(tax_label, " [", taxon, "]"), tax_label
    )
  ) %>%
  ungroup() %>%
  mutate(
    row_id = factor(row_id, levels = rev(row_id)),
    display_label = stringr::str_wrap(display_label, width = 43)
  )

label_map <- setNames(
  ann_top$display_label, as.character(ann_top$row_id)
)

readr::write_csv(
  ann_top,
  file.path(DIR_TABLES, paste0(
    "404_top_indicator_taxa_used_reliable", RELIABILITY_USED, ".csv"
  ))
)

selection_audit <- ann %>%
  count(gradient, direction, name = "n_available") %>%
  left_join(
    ann_top %>% count(gradient, direction, name = "n_displayed"),
    by = c("gradient", "direction")
  )

readr::write_csv(
  selection_audit,
  file.path(DIR_TABLES, "404_panelB_selection_audit.csv")
)

make_indicator_plot <- function(dat, faceted = TRUE) {
  p <- ggplot(
    dat,
    aes(x = cp_plot, y = row_id,
        colour = direction, fill = direction, shape = direction)
  ) +
    geom_segment(
      aes(x = cp_05, xend = cp_95, yend = row_id),
      linewidth = 0.65, alpha = 0.30, show.legend = FALSE
    ) +
    geom_segment(
      aes(x = cp_10, xend = cp_90, yend = row_id),
      linewidth = 1, alpha = 0.90, show.legend = FALSE
    ) +
    geom_point(size = 2.5, stroke = 0.5, colour = "black") +
    scale_y_discrete(
      labels = function(z) unname(label_map[z]),
      expand = expansion(add = 0.65)
    ) +
    scale_x_continuous(
      breaks = pretty_breaks(n = 5),
      expand = expansion(mult = c(0.04, 0.04))
    ) +
    scale_colour_manual(values = pal_dir, drop = FALSE) +
    scale_fill_manual(values = pal_dir, drop = FALSE) +
    scale_shape_manual(values = shape_dir, drop = FALSE) +
    labs(
      x = "Taxon-specific change point (environmental axis units)",
      y = NULL, fill = "Response"
    ) +
    guides(
      colour = "none", shape = "none",
      fill = guide_legend(
        nrow = 1,
        override.aes = list(
          shape = unname(shape_dir), colour = "black", size = 2.5
        )
      )
    ) +
    theme_publication(base_size = 9) +
    theme(
      axis.text.y = element_text(size = 8, lineheight = 1),
      axis.text.x = element_text(size = 8),
      strip.text.y = element_text(angle = 0, size = 8),
      strip.text = element_text(face = "bold"),
      panel.spacing.y = unit(3, "mm"),
      legend.position = "bottom",
      plot.margin = margin(5, 5, 5, 5, unit = "mm")
    )

  if (faceted) {
    p <- p + facet_grid(
      gradient_label ~ .,
      scales = "free_y", space = "free_y", switch = "y"
    ) +
      theme(
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, size = 8)
      )
  }
  p
}

# Una columna para evitar repetir todos los taxa en todas las facetas.
# El alto se adapta al numero de filas y a las etiquetas multilínea.
label_lines <- sum(stringr::str_count(ann_top$display_label, "\n") + 1L)
height_cm <- max(
  24,
  5 + 0.45 * label_lines + 0.8 * n_distinct(ann_top$gradient)
)

p_taxa <- make_indicator_plot(ann_top)

save_plot(
  p_taxa,
  file.path(DIR_PLOTS, paste0(
    "404_top_indicator_taxa_nature_reliable", RELIABILITY_USED
  )),
  width_cm = 23,
  height_cm = height_cm
)

for (g in unique(ann_top$gradient)) {
  dat <- ann_top %>% filter(gradient == g)
  p_g <- make_indicator_plot(dat, faceted = FALSE)
  n_lines <- sum(stringr::str_count(dat$display_label, "\n") + 1L)

  save_plot(
    p_g,
    file.path(DIR_PLOTS, paste0(
      "404_panelB_", g, "_reliable", RELIABILITY_USED
    )),
    width_cm = 18,
    height_cm = max(8, 4 + 0.55 * n_lines)
  )
}

writeLines(
  c(
    "Taxon-specific environmental change points.",
    paste0(
      "Up to three taxa per environmental gradient and response direction ",
      "were selected by descending absolute TITAN2 z-score among indicators ",
      "with purity and reliability >= ", as.numeric(RELIABILITY_USED) / 100,
      ". Ties were resolved by taxon identifier."
    ),
    "Points show observed change points; point size is constant.",
    paste(
      "Light and dark segments represent bootstrap percentiles 5-95 and",
      "10-90, respectively (central 90% and 80% intervals)."
    ),
    paste(
      "The displayed set is an illustrative selection.",
      "All qualifying associations are retained in the indicator source table."
    ),
    paste(
      "Environmental axes retain their fitted units.",
      "Equal numerical values across different axes do not imply",
      "equivalent environmental conditions."
    )
  ),
  file.path(DIR_TABLES, "404_panelB_legend.txt")
)

msg("Panel B: ", nrow(ann_top), " associations displayed.")

# ==============================================================================
# Index
# ==============================================================================

fig_tbl <- tibble(
  figure = c(
    "threshold_forest",
    "indicator_cp_ridge",
    "top_indicator_taxa"
  ),
  file_png = c(
    file.path(DIR_PLOTS, "404_threshold_forest_nature.png"),
    file.path(DIR_PLOTS, paste0("404_indicator_cp_ridge_nature_reliable", RELIABILITY_USED, ".png")),
    file.path(DIR_PLOTS, paste0("404_top_indicator_taxa_nature_reliable", RELIABILITY_USED, ".png"))
  ),
  exists = file.exists(file_png),
  interpretation = c(
    "Main summary of filtered TITAN2 community thresholds by environmental gradient and response direction.",
    "Density/synchrony of taxon-specific change points across gradients.",
    "Taxonomic identity of strongest annotated indicator taxa contributing to thresholds."
  )
)

readr::write_csv(
  fig_tbl,
  file.path(DIR_TABLES, "404_nature_style_figures_index.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(DIR_TABLES, "404_sessionInfo.txt")
)

msg("Done.")
msg("Figures written to: ", DIR_PLOTS)
msg("Index: ", file.path(DIR_TABLES, "404_nature_style_figures_index.csv"))
