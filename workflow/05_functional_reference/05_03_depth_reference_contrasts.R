#!/usr/bin/env Rscript
# ============================================================
# 05_03_depth_reference_contrasts.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Locality/depth matched conserved contrasts relative to Degraded; Fig. 3c,d.
# Inputs (source expressions; complete list in docs/contracts/05_03_depth_reference_contrasts.json):
#   z <- read.delim(path, colClasses = "character", check.names = FALSE)
#   r <- read_tsv(files[1], c("status", "pseudocount", "samples", "profiles"))
#   s <- read_tsv(files[2], c(
#   p <- read_tsv(files[3], c(
# Outputs (source expressions; complete list in contract):
#   ggsave(
#   write.table(
#   write_tsv(d, "604_distances_and_local_contrasts.tsv")
#   write_tsv(counts, "604_profile_counts_by_cell.tsv")
#   write_tsv(summary_cells, "604_descriptive_stage_means.tsv")
#   write_tsv(
#   writeLines(OUT, file.path(BASE, paste0("LATEST_", STEM, ".txt")))
# Algorithmic provenance:
# Depth/locality-matched conserved distances and restored-minus-Degraded contrasts.
#   McArdle & Anderson (2001), doi:10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2.
# Source SHA-256: c1f2d5ab09569965079a7ce96392d556afa8f450ed415bdf916e86d93251d84b
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================

suppressPackageStartupMessages(library(ggplot2))

BASE <- "/home/fjbalvino/Tipping_points/resultados_finales"
RUN602 <- file.path(BASE,
  "602_validar_composicion_funcional_splitplot_20260914_230452")
STEM <- "604_distancias_funcionales_por_profundidad"
LOC <- c("Carmen", "Cozumel", "Tuxpan")
ONT <- c("KEGG_ko", "PFAM")
DEPTH <- c(5, 20, 40)
STAGES <- c("Degraded", "Early restoration",
            "Intermediate restoration", "Advanced restoration")
PAL <- setNames(c("#E41A1C", "#FDAE61", "#FFF7A8", "#A6DDA0"), STAGES)

check <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
}

read_tsv <- function(path, cols) {
  check(file.exists(path), paste("Falta:", path))
  z <- read.delim(path, colClasses = "character", check.names = FALSE)
  check(all(cols %in% names(z)), paste("Columnas incompletas:", path))
  z
}

files <- file.path(RUN602, "tables", c(
  "602_run_summary.tsv",
  "602_reference_distance_by_sample.tsv",
  "602_reference_distance_by_profile.tsv"
))

r <- read_tsv(files[1], c("status", "pseudocount", "samples", "profiles"))
check(
  nrow(r) == 1L && r$status == "PASS" &&
  as.numeric(r$pseudocount) == 1 &&
  as.integer(r$samples) == 51L &&
  as.integer(r$profiles) == 17L,
  "Se requiere la corrida PASS con pc=1, 51 muestras y 17 perfiles."
)

meta <- c("sample_id", "profile_id", "locality", "restoration4", "depth_cm")
s <- read_tsv(files[2], c(
  meta, "ontology", "distance_to_local_depth_Conserved", "used_for_inference"
))
p <- read_tsv(files[3], c(
  "ontology", "profile_id", "mean_distance_to_Conserved"
))

check(
  !anyNA(s[c(meta, "ontology")]) &&
  all(vapply(s[c(meta, "ontology")], function(x)
    all(nzchar(trimws(x))), logical(1))),
  "Hay identificadores o metadata vacios."
)

s$depth_cm <- as.numeric(s$depth_cm)
s$distance <- as.numeric(s$distance_to_local_depth_Conserved)
s$used <- toupper(trimws(s$used_for_inference)) == "TRUE"

check(
  all(toupper(trimws(s$used_for_inference)) %in% c("TRUE", "FALSE")) &&
  setequal(s$ontology, ONT) && nrow(s) == 102L,
  "Capas o indicadores invalidos."
)

for (o in ONT) {
  z <- s[s$ontology == o, ]
  design <- unique(z[c("profile_id", "locality", "restoration4")])

  check(
    nrow(z) == 51L && !anyDuplicated(z$sample_id) &&
    nrow(design) == 17L && !anyDuplicated(design$profile_id) &&
    setequal(design$locality, LOC) &&
    all(design$restoration4 %in% c(STAGES, "Conserved")),
    paste("Diseno incorrecto:", o)
  )

  check(
    all(vapply(split(z$depth_cm, z$profile_id), function(d)
      length(d) == 3L && setequal(d, DEPTH), logical(1))),
    "Perfil incompleto."
  )

  for (stage in c("Conserved", "Degraded")) {
    a <- design[design$restoration4 == stage, ]
    check(
      nrow(a) == 3L && setequal(a$locality, LOC),
      paste("Se requiere un perfil por localidad para", stage)
    )
  }

  check(
    all(z$used == (z$restoration4 != "Conserved")) &&
    all(is.finite(z$distance[z$used]) & z$distance[z$used] >= 0) &&
    all(is.na(z$distance[!z$used])),
    "Distancias o autorreferencias invalidas."
  )
}

a <- s[s$ontology == ONT[1], meta]
b <- s[s$ontology == ONT[2], meta]
check(setequal(a$sample_id, b$sample_id), "IDs distintos entre KO y PFAM.")
b <- b[match(a$sample_id, b$sample_id), ]
check(all(a == b), "Metadata distinta entre KO y PFAM.")

d <- s[s$used, c(meta, "ontology", "distance")]

profile_key <- function(x) paste(x$ontology, x$profile_id, sep = "|")
means <- aggregate(distance ~ ontology + profile_id, d, mean)

check(
  nrow(p) == 28L && !anyDuplicated(profile_key(p)) &&
  setequal(profile_key(p), profile_key(means)),
  "Perfiles distintos de 602."
)

check(
  all(abs(
    as.numeric(p$mean_distance_to_Conserved) -
    means$distance[match(profile_key(p), profile_key(means))]
  ) < 1e-7),
  "Las distancias no reproducen las medias originales de 602."
)

# Cada contraste usa Degraded y Conserved de la MISMA localidad y profundidad.
stratum_key <- function(x) {
  paste(x$ontology, x$locality, x$depth_cm, sep = "|")
}

deg <- d[d$restoration4 == "Degraded", ]
ref <- s[s$restoration4 == "Conserved", ]

check(
  nrow(deg) == 18L && !anyDuplicated(stratum_key(deg)) &&
  nrow(ref) == 18L && !anyDuplicated(stratum_key(ref)),
  "Anclas duplicadas."
)

id_deg <- match(stratum_key(d), stratum_key(deg))
id_ref <- match(stratum_key(d), stratum_key(ref))
check(!anyNA(id_deg) && !anyNA(id_ref), "Falta un ancla local por profundidad.")

d$degraded_sample_id <- deg$sample_id[id_deg]
d$degraded_profile_id <- deg$profile_id[id_deg]
d$degraded_distance <- deg$distance[id_deg]
d$conserved_sample_id <- ref$sample_id[id_ref]
d$conserved_profile_id <- ref$profile_id[id_ref]
d$delta_distance <- d$distance - d$degraded_distance
d$is_degraded_baseline <- d$restoration4 == "Degraded"

check(
  all(is.finite(d$delta_distance)) &&
  all(d$delta_distance[d$is_degraded_baseline] == 0),
  "Contraste invalido."
)

d$ontology <- factor(d$ontology, levels = ONT)
d$locality <- factor(d$locality, levels = LOC)
d$restoration4 <- factor(d$restoration4, levels = STAGES)
d$depth <- factor(d$depth_cm, levels = DEPTH, labels = paste(DEPTH, "cm"))

d <- d[order(
  d$ontology, d$locality, d$depth_cm, d$restoration4, d$profile_id
), ]

d$x <- as.integer(d$restoration4)
cells <- split(seq_len(nrow(d)), interaction(
  d$ontology, d$locality, d$depth, d$restoration4, drop = TRUE
))

for (ii in cells) {
  if (length(ii) > 1L) {
    d$x[ii] <- d$x[ii] + seq(-0.18, 0.18, length.out = length(ii))
  }
}

counts <- aggregate(
  profile_id ~ ontology + locality + depth + restoration4, d, length
)
names(counts)[5] <- "n_profiles"
counts$x <- as.integer(counts$restoration4)

check(
  nrow(d) == 84L && nrow(counts) == 72L &&
  all(counts$n_profiles %in% 1:2),
  "Se esperan 42 distancias por capa y 72 celdas entre ambas capas."
)

summary_cells <- aggregate(
  cbind(distance, delta_distance) ~
    ontology + locality + depth_cm + restoration4,
  d, mean
)
names(summary_cells)[5:6] <- c("mean_distance", "mean_delta_distance")

OUT <- file.path(
  BASE, paste0(STEM, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
check(!dir.exists(OUT), "Ya existe la carpeta de salida.")

for (folder in c("figures", "tables", "provenance")) {
  dir.create(file.path(OUT, folder), recursive = TRUE, showWarnings = FALSE)
}

cat("Auditoria OK: pc=1; 51 muestras/17 perfiles; 3 referencias y 14 comparados.\n")
cat("Output:", OUT, "\n")

separate_sites <- function(g) {
  gt <- ggplotGrob(g)
  panels <- gt$layout[grepl("^panel", gt$layout$name), ]
  cols <- sort(unique(panels$l))
  check(length(cols) == 3L, "Se esperaban tres columnas de localidades.")

  for (j in 1:2) {
    left <- max(panels$r[panels$l == cols[j]]) + 1L
    right <- cols[j + 1L] - 1L
    check(left <= right, "Falta espacio entre localidades.")

    gt <- gtable::gtable_add_grob(
      gt,
      grid::segmentsGrob(
        x0 = .5, x1 = .5, y0 = 0, y1 = 1,
        gp = grid::gpar(col = "grey65", lwd = .7)
      ),
      t = min(panels$t), b = max(panels$b),
      l = left, r = right, clip = "off",
      name = paste0("site_separator_", j)
    )
  }
  gt
}

for (o in ONT) {
  for (metric in c("distance", "delta_distance")) {
    z <- d[d$ontology == o, ]
    z$value <- z[[metric]]
    cc <- counts[counts$ontology == o, ]
    delta <- metric == "delta_distance"

    caption <- paste0(
      "17 profiles: 3 local Conserved references and 14 compared; ",
      "n denotes profiles per cell. Depths belong to the same profiles."
    )

    if (delta) {
      caption <- paste0(
        caption,
        "\nNegative values: closer to Conserved than local Degraded at the same depth. ",
        "Degraded zeros define the baseline."
      )
    }

  # FIG3_CD_LEGIBLE_V1
  g <- ggplot(z, aes(x, value))
  if (delta) g <- g + geom_hline(
    yintercept = 0, linetype = "dashed",
    colour = "grey40", linewidth = .55)

  g <- g +
    geom_point(aes(fill = restoration4, shape = locality),
      size = 4.2, stroke = .8, colour = "black") +
    facet_grid(depth ~ locality) +
    scale_x_continuous(
      breaks = 1:4,
      labels = c("Deg.", "Early", "Inter.", "Adv."),
      limits = c(.6, 4.4),
      expand = expansion(mult = 0)) +
    scale_fill_manual(values = PAL, drop = FALSE) +
    scale_shape_manual(
      values = c(Carmen = 21, Cozumel = 22, Tuxpan = 24)) +
    labs(
      x = "Restoration stage",
      caption = NULL,
      y = if (delta) expression(Delta*" distance to Conserved")
          else "Distance to Conserved") +
    theme_classic(base_size = 18, base_family = "sans") +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text.x = element_text(face = "bold", size = 20),
      strip.text.y = element_text(
        face = "bold", size = 18, angle = 0),
      panel.grid.major.y = element_line(
        colour = "grey92", linewidth = .35),
      panel.spacing.x = grid::unit(10, "mm"),
      panel.spacing.y = grid::unit(6, "mm"),
      axis.text = element_text(size = 16, colour = "black"),
      axis.text.x = element_text(
        angle = 35, hjust = 1, vjust = 1),
      axis.title.x = element_text(
        size = 18, margin = margin(t = 8)),
      axis.title.y = element_text(
        size = 18, margin = margin(r = 8)),
      plot.margin = margin(8, 8, 8, 8))

  if (delta) {
    # Rango comun para KO y PFAM.
    lim <- max(10, 10 * ceiling(max(abs(d$delta_distance)) / 10))
    g <- g + scale_y_continuous(
      limits = c(-lim, lim),
      breaks = function(x) pretty(x, n = 3),
      expand = expansion(mult = c(.05, .05)))
  } else {
    g <- g + scale_y_continuous(
      expand = expansion(mult = c(.06, .08)))
  }
  g <- separate_sites(g)
    path <- file.path(
      OUT, "figures", paste0("604_", o, "_", metric, "_by_depth")
    )

    ggsave(
      paste0(path, ".png"), g,
      width = 20, height = 17, units = "cm", dpi = 600, bg = "white"
    )

    if (capabilities("cairo")) {
      ggsave(
        paste0(path, ".pdf"), g,
        width = 20, height = 17, units = "cm", device = cairo_pdf
      )
    } else {
      ggsave(
        paste0(path, ".pdf"), g,
        width = 20, height = 17, units = "cm", useDingbats = FALSE
      )
    }

    cat("Guardado:", basename(path), "[PNG/PDF]\n")
  }
}

write_tsv <- function(z, name, folder = "tables") {
  write.table(
    z, file.path(OUT, folder, name), sep = "\t",
    quote = FALSE, row.names = FALSE, na = "NA"
  )
}

write_tsv(d, "604_distances_and_local_contrasts.tsv")
write_tsv(counts, "604_profile_counts_by_cell.tsv")
write_tsv(summary_cells, "604_descriptive_stage_means.tsv")
write_tsv(
  data.frame(input = files, md5 = unname(tools::md5sum(files))),
  "604_inputs.tsv", "provenance"
)

capture.output(
  sessionInfo(), file = file.path(OUT, "provenance", "sessionInfo.txt")
)

self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (length(self) == 1L && file.exists(self)) {
  check(
    file.copy(self, file.path(OUT, "provenance")),
    "No se pudo guardar el script."
  )
}

write_tsv(
  data.frame(
    status = "PASS", source_run = RUN602, pseudocount = 1,
    samples = 51, profiles = 17, reference_profiles = 3,
    compared_profiles = 14, distances_per_ontology = 42,
    analysis = "descriptive_depth_diagnostic",
    new_hypothesis_tests = 0, HI_HFR_used = FALSE
  ),
  "604_run_summary.tsv"
)

writeLines(OUT, file.path(BASE, paste0("LATEST_", STEM, ".txt")))
cat("PASS\nRevision descriptiva; no se han calculado nuevas pruebas de hipotesis.\n")
