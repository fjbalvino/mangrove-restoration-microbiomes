#!/usr/bin/env Rscript
# ============================================================
# 05_02_partial_dbRDA.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Constrained ordination of the frozen 602 functional matrices; Fig. 3a,b.
# Inputs (source expressions; complete list in docs/contracts/05_02_partial_dbRDA.json):
#   check <- read.delim(inputs[["summary"]], check.names = FALSE)
#   prep <- readRDS(inputs[["preprocessed"]])
#   previous <- readRDS(inputs[["analysis"]])
# Outputs (source expressions; complete list in contract):
#   writeLines(OUT, file.path(BASE, paste0("LATEST_ATTEMPT_", SCRIPT_ID, ".txt")))
#   write.table(
#   write_tsv(md[, required], "603_metadata_used.tsv")
#   write_tsv(
#   saveRDS(perms, file.path(OUT, "rds", "603_permutations.rds"))
#   ggsave(
#   write_tsv(result_tbl, "603_partial_dbRDA_results.tsv")
#   write_tsv(do.call(rbind, scores_out), "603_sample_scores.tsv")
#   write_tsv(do.call(rbind, axes_out), "603_axis_variance.tsv")
#   saveRDS(models, file.path(OUT, "rds", "603_partial_dbRDA_models.rds"))
# Algorithmic provenance:
# Partial dbRDA of the frozen 602 functional matrices, conditioning on locality and depth.
#   McArdle & Anderson (2001), doi:10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2.
# Source SHA-256: 3995341c2c2a70d0c3683a4c842dd710a0ed669dfcf3b6b6fb2a23ab57c3ed2a
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

library(vegan)
library(ggplot2)

BASE <- "/home/fjbalvino/Tipping_points/resultados_finales"
RUN602 <- file.path(
  BASE,
  "602_validar_composicion_funcional_splitplot_20260914_230452"
)
SCRIPT_ID <- "603_dbRDA_funcional_parcial_restauracion"
SEED <- 123L

inputs <- c(
  preprocessed = file.path(
    RUN602, "rds", "602_primary_preprocessed_data.rds"
  ),
  analysis = file.path(
    RUN602, "rds", "602_analysis_results.rds"
  ),
  summary = file.path(RUN602, "tables", "602_run_summary.tsv")
)
stopifnot(all(file.exists(inputs)))

check <- read.delim(inputs[["summary"]], check.names = FALSE)
stopifnot(
  nrow(check) == 1L, check$status == "PASS",
  check$pseudocount == 1, check$samples == 51L,
  check$profiles == 17L
)

prep <- readRDS(inputs[["preprocessed"]])
previous <- readRDS(inputs[["analysis"]])
md <- prep$metadata

stages <- c(
  "Degraded", "Early restoration", "Intermediate restoration",
  "Advanced restoration", "Conserved"
)
localities <- c("Carmen", "Cozumel", "Tuxpan")

required <- c(
  "sample_id", "profile_id", "locality", "depth_cm", "restoration4"
)
stopifnot(
  all(required %in% names(md)),
  !anyNA(md[, required]),
  nrow(md) == 51L,
  !anyDuplicated(md$sample_id),
  identical(
    as.character(md$sample_id),
    as.character(previous$metadata$sample_id)
  ),
  all(md$restoration4 %in% stages),
  all(md$locality %in% localities)
)

md$locality <- factor(md$locality, levels = localities)
md$restoration4 <- factor(md$restoration4, levels = stages)
md$depth_cm <- as.numeric(as.character(md$depth_cm))
md$depth_f <- factor(md$depth_cm, levels = c(5, 20, 40))
rownames(md) <- md$sample_id
stopifnot(!anyNA(md$depth_f))

profiles <- split(seq_len(nrow(md)), md$profile_id, drop = TRUE)
stopifnot(
  length(profiles) == 17L,
  all(vapply(profiles, function(i) {
    length(i) == 3L &&
      setequal(md$depth_cm[i], c(5, 20, 40)) &&
      length(unique(md$locality[i])) == 1L &&
      length(unique(md$restoration4[i])) == 1L
  }, logical(1)))
)

design <- model.matrix(
  ~ locality + depth_f + restoration4, data = md
)
stopifnot(qr(design)$rank == ncol(design))

perms <- as.matrix(previous$permutations)
stopifnot(
  ncol(perms) == nrow(md),
  nrow(perms) > 0L,
  nrow(unique(as.data.frame(perms))) == nrow(perms)
)
valid_perm <- apply(perms, 1L, function(p) {
  identical(sort(as.integer(p)), seq_len(nrow(md))) &&
    all(md$locality[p] == md$locality) &&
    all(md$depth_cm[p] == md$depth_cm) &&
    all(vapply(profiles, function(i) {
      length(unique(md$profile_id[p[i]])) == 1L
    }, logical(1)))
})
stopifnot(all(valid_perm))

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S", tz = "UTC")
OUT <- file.path(BASE, paste0(SCRIPT_ID, "_", stamp))
if (dir.exists(OUT)) stop("Output directory already exists")
for (sub in c("tables", "figures", "rds", "provenance")) {
  dir.create(file.path(OUT, sub), recursive = TRUE)
}
writeLines(OUT, file.path(BASE, paste0("LATEST_ATTEMPT_", SCRIPT_ID, ".txt")))

write_tsv <- function(x, name) {
  write.table(
    x, file.path(OUT, "tables", name),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}
write_tsv(md[, required], "603_metadata_used.tsv")
write_tsv(
  data.frame(
    permutation = seq_len(nrow(perms)), valid = valid_perm
  ),
  "603_permutation_audit.tsv"
)
saveRDS(perms, file.path(OUT, "rds", "603_permutations.rds"))

cat("Output:", OUT, "\n")
cat("Samples:", nrow(md), "| Profiles:", length(profiles), "\n")
cat("Unique permutations:", nrow(perms), "\n")

pal <- setNames(
  c("#E41A1C", "#FDAE61", "#FFF7A8", "#A6DDA0", "#2B83BA"),
  stages
)
shapes <- c(Carmen = 21, Cozumel = 22, Tuxpan = 24)
sizes <- c("5" = 2, "20" = 2.7, "40" = 3.45)

results <- list()
scores_out <- list()
axes_out <- list()
models <- list()

for (ontology in c("KEGG_ko", "PFAM")) {
  item <- prep[[ontology]]
  m <- item$abundance
  stopifnot(
    is.matrix(m),
    identical(rownames(m), md$sample_id),
    all(is.finite(m)), all(m >= 0),
    all(rowSums(m) > 0)
  )

  for (pc in c(1, 0.5)) {
    role <- if (pc == 1) "primary" else "pseudocount_sensitivity"
    key <- paste0(ontology, "_pc", pc)
    cat("\nRunning:", key, "\n")

    z <- log(m + pc)
    z <- sweep(z, 1L, rowMeans(z), "-")
    d <- dist(z)

    if (pc == 1) {
      stopifnot(
        max(abs(as.numeric(d) - as.numeric(item$distance))) < 1e-8
      )
    }

    fit <- dbrda(
      d ~ restoration4 + Condition(locality + depth_f),
      data = md
    )
    set.seed(SEED)
    test <- anova(fit, permutations = perms, model = "reduced")
    adj <- RsquareAdj(fit)

    total <- fit$tot.chi
    constrained <- fit$CCA$tot.chi
    residual <- fit$CA$tot.chi
    r2_unique <- constrained / total

    expected_id <- if (pc == 1) {
      "primary_common_CLR_pc1"
    } else {
      "sensitivity_common_CLR_pc0.5"
    }

    old <- previous$permanova[
      previous$permanova$ontology == ontology &
        previous$permanova$analysis_id == expected_id &
        previous$permanova$stage_model == "factor",
      , drop = FALSE
    ]
    old_disp <- previous$betadisper[
      previous$betadisper$ontology == ontology &
        previous$betadisper$analysis_id == expected_id,
      , drop = FALSE
    ]
    stopifnot(nrow(old) == 1L, nrow(old_disp) == 1L)

    F_obs <- unname(test$F[1L])
    P_obs <- unname(test$`Pr(>F)`[1L])
    stopifnot(all(is.finite(c(F_obs, P_obs, r2_unique, adj$adj.r.squared))))

    same_F <- abs(F_obs - old$pseudo_f) < 1e-7
    same_R2 <- abs(r2_unique - old$r2) < 1e-7
    same_P <- abs(P_obs - old$p_value) < 1e-10

    if (!same_F || !same_R2) {
      stop("Unexpected F/R2 mismatch with 602: ", key)
    }
    if (!same_P) {
      warning("Permutation P differs from 602; inspect comparison: ", key)
    }

    results[[key]] <- data.frame(
      ontology = ontology,
      role = role,
      pseudocount = pc,
      n_features = ncol(m),
      n_samples = nrow(md),
      n_profiles = length(profiles),
      df_stage = test$Df[1L],
      df_residual = test$Df[2L],
      pseudo_f = F_obs,
      p_value = P_obs,
      r2_unique_fraction_of_total = r2_unique,
      adjusted_r2_unique_fraction_of_total = adj$adj.r.squared,
      r2_fraction_after_conditioning = constrained / (constrained + residual),
      n_permutations = nrow(perms),
      minimum_attainable_p = 1 / (nrow(perms) + 1),
      prior_permanova_p = old$p_value,
      prior_dispersion_p = old_disp$p_value,
      concordant_F = same_F,
      concordant_R2 = same_R2,
      concordant_P = same_P
    )
    print(results[[key]])

    eig <- eigenvals(fit, model = "constrained")
    stopifnot(length(eig) >= 2L, all(eig[1:2] > 0))
    axes_out[[key]] <- data.frame(
      ontology = ontology,
      pseudocount = pc,
      axis = names(eig),
      eigenvalue = as.numeric(eig),
      percent_constrained = 100 * as.numeric(eig) / constrained,
      percent_total = 100 * as.numeric(eig) / total
    )

    # Response-based sample scores; not fitted LC scores.
    xy <- scores(fit, display = "sites", choices = 1:2, scaling = 1)
    xy <- xy[md$sample_id, , drop = FALSE]
    stopifnot(nrow(xy) == 51L, all(is.finite(xy)))
    plotdata <- data.frame(
      md[, required],
      axis1 = xy[, 1L],
      axis2 = xy[, 2L],
      ontology = ontology,
      pseudocount = pc
    )
    plotdata$depth_cm <- factor(plotdata$depth_cm, levels = c(5, 20, 40))
    scores_out[[key]] <- plotdata
    models[[key]] <- fit

    if (pc == 1) {
      percent <- 100 * eig[1:2] / constrained

      p <- ggplot(plotdata, aes(axis1, axis2)) +
        geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
        geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
        stat_ellipse(
          aes(group = restoration4, fill = restoration4),
          geom = "polygon", type = "norm", level = 0.68,
          alpha = 0.075, color = NA, show.legend = FALSE
        ) +
        stat_ellipse(
          aes(group = restoration4, color = restoration4),
          type = "norm", level = 0.68,
          linewidth = 0.55, show.legend = FALSE
        ) +
        geom_point(
          aes(fill = restoration4, shape = locality, size = depth_cm),
          color = "black", stroke = 0.55, alpha = 0.9
        ) +
        scale_fill_manual(values = pal, limits = stages, drop = FALSE) +
        scale_color_manual(values = pal, limits = stages, drop = FALSE) +
        scale_shape_manual(values = shapes, drop = FALSE) +
        scale_size_manual(values = sizes, drop = FALSE) +
        coord_equal() +
        labs(
          x = sprintf("dbRDA1 (%.1f%% of constrained variation)", percent[1]),
          y = sprintf("dbRDA2 (%.1f%% of constrained variation)", percent[2]),
          fill = "Restoration stage",
          shape = "Locality",
          size = "Depth (cm)"
        ) +
        guides(
          color = "none",
          fill = guide_legend(
            order = 1, override.aes = list(shape = 21, size = 3)
          ),
          shape = guide_legend(
            order = 2, override.aes = list(fill = "white", size = 3)
          ),
          size = guide_legend(
            order = 3, override.aes = list(shape = 21, fill = "grey70")
          )
        ) +
        theme_classic(base_size = 11) +
        theme(
          legend.position = "right",
          panel.grid.major = element_line(color = "grey93", linewidth = 0.3)
        )

      for (ext in c("png", "pdf")) {
        ggsave(
          file.path(
            OUT, "figures",
            paste0("603_", ontology, "_partial_dbRDA_pc1.", ext)
          ),
          p, width = 24, height = 17, units = "cm",
          dpi = 400, bg = "white"
        )
      }
    }
  }
}

result_tbl <- do.call(rbind, results)
write_tsv(result_tbl, "603_partial_dbRDA_results.tsv")
write_tsv(do.call(rbind, scores_out), "603_sample_scores.tsv")
write_tsv(do.call(rbind, axes_out), "603_axis_variance.tsv")
saveRDS(models, file.path(OUT, "rds", "603_partial_dbRDA_models.rds"))

write.table(
  data.frame(
    input = names(inputs),
    path = unname(normalizePath(inputs)),
    md5 = unname(tools::md5sum(inputs))
  ),
  file.path(OUT, "provenance", "603_inputs.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUT, "provenance", "603_sessionInfo.txt")
)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg)) {
  source <- sub("^--file=", "", script_arg[1])
  if (file.exists(source)) {
    invisible(file.copy(source, file.path(OUT, "provenance")))
  }
}

writeLines(c(
  "Partial dbRDA of functional gene-potential composition.",
  "Model: Aitchison ~ restoration4 + Condition(locality + depth_f).",
  "51 samples, 17 complete profiles; all restoration stages retained.",
  "Primary pseudocount 1; sensitivity 0.5; unchanged feature filters.",
  "Permutations reused from 602: whole profiles within locality.",
  "Plotted points: response-based site scores, scaling 1, not fitted LC scores.",
  "Axis labels: percentage of constrained inertia, not total variation.",
  "Ellipses: descriptive normal ellipses at 68%; not tests of separation.",
  "Adjusted R2: conventional vegan sample-level adjustment, descriptive;",
  "profile dependence is addressed by the restricted permutation test.",
  "Previous dispersion tests are reported; conditioning does not establish",
  "homogeneous dispersion or causal restoration effects."
), file.path(OUT, "603_analysis_notes.txt"))

write_tsv(
  data.frame(
    status = "PASS", samples = 51L, profiles = 17L,
    models = nrow(result_tbl),
    permutations = nrow(perms),
    primary_pseudocount = 1,
    sensitivity_pseudocount = 0.5
  ),
  "603_run_summary.tsv"
)
writeLines(OUT, file.path(BASE, paste0("LATEST_", SCRIPT_ID, ".txt")))

cat("\nPASS\n")
cat("Results:", file.path(OUT, "tables", "603_partial_dbRDA_results.tsv"), "\n")
cat("Figures:", file.path(OUT, "figures"), "\n")
