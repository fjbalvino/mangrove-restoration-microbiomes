#!/usr/bin/env Rscript
# ============================================================
# 06_02_MHI_mixed_models.R
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Corrected 003b MHI, REML/lmerTest, profile random intercept; separate MHI multiplicity families.
# Inputs (source expressions; complete list in docs/contracts/06_02_MHI_mixed_models.json):
#   read_tsv <- function(p) read.delim(p, check.names = FALSE)
#   if (!nzchar(source_run)) source_run <- trimws(readLines(pointer, n = 1L))
#   s <- read_tsv(paths[1])
#   m <- read.csv(paths[2], check.names = FALSE)
#   selection <- read_tsv(paths[4])
#   z <- tryCatch(read.delim(con, check.names = FALSE), finally = close(con))
# Outputs (source expressions; complete list in contract):
#   write_tsv <- function(x, p) data.table::fwrite(x, p, sep = "\t", na = "NA")
#   writeLines(capture.output(sessionInfo()), file.path(run, "logs", "603b_sessionInfo.txt"))
#   write_tsv(data.frame(input = paths, md5 = unname(hashes)), file.path(run, "tables", "603b_inputs.tsv"))
#   write_tsv(m, file.path(run, "tables", "603b_metadata_used.tsv"))
#   write_tsv(data.frame(MHI_mean = mhi_mean, MHI_sd = mhi_sd, standardization_n = 51),
#   saveRDS(part, file.path(run, "checkpoints", sprintf("batch_%04d.rds", h)))
#   write_tsv(r, file.path(run, "tables", "603b_all_results.tsv"))
#   write_tsv(diagnostics, file.path(run, "tables", "603b_fit_diagnostics.tsv"))
#   write_tsv(families, file.path(run, "tables", "603b_test_families.tsv"))
#   write_tsv(data.frame(status = status, mode = mode, samples = 51, profiles = 17,
# Algorithmic provenance:
# REML mixed models for MHI, locality, depth and profile intercept; Satterthwaite and BH.
#   Bates et al. (2015), doi:10.18637/jss.v067.i01; Kuznetsova et al. (2017), doi:10.18637/jss.v082.i13; Benjamini & Hochberg (1995), doi:10.1111/j.2517-6161.1995.tb02031.x.
# Source SHA-256: 3bbf650e53d676f6fd3a2ca0d9e32394fed4bb0e15786a2fcaf677d48cad567a
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================

# Figure 4: exploratory protein associations; 51 samples, 17 profiles.
# Does not average depths, refit MHI, re-centre CLR, or run enrichment.
options(stringsAsFactors = FALSE, warn = 1)
NAME <- "603b_asociaciones_MHI_por_profundidad"
BASE <- "/home/fjbalvino/Tipping_points/resultados_finales"
argv <- commandArgs(TRUE)
stopifnot(length(argv) %% 2L == 0L)
if (length(argv)) stopifnot(all(argv[seq.int(1L, length(argv), by = 2L)] %in%
  c("--mode", "--out_root", "--workers", "--input_run")))
arg <- function(key, default) {
  i <- which(argv == key)
  if (!length(i)) return(default)
  stopifnot(length(i) == 1L, i < length(argv))
  argv[i + 1L]
}
mode <- arg("--mode", "pilot")
root <- arg("--out_root", BASE)
workers <- as.integer(arg("--workers", "2"))
stopifnot(mode %in% c("pilot", "full"), is.finite(workers), workers >= 1L)
for (p in c("lme4", "lmerTest", "data.table", "emmeans")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Falta paquete: ", p)
}
data.table::setDTthreads(1L)
read_tsv <- function(p) read.delim(p, check.names = FALSE)
write_tsv <- function(x, p) data.table::fwrite(x, p, sep = "\t", na = "NA")
need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
flag <- function(x, expected) identical(tolower(as.character(x)), expected)
pointer <- file.path(root, "LATEST_601b_preparar_CLR_proteico_por_muestra.txt")
source_run <- arg("--input_run", "")
if (!nzchar(source_run)) source_run <- trimws(readLines(pointer, n = 1L))
source_run <- normalizePath(source_run, mustWork = TRUE)
paths <- file.path(source_run, "tables", c(
  "601b_run_summary.tsv", "601b_metadata_51samples.csv",
  "601b_top_100000_CLR_features_x_samples.tsv.gz",
  "601b_selected_features_top_100000.tsv"))
need(all(file.exists(paths)), "Faltan salidas de 601b")
hashes <- tools::md5sum(paths)
s <- read_tsv(paths[1])
need(nrow(s) == 1L && s$status == "PASS" && s$samples == 51L &&
       s$profiles == 17L && s$selected_features == 100000L &&
       s$catalog_features == 1475487L && s$pseudocount == 1 &&
       flag(s$depth_averaged, "false") &&
       flag(s$selection_uses_HI_MHI_HFR_stage, "false") &&
       s$CLR_center == "full_filtered_catalog" &&
       s$selection == "global_CLR_variance", "601b no coincide con el PASS esperado")
m <- read.csv(paths[2], check.names = FALSE)
required <- c("sample_id", "profile_id", "locality", "restoration4", "depth_cm", "MHI_local")
need(all(required %in% names(m)), "Faltan campos de metadata")
need(nrow(m) == 51L && !anyNA(m[, required]) && !anyDuplicated(m$sample_id),
     "Metadata incompleta o duplicada")
for (k in c("sample_id", "profile_id", "locality", "restoration4")) {
  need(all(nzchar(trimws(m[[k]]))), paste("Campo vacio:", k))
}
m$depth_cm <- as.numeric(m$depth_cm)
m$MHI_local <- as.numeric(m$MHI_local)
need(all(is.finite(m$MHI_local)) && all(abs(m$MHI_local) <= 1), "MHI invalido")
profiles <- split(seq_len(nrow(m)), m$profile_id)
need(length(profiles) == 17L && all(vapply(profiles, function(i) {
  identical(sort(m$depth_cm[i]), c(5, 20, 40)) &&
    length(unique(m$locality[i])) == 1L &&
    length(unique(m$restoration4[i])) == 1L
}, logical(1))), "Se requieren 17 perfiles completos y consistentes")
need(setequal(m$locality, c("Carmen", "Cozumel", "Tuxpan")), "Localidades incompatibles")
m$depth <- factor(m$depth_cm, levels = c(5, 20, 40))
m$locality <- factor(m$locality, levels = c("Carmen", "Cozumel", "Tuxpan"))
m$profile_id <- factor(m$profile_id)
contrasts(m$depth) <- contr.treatment(3L)
contrasts(m$locality) <- contr.treatment(3L)
mhi_mean <- mean(m$MHI_local)
mhi_sd <- sd(m$MHI_local)
need(is.finite(mhi_sd) && mhi_sd > 0, "MHI sin variacion")
m$MHI_z <- (m$MHI_local - mhi_mean) / mhi_sd
forms <- list(common = y ~ MHI_z + depth + locality + (1 | profile_id),
              interaction = y ~ MHI_z * depth + locality + (1 | profile_id))
for (f in forms) {
  X <- model.matrix(lme4::nobars(f), transform(m, y = 0))
  need(qr(X)$rank == ncol(X), "Modelo con rango deficiente")
}
selection <- read_tsv(paths[4])
need(nrow(selection) == 100000L && !anyDuplicated(selection$feature_id) &&
       all(selection$rank == seq_len(100000L)) &&
       all(is.finite(selection$sd_clr) & selection$sd_clr > 0), "Seleccion invalida")
cat("Leyendo CLR preparado por 601b; no se recalcula ni se promedian profundidades.\n")
con <- gzfile(paths[3], "rt")
z <- tryCatch(read.delim(con, check.names = FALSE), finally = close(con))
need(nrow(z) == 100000L && identical(names(z)[1], "feature_id") &&
       length(names(z)) == 52L && !anyDuplicated(names(z)) &&
       setequal(names(z)[-1], m$sample_id) &&
       identical(as.character(z$feature_id), as.character(selection$feature_id)),
     "IDs o dimensiones CLR incompatibles")
Y <- as.matrix(z[, m$sample_id, drop = FALSE])
need(is.numeric(Y) && all(is.finite(Y)), "CLR no numerico/no finito")
rm(z)
need(max(abs(apply(Y, 1L, sd) - selection$sd_clr)) < 1e-7, "SD del CLR discrepante")

# Contrast vectors built from the actual design, not coefficient-position guesses.
design <- function(depth_value, x) {
  d <- m[1L, , drop = FALSE]
  d$depth <- factor(depth_value, levels = levels(m$depth))
  contrasts(d$depth) <- contrasts(m$depth)
  d$MHI_z <- x
  model.matrix(~ MHI_z * depth + locality, d)
}
slopes <- do.call(rbind, lapply(c("5", "20", "40"), function(d) design(d, 1) - design(d, 0)))
interaction_L <- rbind(slopes[2L, ] - slopes[1L, ], slopes[3L, ] - slopes[1L, ])
capture <- function(expr) {
  notes <- character()
  ans <- tryCatch(withCallingHandlers(expr,
    warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleMessage") }),
    error = function(e) e)
  list(value = ans, notes = unique(notes))
}
fit_one <- function(y, which_model) {
  d <- m
  d$y <- y
  a <- capture(lmerTest::lmer(forms[[which_model]], data = d, REML = TRUE,
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))
  if (inherits(a$value, "error")) return(list(fit = NULL, status = "error", singular = NA,
    notes = conditionMessage(a$value), icc = NA_real_))
  f <- a$value
  singular <- lme4::isSingular(f, tol = 1e-4)
  msgs <- unique(c(a$notes, unlist(f@optinfo$conv$lme4$messages)))
  nonboundary <- msgs[!grepl("boundary.*singular", msgs, ignore.case = TRUE)]
  ok <- all(f@optinfo$conv$opt == 0) && !length(nonboundary) &&
    all(is.finite(lme4::fixef(f)))
  v <- as.data.frame(lme4::VarCorr(f))$vcov
  list(fit = f, status = if (!ok) "review_warning" else if (singular) "singular" else "ok",
       singular = singular, notes = paste(msgs, collapse = " | "), icc = v[1L] / sum(v))
}
contrast_row <- function(a, L, family, depth = NA_integer_, joint = FALSE) {
  row <- data.frame(family = family, depth_cm = depth, estimate = NA_real_, se = NA_real_,
    df_num = if (joint) 2 else 1, df_den = NA_real_, statistic = NA_real_,
    ci_low = NA_real_, ci_high = NA_real_, p = NA_real_, status = a$status,
    singular = a$singular, icc = a$icc, notes = a$notes)
  if (is.null(a$fit)) return(row)
  a_test <- capture(if (joint) lmerTest::contestMD(a$fit, L, ddf = "Satterthwaite") else
    lmerTest::contest1D(a$fit, L, ddf = "Satterthwaite", confint = TRUE, level = 0.95))
  t <- a_test$value
  if (inherits(t, "error")) {
    row$status <- "test_error"
    row$notes <- paste(row$notes, conditionMessage(t), sep = " | ")
    return(row)
  }
  if (joint) {
    row[, c("df_num", "df_den", "statistic", "p")] <-
      t[, c("NumDF", "DenDF", "F value", "Pr(>F)"), drop = FALSE]
  } else {
    row[, c("estimate", "se", "df_den", "statistic", "ci_low", "ci_high", "p")] <-
      t[, c("Estimate", "Std. Error", "df", "t value", "lower", "upper", "Pr(>|t|)"), drop = FALSE]
  }
  if (length(a_test$notes)) {
    row$status <- "test_warning"
    row$notes <- paste(c(row$notes, a_test$notes), collapse = " | ")
  }
  values <- if (joint) unlist(row[c("statistic", "df_den", "p")]) else
    unlist(row[c("estimate", "se", "df_den", "statistic", "ci_low", "ci_high", "p")])
  if (any(!is.finite(values)) || row$p < 0 || row$p > 1 || row$df_den <= 0) row$status <- "invalid_test"
  row
}
analyse <- function(i) {
  common <- fit_one(Y[i, ], "common")
  k <- ncol(model.matrix(~ MHI_z + depth + locality, m))
  L <- numeric(k)
  L[match("MHI_z", colnames(model.matrix(~ MHI_z + depth + locality, m)))] <- 1
  interaction <- fit_one(Y[i, ], "interaction")
  if (!is.null(interaction$fit)) need(identical(names(lme4::fixef(interaction$fit)),
                                              colnames(slopes)), "Orden de contrastes incompatible")
  rows <- list(contrast_row(common, L, "MHI_common"),
               contrast_row(interaction, interaction_L, "MHI_by_depth", joint = TRUE))
  for (j in 1:3) rows[[length(rows) + 1L]] <-
    contrast_row(interaction, slopes[j, ], "MHI_depth_slopes", c(5L, 20L, 40L)[j])
  ans <- data.table::rbindlist(rows)
  ans$feature_id <- selection$feature_id[i]
  ans$selection_rank <- selection$rank[i]
  ans$estimate_sd <- ans$estimate / selection$sd_clr[i]
  ans$ci_low_sd <- ans$ci_low / selection$sd_clr[i]
  ans$ci_high_sd <- ans$ci_high / selection$sd_clr[i]
  ans
}

# Independent contrast check with emmeans before any production fitting.
set.seed(603)
u <- rnorm(17, sd = 3)
test_y <- 1.2 * m$MHI_z + 0.8 * m$MHI_z * (m$depth_cm == 20) -
  0.5 * m$MHI_z * (m$depth_cm == 40) + u[as.integer(m$profile_id)] + rnorm(51, sd = 0.5)
tf <- fit_one(test_y, "interaction")
need(!is.null(tf$fit) && tf$status == "ok", "Fallo en prueba sintetica del modelo")
et <- as.data.frame(emmeans::emtrends(tf$fit, ~ depth, var = "MHI_z",
  data = transform(m, y = test_y), lmer.df = "satterthwaite"))
need(identical(as.character(et$depth), c("5", "20", "40")), "Orden inesperado emmeans")
b <- lme4::fixef(tf$fit)
V <- as.matrix(vcov(tf$fit))
need(max(abs(drop(slopes %*% b) - et$MHI_z.trend)) < 1e-7 &&
       max(abs(sqrt(diag(slopes %*% V %*% t(slopes))) - et$SE)) < 1e-7,
     "Contrastes no coinciden con emmeans")
for (j in 1:3) need(contrast_row(tf, slopes[j, ], "selftest")$status == "ok", "Fallo contest1D")
need(contrast_row(tf, interaction_L, "selftest", joint = TRUE)$status == "ok", "Fallo contestMD")
cat("Autocomprobacion de modelos y contrastes: PASS\n")
idx <- if (mode == "pilot") unique(as.integer(round(seq(1, nrow(Y), length.out = 500)))) else seq_len(nrow(Y))
run <- file.path(root, paste0(NAME, "_", mode, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
need(!dir.exists(run), "Carpeta de salida ya existe")
dir.create(file.path(run, "tables"), recursive = TRUE)
dir.create(file.path(run, "checkpoints"))
dir.create(file.path(run, "logs"))
writeLines(capture.output(sessionInfo()), file.path(run, "logs", "603b_sessionInfo.txt"))
self <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(self)) file.copy(sub("^--file=", "", self[1]), file.path(run, "logs", paste0(NAME, ".R")))
write_tsv(data.frame(input = paths, md5 = unname(hashes)), file.path(run, "tables", "603b_inputs.tsv"))
write_tsv(m, file.path(run, "tables", "603b_metadata_used.tsv"))
write_tsv(data.frame(MHI_mean = mhi_mean, MHI_sd = mhi_sd, standardization_n = 51),
          file.path(run, "tables", "603b_scaling.tsv"))
cat("Output:", run, "\nMode:", mode, "| Features:", length(idx), "| Workers:", workers, "\n")
started <- proc.time()[3]
chunks <- split(idx, ceiling(seq_along(idx) / 100L))
all_results <- vector("list", length(chunks))
for (h in seq_along(chunks)) {
  ii <- chunks[[h]]
  out <- if (workers == 1L) lapply(ii, analyse) else
    parallel::mclapply(ii, analyse, mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
  need(all(vapply(out, is.data.frame, logical(1))), "Fallo de proceso; revisar log y checkpoints")
  part <- data.table::rbindlist(out)
  saveRDS(part, file.path(run, "checkpoints", sprintf("batch_%04d.rds", h)))
  all_results[[h]] <- part
  done <- min(h * 100L, length(idx))
  elapsed <- proc.time()[3] - started
  cat(format(Sys.time(), "%H:%M:%S"), "|", done, "/", length(idx),
      "proteinas | min:", round(elapsed / 60, 1),
      "| estimacion 100000:", round(elapsed / done * 100000 / 3600, 2), "h\n")
  flush.console()
}
r <- data.table::rbindlist(all_results)
need(nrow(r) == 5L * length(idx), "Numero de resultados inesperado")
# Do not select features by pilot significance. FDR only over complete families.
r$usable <- r$status == "ok" & is.finite(r$p)
r$q <- NA_real_
if (mode == "full") {
  for (fam in unique(r$family)) {
    jj <- which(r$family == fam)
    pp <- ifelse(r$usable[jj], r$p[jj], 1)
    qq <- p.adjust(pp, method = "BH")
    r$q[jj] <- ifelse(r$usable[jj], qq, NA_real_)
  }
}
r$direction <- ifelse(is.na(r$estimate), NA_character_,
                      ifelse(r$estimate > 0, "higher_MHI", "lower_MHI"))
r$passes_q_0_10 <- !is.na(r$q) & r$q <= 0.10
diagnostics <- r[, .(n_tests = .N), by = .(family, depth_cm, status, singular)]
families <- r[, .(n_tests = .N, n_usable = sum(usable),
                  n_q_0_10 = if (mode == "full") sum(passes_q_0_10) else NA_integer_), by = family]
write_tsv(r, file.path(run, "tables", "603b_all_results.tsv"))
write_tsv(diagnostics, file.path(run, "tables", "603b_fit_diagnostics.tsv"))
write_tsv(families, file.path(run, "tables", "603b_test_families.tsv"))
print(diagnostics)
need(identical(unname(hashes), unname(tools::md5sum(paths))), "Entradas cambiaron durante ejecucion")
need(all(families$n_usable > 0), "Ninguna prueba utilizable en al menos una familia; revisar")
status <- if (mode == "pilot") "PILOT_COMPLETE" else "PASS"
write_tsv(data.frame(status = status, mode = mode, samples = 51, profiles = 17,
  features = length(idx), catalog_features = 1475487, pseudocount = 1,
  depth_averaged = FALSE, MHI_refitted = FALSE, fit_method = "REML",
  test_method = "Satterthwaite", random_effect = "profile_intercept",
  FDR_computed = mode == "full", q_cutoff = 0.10,
  non_ok_tests = sum(!r$usable), elapsed_minutes = (proc.time()[3] - started) / 60),
  file.path(run, "tables", "603b_run_summary.tsv"))
suffix <- if (mode == "pilot") "_pilot" else ""
latest <- file.path(root, paste0("LATEST_", NAME, suffix, ".txt"))
tmp <- paste0(latest, ".", Sys.getpid(), ".tmp")
writeLines(run, tmp)
need(file.rename(tmp, latest), "No se pudo publicar puntero")
cat("\n", status, "\nPointer:", latest, "\n", sep = "")
if (mode == "pilot") cat("Piloto computacional: sin FDR, sin seleccion de hits; revisar diagnosticos.\n")
