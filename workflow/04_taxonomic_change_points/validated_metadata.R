# Validate environmental values and sample design against the canonical 003b file.
# CSV byte layout is not a scientific identifier: compare values by sample ID.
titan_script_dir <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a) != 1L) stop("Run this script with Rscript.")
  dirname(normalizePath(sub("^--file=", "", a), mustWork = TRUE))
}

canonical_titan_metadata <- function(script_dir) {
  normalizePath(file.path(script_dir, "../../source_data/indices/metadata_003b.csv"),
                mustWork = TRUE)
}

validate_titan_metadata <- function(meta, canonical_path, audit_path = NULL,
                                    sample_col = "sample_id") {
  ref <- read.csv(canonical_path, check.names = FALSE, stringsAsFactors = FALSE)
  gradients <- c("MHI_local", "moisture_stress_PC1", "nutrients_redox_PC1",
                 "physicochemical_PC1", "vegetation_landscape_PC1", "water_inundation_PC1")
  design <- c("locality", "depth_cm", "restoration4", "lat_block")
  required <- c(sample_col, gradients, design)
  if (!all(required %in% names(meta))) stop("Missing canonical metadata columns: ",
                                          paste(setdiff(required, names(meta)), collapse = ", "))
  ids <- as.character(meta[[sample_col]])
  if (length(ids) != 51L || anyNA(ids) || anyDuplicated(ids) ||
      !setequal(ids, ref$sample_id)) stop("Exactly the 51 canonical sample IDs are required.")
  ref <- ref[match(ids, ref$sample_id), , drop = FALSE]
  for (v in design) {
    if (!identical(as.character(meta[[v]]), as.character(ref[[v]])))
      stop("Sample design differs from canonical 003b: ", v)
  }
  audit <- do.call(rbind, lapply(gradients, function(v) {
    observed <- as.numeric(meta[[v]])
    delta <- observed - ref[[v]]
    data.frame(sample_id = ids, gradient = v, used = observed,
               canonical = ref[[v]], difference = delta,
               pass = is.finite(delta) & abs(delta) <= 1e-12)
  }))
  if (!is.null(audit_path)) write.csv(audit, audit_path, row.names = FALSE)
  if (!all(audit$pass)) stop("Environmental values differ from canonical 003b. ",
                            "Use the validated metadata; do not rename historical indices.")
  invisible(audit)
}
