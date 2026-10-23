# Preserve the serial TITAN2 2.4.4 RNG stream when scheduling global gradients
# independently. This advances only the documented RNG calls, not the statistics.
# Every completed fit must have the expected final RNG state.
advance_titan_rng <- function(env, minSplt = 5L, numPerm = 999L, nBoot = 500L) {
  stopifnot(as.character(packageVersion("TITAN2")) == "2.4.4",
            identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection")))
  n <- length(env)
  partitions <- n - 2L * minSplt + 1L
  rank(env, ties.method = "random")
  advance_permutations <- function() {
    for (j in seq_len(partitions)) for (k in seq_len(numPerm)) sample.int(n)
    invisible(NULL)
  }
  advance_permutations()
  for (i in seq_len(nBoot)) {
    runif(1, 0, 100)
    idx <- sample(seq_len(n), replace = TRUE)
    rank(env[idx], ties.method = "random")
    advance_permutations()
  }
  invisible(.Random.seed)
}
