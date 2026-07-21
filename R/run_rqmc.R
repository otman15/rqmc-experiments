source("rsobol_precomputed.R")


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

config <- list(
  fn = "fiftysobol.col",
  s = 4,
  k = 8,
  M = 32,
  replications = 1000,
  first_seed = 1
)


# ------------------------------------------------------------
# Integrand
#
# Replace the body of this function with the function used
# in your experiment.
# ------------------------------------------------------------

f <- function(x) {
  prod(x)
}


# ------------------------------------------------------------
# Seed-independent construction
#
# This section is executed once, before the timed loop.
# ------------------------------------------------------------

# Construct all original, unrandomized Sobol points in bit form.
thebits <- .rsobol.sobolbits(
  fn = config$fn,
  m = config$k,
  s = config$s,
  M = config$M
)

# Construct all prefix indices used by the NUS permutations.
prefix_indices <- .rsobol.prefixindices(
  thebits = thebits,
  m = config$k,
  s = config$s
)


# ------------------------------------------------------------
# Timed RQMC replications
# ------------------------------------------------------------

estimates <- numeric(config$replications)

start_time <- proc.time()[["elapsed"]]

for (r in seq_len(config$replications)) {
  seed <- config$first_seed + r - 1
  
  # Generate one independent NUS randomization.
  x <- rsobol.precomputed(
    thebits = thebits,
    prefix_indices = prefix_indices,
    m = config$k,
    s = config$s,
    M = config$M,
    seed = seed
  )
  
  # Compute one RQMC estimate.
  estimates[r] <- mean(apply(x, 1, f))
}

elapsed_time <-
  proc.time()[["elapsed"]] - start_time


# ------------------------------------------------------------
# Results
# ------------------------------------------------------------
cat(
  "s",
  config$s,"\n",
  "k:",
  config$k,"\n"
)

cat(
  "m/rep",
  config$replications, "\n",
  "time:",
  elapsed_time,
  "seconds\n"
)

cat(
  "RQMC estimate:",
  mean(estimates),
  "\n",
  "RQMC var:",
  var(estimates),
  "\n"
)