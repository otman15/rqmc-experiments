source("rsobol_precomputed.R")
source("integrands.R")
source("rqmc_histogram.R")

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

config <- list(
  fn = "fiftysobol.col",
  s = 4,
  k = 12,
  M = 32,
  replications = 1000,
  first_seed = 1,
  integrand = mc2
)

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
  estimates[r] <- mean(config$integrand(x))
}

elapsed_time <-
  proc.time()[["elapsed"]] - start_time


# ------------------------------------------------------------
# Results
# ------------------------------------------------------------
cat(
  "s: ",
  config$s,"  ",
  "k:",
  config$k,"  ",
  "m/rep: ",
  config$replications, "\n \n"
)

cat(
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



#########################hist
# 
plot_rqmc_histogram(
  estimates = estimates,
  modelname = "mc2",
  s = config$s,
  m = config$k,
  R = config$replications
)