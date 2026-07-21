source("rsobol_precomputed.R")
source("integrands.R")


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

config <- list(
  fn = "fiftysobol.col",
  s = 2,
  k = 12,
  M = 30,
  replications = 10000,
  first_seed = 1,
  integrand = sumueu
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



#########################hist

mean_estimate <- mean(estimates)
centered <- estimates - mean_estimate

variance <- var(estimates)

m2 <- mean(centered^2)
skewness <- mean(centered^3) / m2^(3 / 2)
kurtosis <- mean(centered^4) / m2^2 - 3

hist(
  estimates,
  breaks = seq(
    min(estimates),
    max(estimates),
    length.out = 101
  ),
  main = paste0(
    "SumUeU — s = ", config$s,
    ", n = 2^", config$m,
    ", R = ", config$replications
  ),
  xlab = "RQMC estimate",
  ylab = "Frequency"
)

legend(
  "topright",
  legend = c(
    sprintf("Variance = %.6e", variance),
    sprintf("Skewness = %.6f", skewness),
    sprintf("Kurtosis = %.6f", kurtosis)
  ),
  bty = "n"
)