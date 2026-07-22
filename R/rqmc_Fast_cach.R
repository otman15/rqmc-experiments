source("Art_Owen_Fast_Cache.R")
source("integrands.R")
source("rqmc_histogram.R")

config <- list(
  m = 12,
  s = 4,
  M = 30,
  replications = 10000,
  first_seed = 1,
  modelname = "mc2",
  integrand = mc2
)

estimates <- numeric(config$replications)
start_time <- proc.time()[["elapsed"]]

for (r in seq_len(config$replications)) {
  seed <- config$first_seed + r - 1
  
  x <- rsobol(
    m = config$m,
    s = config$s,
    M = config$M,
    seed = seed
  )
  
  estimates[r] <- mean(config$integrand(x))
}

elapsed_time <- proc.time()[["elapsed"]] - start_time

cat(
  "s: ",
  config$s,"  ",
  "k:",
  config$m,"  ",
  "m/rep: ",
  config$R, "\n \n"
)

cat("Time needed:", elapsed_time, "seconds\n")
cat("RQMC estimate:", mean(estimates), "\n")
cat("RQMC var:", var(estimates), "\n")

plot_rqmc_histogram(
  estimates = estimates,
  modelname = config$modelname,
  s = config$s,
  m = config$m,
  R = config$replications
)