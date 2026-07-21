source("Art_Owen_Rcode.R")
m <- 4
s <- 4
R <- 1000

f <- function(x) {
  prod(x)
}

estimates <- numeric(R)
start_time <- proc.time()

for (r in 1:R) {
  x <- rsobol(m = m, s = s, seed = r)
  # One RQMC estimate from this randomized point set.
  function_values <- apply(x, 1, f)
  estimates[r] <- mean(function_values)

}

elapsed_time <- proc.time() - start_time


cat("Time needed:", elapsed_time["elapsed"], "seconds\n")

cat("RQMC estimate:", mean(estimates), "\n")
cat("RQMC var:", var(estimates), "\n")