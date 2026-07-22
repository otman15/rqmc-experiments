source("Art_Owen_Rcode.R")
source("integrands.R")
source("rqmc_histogram.R")

m <- 10
s <- 5
R <- 1000


#f <- function(x) {
#  prod(x)
#}

f=mc2
estimates <- numeric(R)
start_time <- proc.time()

for (r in 1:R) {
  x <- rsobol(m = m, s = s, seed = r)
  # One RQMC estimate from this randomized point set.
  estimates[r] <- mean(f(x))
  #estimates[r] <- mean(function_values)

}

elapsed_time <- proc.time() - start_time

cat(
  "s: ",
  s,"  ",
  "k:",
  m,"  ",
  "m/rep: ",
  R, "\n \n"
)


cat("Time :", elapsed_time["elapsed"], "seconds\n")
cat("RQMC estimate:", mean(estimates), "\n")
cat("RQMC var:", var(estimates), "\n")

##################hist

plot_rqmc_histogram(
  estimates = estimates,
  modelname = "SumUeU",
  s = s,
  m=m,
  R = R,
  elapsed_time = elapsed_time[["elapsed"]]
)