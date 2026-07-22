source("Art_Owen_Fast.R")
source("integrands.R")
source("rqmc_histogram.R")

k_values <- c(8, 10, 12)
s_values <- c(2, 4, 8, 16)

R <- 10000
M <- 30

f <- mc2
func_name <- "mc2"

for (s in s_values) {
  for (k in k_values) {
    
    estimates <- numeric(R)
    start_time <- proc.time()
    
    for (r in 1:R) {
      x <- rsobol(m = k, s = s, M = M, seed = r)
      estimates[r] <- mean(f(x))
    }
    
    time_used <- proc.time() - start_time
    cup_time = time_used[["user.self"]] + time_used[["sys.self"]]
    
    cat(
      func_name,
      "  s: ",
      s, "  ",
      "k: ",
      k, "  ",
      "m: ",
      R, "\n\n"
    )
    
    cat("cpu_Time:", cup_time , "seconds\n")
    cat("RQMC estimate:", mean(estimates), "\n")
    cat("RQMC var:", var(estimates), "\n\n")
  }
}



##################hist

# plot_rqmc_histogram(
#   estimates = estimates,
#   modelname = "mc2",
#   s = s,
#   m=m,
#   R = R,
#   elapsed_time = elapsed_time[["elapsed"]] 
# )