source("integrands.R")

s_values <- c(2, 4, 8, 16)
k_values <- c(8, 10, 12, 14) # n = 2^k

M <- 32
first_seed <- 1
integrand <- mc2

R_values <- c(1, 5, 10, 50, 100)
num_runs <- 10 #number runs for the median time

results <- data.frame(
  method = character(),
  s = integer(),
  k = integer(),
  m = numeric(),
  median_time = numeric()
)

for (s in s_values) {
  for (k in k_values) {

    for (i in seq_along(R_values)) {
      R <- R_values[i]
      cpu_times <- numeric(num_runs)

      for (run in seq_len(num_runs)) {
        source("Art_Owen_Fast.R")

        estimates <- numeric(R)

        start_time <- proc.time()

        for (r in seq_len(R)) {
          points <- rsobol(
            m = k,
            s = s,
            M = M,
            seed = first_seed + r - 1
          )

          estimates[r] <- mean(integrand(points))
        }

        time_used <- proc.time() - start_time

        cpu_times[run] <-
          time_used[["user.self"]] +
          time_used[["sys.self"]]
      }

      results <- rbind(
        results,
        data.frame(
          method = "ART-Owen",
          s = s,
          k = k,
          m = R,
          median_time = median(cpu_times)
        )
      )
    }

    print(paste0("s: ", s, "  k: ", k));
  }
}

write.csv(
  results,
  "R-results/owen_time_by_rep.csv",
  row.names = FALSE,
  quote = FALSE
)

print(results)
