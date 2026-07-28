s_values <- c(2, 4, 8)
k_values <- c(8, 10, 12)

regression_results <- data.frame(
  s = integer(),
  k = integer(),
  a = numeric(),
  b = numeric(),
  r_squared = numeric()
)

for (s in s_values) {
  for (k in k_values) {

    input_file <- paste0(
      "R-results/owen_time_by_rep_s",
      s,
      "_k",
      k,
      ".csv"
    )

    data <- read.csv(input_file)

    regression <- lm(median_cpu_time ~ m, data = data)

    regression_results <- rbind(
    regression_results,
    data.frame(
        s = s,
        k = k,
        a = coef(regression)[["m"]],
        b = coef(regression)[["(Intercept)"]],
        r_squared = summary(regression)$r.squared
    )
    )
  }
}

write.csv(
  regression_results,
  "../results/owen_time_regression_results.csv",
  row.names = FALSE
)

print(regression_results)