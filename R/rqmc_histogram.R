plot_rqmc_histogram <- function(estimates, modelname, s, m, R, elapsed_time) {
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
      modelname,
      " — s = ", s,
      ", n = 2^", m,
      ", m = ", R,
      ", time = ", sprintf("%.1f", elapsed_time), " s"

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
}