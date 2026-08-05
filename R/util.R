
compute_excess_kurtosis <- function(x, type = 2) {
  if (requireNamespace("e1071", quietly = TRUE)) {
    return(e1071::kurtosis(x, type = type))
  }

  x <- as.numeric(x)
  x <- x[!is.na(x)]
  n <- length(x)

  if (n < 4L) {
    return(NA_real_)
  }

  centered <- x - mean(x)
  m2 <- mean(centered^2)
  m4 <- mean(centered^4)

  if (m2 == 0) {
    return(NA_real_)
  }

  biased_kurtosis <- m4 / (m2^2) - 3

  if (type == 2) {
    return(((n - 1) / ((n - 2) * (n - 3))) * ((n + 1) * biased_kurtosis + 6))
  }

  biased_kurtosis
}


compute_skewness <- function(x, type = 2) {
  if (requireNamespace("e1071", quietly = TRUE)) {
    return(e1071::skewness(x, type = type))
  }

  x <- as.numeric(x)
  x <- x[!is.na(x)]
  n <- length(x)

  if (n < 3L) {
    return(NA_real_)
  }

  centered <- x - mean(x)
  m2 <- mean(centered^2)

  if (m2 == 0) {
    return(NA_real_)
  }

  biased_skewness <- mean(centered^3) / m2^(3 / 2)

  if (type == 2) {
    return(sqrt(n * (n - 1)) / (n - 2) * biased_skewness)
  }

  if (type == 1) {
    return(biased_skewness)
  }

  stop("type must be 1 or 2 when e1071 is unavailable")
}





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
