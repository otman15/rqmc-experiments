source("integrands.R")
source("util.R")


run_rqmc_replications <- function(integrand, k, s, m, M, first_seed) {
  source("Art_Owen_Fast.R")


  start_time <- proc.time()
  estimates <- numeric(m)

  for (replication in seq_len(m)) {
    points <- rsobol(
      m = k,
      s = s,
      M = M,
      seed = first_seed + replication - 1
    )

    estimates[replication] <- mean(integrand(points))
  }

  time_used <- proc.time() - start_time
  cpu_time <- time_used[["user.self"]] + time_used[["sys.self"]]

  list(estimates = estimates, cpu_time = cpu_time)
}


run_timing_experiments <- function(
  models,
  s_values,
  k_values,
  timing_m_values,
  M,
  first_seed,
  num_runs,
  output_dir
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (model in models) {
    results <- data.frame(
      method = character(),
      s = integer(),
      k = integer(),
      m = numeric(),
      median_time = numeric()
    )

    for (s in s_values) {
      for (k in k_values) {
        for (m in timing_m_values) {
          cpu_times <- numeric(num_runs)

          for (run in seq_len(num_runs)) {
            run_first_seed <- first_seed + (run - 1) * m
            batch <- run_rqmc_replications(
              model$integrand,
              k,
              s,
              m,
              M,
              run_first_seed
            )
            cpu_times[run] <- batch$cpu_time
          }

          median_time <- median(cpu_times)
          results <- rbind(
            results,
            data.frame(
              method = "ART-Owen",
              s = s,
              k = k,
              m = m,
              median_time = median_time
            )
          )

          cat(sprintf(
            "model=%s, s=%d, k=%d, m=%d, median CPU time=%.6f seconds\n",
            model$name,
            s,
            k,
            m,
            median_time
          ))
        }
      }
    }

    output_path <- file.path(
      output_dir,
      paste0("owen_", model$name, "_time_rep-", num_runs, ".csv")
    )
    write.csv(results, output_path, row.names = FALSE, quote = FALSE)

    print(results)
    cat("Results written to", output_path, "\n")
  }
}


run_statistics_experiments <- function(
  models,
  method,
  s_values,
  k_values,
  statistics_m,
  M,
  first_seed,
  output_dir
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (model in models) {
    for (s in s_values) {
      for (k in k_values) {
        batch <- run_rqmc_replications(
          model$integrand,
          k,
          s,
          statistics_m,
          M,
          first_seed
        )
        estimates <- batch$estimates
        mean_estimate <- mean(estimates)
        variance <- if (statistics_m > 1) var(estimates) else NA_real_
        skewness <- compute_skewness(estimates, type = 2)
        excess_kurtosis <- compute_excess_kurtosis(estimates, type = 2)

        cat(sprintf(
          paste0(
            "model=%s, s=%d, k=%d, m=%d, mean=%.16e, variance=%.16e, ",
            "skewness=%.16e, excess kurtosis=%.16e, CPU time=%.6f seconds\n"
          ),
          model$name,
          s,
          k,
          statistics_m,
          mean_estimate,
          variance,
          skewness,
          excess_kurtosis,
          batch$cpu_time
        ))

        output_path <- file.path(
          output_dir,
          paste0(model$name, "-", s, "-", method, "-", k, "-", statistics_m, ".dat")
        )
        writeLines(
          formatC(sort(estimates), format = "e", digits = 17, decimal.mark = "."),
          output_path
        )
        cat("Results written to", output_path, "\n")
      }
    }
  }
}


main <- function() {
  TIMING <- "TIMING"
  STATISTICS <- "STATISTICS"
  mode <- TIMING

  method <- "ART-Owen"
  models <- list(
    list(name = "MC2", integrand = mc2),
    list(name = "SumUeU", integrand = sumueu),
    list(
      name = "SmoothPerB4",
      integrand = function(x) smoothperb4(x, expw = 1.0)
    ),
    list(name = "Polynomial", integrand = polynomial)
  )

  models_Timing <- list(list(name = "MC2", integrand = mc2))
  
  s_values <- c(2, 4, 8, 16, 32)
  k_values <- c(8, 10, 12, 14, 16)
  timing_m_values <- c(1, 2, 5, 10, 20, 30, 50)
  statistics_m <- 1000
  M <- 32
  first_seed <- 1
  num_runs <- 11

  timing_output_dir <- "R-results"
  statistics_output_dir <- file.path(timing_output_dir, "dat")

  if (mode == TIMING) {
    run_timing_experiments(
      models_Timing,
      s_values,
      k_values,
      timing_m_values,
      M,
      first_seed,
      num_runs,
      timing_output_dir
    )
  } else if (mode == STATISTICS) {
    run_statistics_experiments(
      models,
      method,
      s_values,
      k_values,
      statistics_m,
      M,
      first_seed,
      statistics_output_dir
    )
  } else {
    stop("mode must be TIMING or STATISTICS")
  }
}


main()
