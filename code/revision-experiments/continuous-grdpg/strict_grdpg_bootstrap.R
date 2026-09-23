#!/usr/bin/env Rscript

# Strict-GRDPG extension of the Section 5.2 coverage/power experiment.
#
# The continuous latent-position model is
#   xi_i ~ Uniform(0, 1),
#   Z_i = (0.69 + 0.05 cos(2 pi xi_i),
#          0.45 sin(2 pi xi_i),
#          0.40 cos(4 pi xi_i)),
#   P_ij = rho_n Z_i' diag(1, 1, -1) Z_j.
#
# Its signed population eigenvalues are 0.47735, 0.10125, and -0.08.
# Hence it is a continuous, genuinely indefinite GRDPG rather than an SBM
# or an ordinary RDPG. Its mean kernel value, 0.4761, nearly matches the
# value 0.47625 in the three-block SBM used in the submitted manuscript.

suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppArmadillo)
  library(RSpectra)
  library(parallel)
})

parse_arguments <- function(args) {
  out <- list(
    mode = "pilot",
    n = NULL,
    nsim = NULL,
    nboot = NULL,
    cores = NULL,
    delta = "-0.3",
    msizes = "none",
    seed = 2601L,
    signals = "null,weak,strong",
    output = file.path("code", "revision-experiments", "continuous-grdpg", "results")
  )

  for (arg in args) {
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(pieces) != 2L || !pieces[1] %in% names(out)) {
      stop("Unknown argument: ", arg)
    }
    key <- pieces[1]
    value <- pieces[2]
    if (key %in% c("mode", "output", "signals", "delta", "msizes")) {
      out[[key]] <- value
    } else if (key %in% c("n", "nsim", "nboot", "cores", "seed")) {
      out[[key]] <- as.integer(value)
    } else {
      out[[key]] <- as.numeric(value)
    }
  }

  if (!out$mode %in% c("pilot", "full")) {
    stop("--mode must be pilot or full")
  }

  defaults <- if (out$mode == "pilot") {
    list(n = 600L, nsim = 20L, nboot = 399L, cores = 4L)
  } else {
    list(n = 3000L, nsim = 400L, nboot = 2000L, cores = 4L)
  }
  for (key in names(defaults)) {
    if (is.null(out[[key]])) out[[key]] <- defaults[[key]]
  }
  out$signals <- strsplit(out$signals, ",", fixed = TRUE)[[1]]
  if (!length(out$signals) || anyDuplicated(out$signals) ||
      !all(out$signals %in% c("null", "weak", "strong"))) {
    stop("--signals must contain distinct values from null,weak,strong")
  }
  for (key in c("n", "nsim", "nboot", "cores", "seed")) {
    if (length(out[[key]]) != 1L || !is.finite(out[[key]]) || out[[key]] < 1L)
      stop("--", key, " must be a positive integer")
  }
  out$delta <- as.numeric(strsplit(out$delta, ",", fixed = TRUE)[[1]])
  if (out$n < 6L || !length(out$delta) || any(!is.finite(out$delta)) ||
      any(out$delta > 0) || anyDuplicated(out$delta))
    stop("Use n >= 6 and distinct finite deltas <= 0")
  out$msizes <- if (out$msizes == "none") integer(0) else
    as.numeric(strsplit(out$msizes, ",", fixed = TRUE)[[1]])
  if (any(!is.finite(out$msizes)) || any(out$msizes != floor(out$msizes)) ||
      any(out$msizes < 6L | out$msizes >= out$n) || anyDuplicated(out$msizes))
    stop("--msizes must be none or distinct integers with 6 <= m < n")
  out$msizes <- as.integer(out$msizes)
  out
}

config <- parse_arguments(commandArgs(trailingOnly = TRUE))
dir.create(config$output, recursive = TRUE, showWarnings = FALSE)
script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- sub("^--file=", "", script_argument[1L])
implementation_fingerprint <- unname(tools::md5sum(script_path))
checkpoint_dir <- file.path(config$output, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
log_connection <- file(file.path(config$output, "run.log"), open = "at")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")
message("Started ", format(Sys.time(), tz = "UTC"), " UTC; PID ", Sys.getpid(),
        "; source fingerprint ", implementation_fingerprint)

save_checkpoint <- function(object, path) {
  temporary <- tempfile(pattern = "checkpoint_", tmpdir = dirname(path), fileext = ".rds")
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) stop("Could not save checkpoint: ", path)
}

Rcpp::cppFunction(
  code = '
  arma::mat simulate_strict_grdpg_cpp(const arma::mat& Z,
                                       const double rho) {
    const int n = Z.n_rows;
    arma::mat A(n, n, arma::fill::zeros);

    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        const double kernel = Z(i, 0) * Z(j, 0)
                            + Z(i, 1) * Z(j, 1)
                            - Z(i, 2) * Z(j, 2);
        const double probability = std::min(1.0, std::max(0.0, rho * kernel));
        A(i, j) = R::rbinom(1.0, probability);
        A(j, i) = A(i, j);
      }
    }
    return A;
  }',
  depends = "RcppArmadillo"
)

latent_positions <- function(n) {
  xi <- runif(n)
  Z <- cbind(
    0.69 + 0.05 * cos(2 * pi * xi),
    0.45 * sin(2 * pi * xi),
    0.40 * cos(4 * pi * xi)
  )
  colnames(Z) <- c("z1", "z2", "z3")
  list(xi = xi, Z = Z)
}

adjacency_spectral_embedding <- function(A, dimension, rho) {
  decomposition <- RSpectra::eigs_sym(A, k = dimension, which = "LM")
  ordering <- order(abs(decomposition$values), decreasing = TRUE)
  values <- decomposition$values[ordering]
  vectors <- decomposition$vectors[, ordering, drop = FALSE]
  embedding <- vectors %*% diag(sqrt(abs(values)), nrow = dimension) / sqrt(rho)
  list(Z = embedding, eigenvalues = values)
}

gaussian_multiplier_draws <- function(influence, number) {
  n <- nrow(influence)
  covariance <- crossprod(influence) / n^2
  decomposition <- eigen(covariance, symmetric = TRUE)
  tolerance <- max(decomposition$values) * 1e-10
  keep <- decomposition$values > tolerance
  if (!any(keep)) {
    return(matrix(0, nrow = number, ncol = ncol(influence)))
  }
  root <- sweep(
    t(decomposition$vectors[, keep, drop = FALSE]),
    1L,
    sqrt(decomposition$values[keep]),
    `*`
  )
  matrix(rnorm(number * sum(keep)), nrow = number) %*% root
}

bootstrap_ols <- function(design, response, number) {
  n <- nrow(design)
  p <- ncol(design)
  quadratic <- crossprod(design) / n
  linear <- colMeans(design * response)
  beta_hat <- solve(quadratic, linear)

  quadratic_influence <- matrix(NA_real_, nrow = n, ncol = p^2)
  for (column in seq_len(p)) {
    for (row in seq_len(p)) {
      index <- row + (column - 1L) * p
      quadratic_influence[, index] <-
        design[, row] * design[, column] - quadratic[row, column]
    }
  }
  linear_influence <- design * response -
    matrix(linear, nrow = n, ncol = p, byrow = TRUE)
  perturbations <- gaussian_multiplier_draws(
    cbind(quadratic_influence, linear_influence),
    number
  )

  beta_boot <- matrix(NA_real_, nrow = p, ncol = number)
  for (b in seq_len(number)) {
    quadratic_b <- quadratic + matrix(perturbations[b, seq_len(p^2)], p, p)
    quadratic_b <- (quadratic_b + t(quadratic_b)) / 2
    linear_b <- linear + perturbations[b, p^2 + seq_len(p)]
    beta_boot[, b] <- tryCatch(
      solve(quadratic_b, linear_b),
      error = function(e) rep(NA_real_, p)
    )
  }

  list(beta_hat = beta_hat, beta_boot = beta_boot)
}

one_repetition <- function(repetition, signal, beta, config) {
  signal_offset <- switch(signal, strong = 100000L, weak = 200000L, null = 300000L)
  set.seed(config$seed + signal_offset + repetition)

  n <- config$n
  rho <- n^config$delta
  latent <- latent_positions(n)
  A <- simulate_strict_grdpg_cpp(latent$Z, rho)
  embedding <- adjacency_spectral_embedding(A, dimension = 3L, rho = rho)

  conventional_covariance <- matrix(c(1, 0.3, 0.3, 1), 2L, 2L)
  conventional <- matrix(rnorm(2L * n), nrow = n, ncol = 2L) %*%
    chol(conventional_covariance)

  latent_design <- cbind(conventional, latent$Z)
  observed_design <- cbind(conventional, embedding$Z)
  response <- drop(latent_design %*% beta) + rnorm(n)

  edge_density <- sum(A) / (n * (n - 1))
  # Principal-angle cosines diagnose recovery without aligning the regression.
  principal_cosines <- svd(crossprod(qr.Q(qr(latent$Z)),
                                     qr.Q(qr(embedding$Z))), nu = 0, nv = 0)$d

  infer_sample <- function(m) {
  # One subset per m per independent graph. Fixed full-graph ASE is reused.
  # A separate seed for each m makes subset selection insensitive to ordering
  # and to the number of full-sample bootstrap draws.
  if (m == n) {
    index <- seq_len(n)
    subset_seed <- NA_integer_
  } else {
    subset_seed <- as.integer((as.double(config$seed) + signal_offset +
      1000003 * repetition + 1009 * m) %% 2147483646 + 1)
    set.seed(subset_seed)
    index <- sample.int(n, m, replace = FALSE)
  }
  design_m <- observed_design[index, , drop = FALSE]
  response_m <- response[index]
  base_row <- data.frame(
    repetition = repetition, signal = signal, n = n, m = m,
    method = if (m == n) "full" else "downsample",
    delta = config$delta, rho = rho, subset_seed = subset_seed,
    edge_density = edge_density, mean_degree = edge_density * (n - 1),
    isolated_fraction = mean(rowSums(A) == 0),
    subspace_min_cosine = min(principal_cosines),
    subspace_mean_squared_cosine = mean(principal_cosines^2),
    eig1 = embedding$eigenvalues[1L], eig2 = embedding$eigenvalues[2L],
    eig3 = embedding$eigenvalues[3L],
    negative_eigenvalues = sum(embedding$eigenvalues < 0),
    design_condition = kappa(crossprod(design_m), exact = TRUE)
  )
  fit <- tryCatch(bootstrap_ols(design_m, response_m, config$nboot),
                  error = function(e) list(error = conditionMessage(e)))
  valid_draws <- if (is.null(fit$error))
    sum(colSums(is.finite(fit$beta_boot)) == length(beta)) else 0L
  if (!is.null(fit$error) || valid_draws < 0.95 * config$nboot) {
    # Preserve failures in the denominator and make the conditional summaries
    # explicit. Never silently drop difficult sparse-network repetitions.
    error_text <- if (!is.null(fit$error)) fit$error else "Too few valid bootstrap draws"
    return(cbind(base_row, data.frame(
      beta1_hat = NA_real_, beta2_hat = NA_real_, beta1_cover = NA,
      beta2_cover = NA, beta1_width = NA_real_, beta2_width = NA_real_,
      reject_network = NA, p_value = NA_real_,
      valid_bootstrap_draws = valid_draws, inference_ok = FALSE,
      inference_error = error_text
    )))
  }
  complete <- colSums(is.finite(fit$beta_boot)) == length(beta)
  beta_boot <- fit$beta_boot[, complete, drop = FALSE]

  ci1 <- quantile(beta_boot[1L, ], c(0.025, 0.975), na.rm = TRUE)
  ci2 <- quantile(beta_boot[2L, ], c(0.025, 0.975), na.rm = TRUE)
  test_boot <- m * colSums(
    (beta_boot[3:5, , drop = FALSE] - fit$beta_hat[3:5])^2
  )
  test_observed <- m * sum(fit$beta_hat[3:5]^2)
  # Follow the manuscript's upper-tail critical-value rejection rule.
  # Report a +1-adjusted upper-tail Monte Carlo p-value; with finite draws,
  # its threshold can differ slightly from the interpolated sample quantile.
  critical_value <- unname(quantile(test_boot, 0.95, na.rm = TRUE))
  p_value <- (1 + sum(test_boot >= test_observed)) / (1 + length(test_boot))

  cbind(base_row, data.frame(
    beta1_hat = fit$beta_hat[1L],
    beta2_hat = fit$beta_hat[2L],
    beta1_cover = beta[1L] >= ci1[1L] && beta[1L] <= ci1[2L],
    beta2_cover = beta[2L] >= ci2[1L] && beta[2L] <= ci2[2L],
    beta1_width = diff(ci1),
    beta2_width = diff(ci2),
    reject_network = test_observed > critical_value,
    p_value = p_value,
    valid_bootstrap_draws = ncol(beta_boot),
    inference_ok = TRUE, inference_error = ""
  ))
  }
  do.call(rbind, lapply(c(n, config$msizes), infer_sample))
}

summarize_results <- function(results) {
  groups <- split(results, interaction(results$delta, results$signal,
                                       results$m, drop = TRUE))
  safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  do.call(rbind, lapply(groups, function(current) {
    data.frame(
      signal = current$signal[1L], delta = current$delta[1L],
      n = current$n[1L], m = current$m[1L], method = current$method[1L],
      repetitions = nrow(current),
      successful_inferences = sum(current$inference_ok),
      failure_rate = mean(!current$inference_ok),
      # These estimates condition on successful inference; show failure_rate.
      coverage_beta1 = safe_mean(current$beta1_cover),
      coverage_beta2 = safe_mean(current$beta2_cover),
      average_width_beta1 = safe_mean(current$beta1_width),
      average_width_beta2 = safe_mean(current$beta2_width),
      rejection_rate = safe_mean(current$reject_network),
      valid_and_cover_beta1 = mean(current$beta1_cover %in% TRUE),
      valid_and_cover_beta2 = mean(current$beta2_cover %in% TRUE),
      average_degree = mean(current$mean_degree),
      average_isolated_fraction = mean(current$isolated_fraction),
      average_subspace_min_cosine = mean(current$subspace_min_cosine),
      average_valid_bootstrap_draws = mean(current$valid_bootstrap_draws),
      mean_beta1_hat = safe_mean(current$beta1_hat),
      mean_beta2_hat = safe_mean(current$beta2_hat)
    )
  }))
}

signals <- list(
  null = c(1, 2, 0, 0, 0),
  strong = c(1, 2, 1, 2, 1),
  weak = c(1, 2, 0, 0.01, 0)
)
signals <- signals[config$signals]

run_signal <- function(signal, run_config) {
  beta <- signals[[signal]]
  message("Running delta=", run_config$delta, ", ", signal,
          " network signal: ", run_config$nsim, " repetitions")
  checkpoint_path <- file.path(checkpoint_dir, paste0(
    "delta_", gsub("-", "m", format(run_config$delta, trim = TRUE)),
    "_", signal, ".rds"))
  # Permit changes in output location or worker count when resuming, but not
  # changes in the statistical configuration or implementation.
  statistical_config <- run_config[setdiff(names(run_config), c("output", "cores"))]
  pieces <- list()
  if (file.exists(checkpoint_path)) {
    prior <- readRDS(checkpoint_path)
    if (!identical(prior$statistical_config, statistical_config) ||
        !identical(prior$fingerprint, implementation_fingerprint))
      stop("Checkpoint configuration/source mismatch; use a new output folder")
    pieces <- prior$pieces
    message("Resuming ", signal, " after ", length(pieces), " repetitions")
  }
  completed <- length(pieces)
  if (completed < run_config$nsim) {
    for (first in seq.int(completed + 1L, run_config$nsim, by = 20L)) {
      ids <- seq.int(first, min(first + 19L, run_config$nsim))
      batch <- parallel::mclapply(
        ids, one_repetition, signal = signal, beta = beta, config = run_config,
        mc.cores = min(run_config$cores, length(ids)),
        mc.preschedule = TRUE, mc.set.seed = FALSE
      )
      failed <- vapply(batch, inherits, logical(1), "try-error")
      if (any(failed)) stop("Simulation worker failed: ", batch[[which(failed)[1L]]])
      pieces <- c(pieces, batch)
      save_checkpoint(list(statistical_config = statistical_config,
                           fingerprint = implementation_fingerprint,
                           pieces = pieces, updated_at = Sys.time()), checkpoint_path)
      message(format(Sys.time(), tz = "UTC"), " UTC: delta=", run_config$delta,
              " ", signal, " ", length(pieces), "/", run_config$nsim, " saved")
    }
  }
  do.call(rbind, pieces)
}

results <- do.call(rbind, lapply(config$delta, function(delta) {
  run_config <- config
  run_config$delta <- delta
  do.call(rbind, lapply(names(signals), run_signal, run_config = run_config))
}))
summary <- summarize_results(results)

tag <- paste0(
  config$mode,
  "_n", config$n,
  "_r", config$nsim,
  "_b", config$nboot,
  "_delta", paste(gsub("-", "m", format(config$delta, trim = TRUE)), collapse = "-"),
  "_ms", if (length(config$msizes)) paste(config$msizes, collapse = "-") else "none",
  "_", paste(config$signals, collapse = "-"), "_seed", config$seed,
  "_v3"
)
write.csv(results, file.path(config$output, paste0("strict_grdpg_repetitions_", tag, ".csv")), row.names = FALSE)
write.csv(summary, file.path(config$output, paste0("strict_grdpg_summary_", tag, ".csv")), row.names = FALSE)
saveRDS(
  list(config = config, implementation = "v3_full_graph_ASE_subset_multiplier",
       fingerprint = implementation_fingerprint,
       session_info = sessionInfo(),
       model = "strict continuous GRDPG(2,1)", results = results, summary = summary),
  file.path(config$output, paste0("strict_grdpg_", tag, ".rds"))
)

print(summary, row.names = FALSE)
message("Completed ", format(Sys.time(), tz = "UTC"), " UTC")
sink(type = "message")
sink()
close(log_connection)
