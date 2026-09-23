#!/usr/bin/env Rscript
# Read-only validation of the completed experiment, followed by table exports.
args <- commandArgs(trailingOnly = TRUE)
source_dir <- if (length(args)) args[1L] else
  "code/revision-experiments/continuous-grdpg/results/formal_n3000_b2000_r400/dense"
output_dir <- if (length(args) >= 2L) args[2L] else
  file.path(dirname(source_dir), "dense-summary")
files <- list.files(source_dir, pattern = "^strict_grdpg_.*[.]rds$", full.names = TRUE)
if (length(files) != 1L) stop("Exactly one completed dense RDS is required")
run <- readRDS(files)
d <- run$results
stopifnot(run$config$n == 3000L, run$config$nboot == 2000L,
          run$config$nsim == 400L, identical(run$config$delta, -0.3),
          length(run$config$msizes) == 0L, nrow(d) == 1200L,
          !anyDuplicated(d[c("signal", "repetition")]),
          all(d$method == "full"), all(d$n == 3000L), all(d$m == 3000L),
          all(d$inference_ok), all(d$valid_bootstrap_draws == 2000L),
          all(is.finite(d$beta1_width)), all(is.finite(d$beta2_width)),
          all(d$beta1_width > 0), all(d$beta2_width > 0))
signals <- c("null", "weak", "strong")
labels <- c(null = "Null", weak = "Weak", strong = "Strong")
metrics <- do.call(rbind, lapply(signals, function(signal) {
  x <- d[d$signal == signal, ]
  stopifnot(nrow(x) == 400L, identical(sort(x$repetition), seq_len(400L)))
  data.frame(signal = signal, repetitions = nrow(x),
    covered_beta1 = sum(x$beta1_cover), coverage_beta1 = sum(x$beta1_cover) / nrow(x),
    average_width_beta1 = sum(x$beta1_width) / nrow(x),
    covered_beta2 = sum(x$beta2_cover), coverage_beta2 = sum(x$beta2_cover) / nrow(x),
    average_width_beta2 = sum(x$beta2_width) / nrow(x),
    rejected = sum(x$reject_network), rejection_rate = sum(x$reject_network) / nrow(x),
    mean_degree = mean(x$mean_degree), mean_beta1_hat = mean(x$beta1_hat),
    mean_beta2_hat = mean(x$beta2_hat))
}))
# Independently recomputed aggregates must agree with the runtime summary.
for (signal in signals) {
  old <- run$summary[run$summary$signal == signal, ]
  fresh <- metrics[metrics$signal == signal, ]
  stopifnot(nrow(old) == 1L)
  for (key in c("coverage_beta1", "coverage_beta2", "average_width_beta1",
                "average_width_beta2", "rejection_rate"))
    stopifnot(abs(old[[key]] - fresh[[key]]) < 1e-12)
}
uncertainty <- do.call(rbind, lapply(signals, function(signal) {
  x <- d[d$signal == signal, ]
  do.call(rbind, lapply(c("beta1_cover", "beta2_cover", "reject_network"), function(key) {
    successes <- sum(x[[key]])
    p <- successes / nrow(x)
    ci <- binom.test(successes, nrow(x))$conf.int
    data.frame(signal = signal, metric = key, successes = successes,
               repetitions = nrow(x), estimate = p,
               monte_carlo_se = sqrt(p * (1 - p) / nrow(x)),
               exact_95_lower = ci[1L], exact_95_upper = ci[2L])
  }))
}))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(metrics, file.path(output_dir, "dense_results.csv"), row.names = FALSE)
write.csv(uncertainty, file.path(output_dir, "dense_monte_carlo_uncertainty.csv"), row.names = FALSE)
display_metrics <- metrics[match(c("strong", "weak", "null"), metrics$signal), ]
setting_labels <- c(strong = "Strong network effect", weak = "Weak network effect",
                    null = "No network effect")
table_rows <- vapply(seq_len(nrow(display_metrics)), function(i) {
  r <- display_metrics[i, ]
  sprintf("| %s | %.4f | %.4f | %.4f |", setting_labels[r$signal],
          r$coverage_beta1, r$coverage_beta2, r$rejection_rate)
}, character(1))
width_rows <- vapply(seq_len(nrow(display_metrics)), function(i) {
  r <- display_metrics[i, ]
  sprintf("| %s | %.4f | %.4f |", setting_labels[r$signal],
          r$average_width_beta1, r$average_width_beta2)
}, character(1))
null_rate <- metrics$rejection_rate[metrics$signal == "null"]
weak_rate <- metrics$rejection_rate[metrics$signal == "weak"]
strong_rate <- metrics$rejection_rate[metrics$signal == "strong"]
note <- c(
  "# Denser continuous GRDPG: completed formal results",
  "",
  "All 1,200 graph repetitions are complete: 400 each under null, weak, and strong network signals.",
  "",
  "n=3000; rho_n=n^(-0.3); 2000 Gaussian multiplier draws per repetition; nominal CI coverage 95%; network-test level 5%.",
  "",
  "| Setting | Coverage beta1 | Coverage beta2 | Power (network effect) |",
  "|---|---:|---:|---:|", table_rows,
  "",
  "In the empirical-power column, the null row reports type-I error (size); weak and strong rows report power.",
  "",
  "## Companion table: mean CI widths",
  "",
  "| Setting | Mean width beta1 | Mean width beta2 |",
  "|---|---:|---:|", width_rows,
  "",
  "## Meeting interpretation",
  "",
  sprintf("- Coverage ranges from %.2f%% to %.2f%% across the six coefficient/signal combinations.",
          100 * min(metrics$coverage_beta1, metrics$coverage_beta2),
          100 * max(metrics$coverage_beta1, metrics$coverage_beta2)),
  sprintf("- The exact-null rejection rate is %.2f%% versus nominal 5%%. The weak-signal rejection rate is %.2f%%; it measures power, not null calibration.",
          100 * null_rate, 100 * weak_rate),
  sprintf("- Strong-signal power is %.2f%% in these simulations. An observed 100%% does not establish population power exactly equal to one.",
          100 * strong_rate),
  "- At a true rate of 95% or 5%, 400 repetitions give a Monte Carlo standard error of about 1.09 percentage points (roughly +/-2.14 points for a 95% normal reference band). Exact binomial intervals for the observed rates are in dense_monte_carlo_uncertainty.csv.",
  "",
  "## Model and limits",
  "",
  "The continuous latent coordinates are Z=(0.69+0.05 cos(2 pi xi), 0.45 sin(2 pi xi), 0.40 cos(4 pi xi)), xi~Uniform(0,1), with signature diag(1,1,-1). Kernel eigenvalues are (0.47735,0.10125,-0.08), so this example is outside SBM and ordinary RDPG.",
  "",
  "The conventional coefficients are (1,2). Conventional covariates have zero means, unit variances and correlation 0.3, and are independent of the latent positions. The network coefficients are (0,0,0), (0,0.01,0), and (1,2,1). The response model is linear with independent N(0,1) errors and no added intercept.",
  "",
  sprintf("The realized average degree is %.2f. ASE uses known rho_n, three eigenpairs ordered by absolute magnitude, and no empirical Procrustes alignment for inference. These results concern this particular model and the known-sparsity implementation. They do not establish performance across general sparse GRDPGs or validate the estimated-sparsity implementation.", mean(d$mean_degree)),
  "",
  "The sparse/downsampling experiment is separate and must be reported from its own results. These dense results do not imply that downsampling will improve sparse inference.",
  "",
  "## Verification and provenance",
  "",
  "- Confirmed 400 distinct repetitions per signal, n=3000, and 2000 valid bootstrap draws in every inference.",
  "- No inference failures or missing interval widths; no repetitions excluded.",
  "- Recomputed coverage, mean widths and rejection rates from individual repetitions and reconciled them to the saved summary.",
  paste0("- Source: `", normalizePath(files), "`."),
  paste0("- Source fingerprint: `", run$fingerprint, "`."),
  "- Prior SBM draft settings and results are not an exact baseline for numerical comparisons because the accepted scaling/testing corrections changed that implementation."
)
writeLines(note, file.path(output_dir, "dense_meeting_summary.md"))
tex_rows <- vapply(seq_len(nrow(display_metrics)), function(i) {
  r <- display_metrics[i, ]
  sprintf("%s & %.4f & %.4f & %.4f \\\\", setting_labels[r$signal],
          r$coverage_beta1, r$coverage_beta2, r$rejection_rate)
}, character(1))
tex <- c(
  "% Requires \\usepackage{booktabs}. Table only; suitable for \\input{}.",
  "\\begin{table}[htbp]", "\\centering", "\\small",
  "\\caption{\\textit{Bootstrap coverage and power under different network signal settings in a continuous GRDPG}}",
  "\\label{tab:strict-grdpg-dense}", "\\renewcommand{\\arraystretch}{1.2}",
  "\\begin{tabular}{lccc}", "\\toprule",
  "Setting & Coverage ($\\beta_1$) & Coverage ($\\beta_2$) & Power (network effect) \\\\",
  "\\midrule", tex_rows, "\\bottomrule", "\\end{tabular}",
  "\\par\\smallskip", "\\begin{minipage}{0.98\\linewidth}\\footnotesize",
  "Entries are empirical coverage probabilities of 95\\% bootstrap confidence intervals for the conventional covariate coefficients $(\\beta_1,\\beta_2)=(1,2)$ and rejection frequencies of the level-$0.05$ test for an overall network effect. Each setting uses 400 independent repetitions with $n=3000$, $\\rho_n=n^{-0.3}$, and 2000 multiplier bootstrap draws per repetition. The last column reports power under the alternatives and type-I error in the no-network-effect setting. The adjacency spectral embedding is normalized using known $\\rho_n$.",
  "\\end{minipage}", "\\end{table}"
)
writeLines(tex, file.path(output_dir, "dense_results_table.tex"))
width_tex_rows <- vapply(seq_len(nrow(display_metrics)), function(i) {
  r <- display_metrics[i, ]
  sprintf("%s & %.4f & %.4f \\\\", setting_labels[r$signal],
          r$average_width_beta1, r$average_width_beta2)
}, character(1))
width_tex <- c(
  "% Requires \\usepackage{booktabs}. Companion table; preserves the width results.",
  "\\begin{table}[htbp]", "\\centering", "\\small",
  "\\caption{\\textit{Average widths of 95\\% bootstrap confidence intervals in a continuous GRDPG}}",
  "\\label{tab:strict-grdpg-dense-width}", "\\renewcommand{\\arraystretch}{1.2}",
  "\\begin{tabular}{lcc}", "\\toprule",
  "Setting & Mean width ($\\beta_1$) & Mean width ($\\beta_2$) \\\\",
  "\\midrule", width_tex_rows, "\\bottomrule", "\\end{tabular}",
  "\\par\\smallskip", "\\begin{minipage}{0.98\\linewidth}\\footnotesize",
  "Entries are average widths of 95\\% bootstrap confidence intervals for the conventional covariate coefficients $(\\beta_1,\\beta_2)=(1,2)$ over 400 independent repetitions per setting. Here $n=3000$, $\\rho_n=n^{-0.3}$, and each repetition uses 2000 multiplier bootstrap draws. The adjacency spectral embedding is normalized using known $\\rho_n$. Coverage and network-effect testing results are reported in Table~\\ref{tab:strict-grdpg-dense}.",
  "\\end{minipage}", "\\end{table}"
)
writeLines(width_tex, file.path(output_dir, "dense_widths_table.tex"))
cat("Completed and verified dense summary saved to", output_dir, "\n")
print(metrics, row.names = FALSE)
