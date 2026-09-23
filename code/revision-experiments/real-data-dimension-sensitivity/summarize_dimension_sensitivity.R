args <- commandArgs(trailingOnly = TRUE)
results_root <- if (length(args) > 0) args[[1]] else "dimension_sensitivity_results"

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

dims <- 2:4
key_coefficients <- c(
  "(Intercept)",
  "GRC_NUMERIC",
  "RETURN(1) Returning to school",
  "TREAT_BINARY",
  "LIVEWB(1) Yes: lives with both parents",
  "Race_overallAsian",
  "Race_overallBlack",
  "Race_overallHispanic",
  "Race_overallWhite",
  "GENDER(1) Boy",
  "SCHID_collapsed35",
  "SCHID_collapsedOther"
)

display_names <- c(
  "(Intercept)" = "Intercept",
  "GRC_NUMERIC" = "Grade",
  "RETURN(1) Returning to school" = "Returning to school",
  "TREAT_BINARY" = "Treatment",
  "LIVEWB(1) Yes: lives with both parents" = "Lives with both parents",
  "Race_overallAsian" = "Race: Asian",
  "Race_overallBlack" = "Race: Black",
  "Race_overallHispanic" = "Race: Hispanic",
  "Race_overallWhite" = "Race: White",
  "GENDER(1) Boy" = "Gender: Boy",
  "SCHID_collapsed35" = "School 35",
  "SCHID_collapsedOther" = "Other schools"
)

coefficient_comparison <- bind_rows(lapply(dims, function(d) {
  path <- file.path(results_root, paste0("d", d), "coefficient_pvalue_table.csv")
  x <- read.csv(path, check.names = FALSE)
  x$Dimension <- d
  x
})) %>%
  filter(Coefficient %in% key_coefficients) %>%
  mutate(
    Variable = unname(display_names[Coefficient]),
    Variable = factor(Variable, levels = unname(display_names[key_coefficients]))
  ) %>%
  arrange(Variable, Dimension)

write.csv(
  coefficient_comparison,
  file.path(results_root, "key_coefficient_sensitivity.csv"),
  row.names = FALSE
)

network_tests <- bind_rows(lapply(dims, function(d) {
  read.csv(file.path(results_root, paste0("d", d), "overall_network_test.csv"))
}))
write.csv(network_tests, file.path(results_root, "overall_network_test_comparison.csv"), row.names = FALSE)

school_effects <- bind_rows(lapply(dims, function(d) {
  x <- read.csv(file.path(results_root, paste0("d", d), "school_network_effect_summary.csv"))
  x$Dimension <- d
  x
}))
write.csv(school_effects, file.path(results_root, "school_network_effect_comparison.csv"), row.names = FALSE)

node_effects <- lapply(dims, function(d) {
  x <- readRDS(file.path(results_root, paste0("d", d), "analysis_results.rds"))$network_effect
  data.frame(
    SCHID = x$SCHID,
    StudentIndex = x$StudentIndex,
    NetworkEffect = x$NetworkEffect
  )
})

stopifnot(all(vapply(node_effects[-1], function(x) {
  identical(x$SCHID, node_effects[[1]]$SCHID) &&
    identical(x$StudentIndex, node_effects[[1]]$StudentIndex)
}, logical(1))))

network_effect_matrix <- do.call(cbind, lapply(node_effects, `[[`, "NetworkEffect"))
colnames(network_effect_matrix) <- paste0("d", dims)
network_effect_correlations <- cor(network_effect_matrix)
write.csv(
  network_effect_correlations,
  file.path(results_root, "student_network_effect_correlations.csv")
)

sensitivity_plot <- ggplot(
  coefficient_comparison %>% filter(!Variable %in% c("Intercept", "School 35", "Other schools")),
  aes(x = factor(Dimension), y = Estimate, group = Variable)
) +
  geom_hline(yintercept = 0, color = "grey60", size = 0.35) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.12, color = "#24527A") +
  geom_line(color = "#24527A", size = 0.7) +
  geom_point(color = "#24527A", size = 1.8) +
  facet_wrap(~Variable, scales = "free_y", ncol = 3) +
  labs(
    x = "Embedding dimension per school",
    y = "Estimate"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")
ggsave(
  file.path(results_root, "key_coefficient_sensitivity.png"),
  sensitivity_plot,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

rms_plot <- ggplot(school_effects, aes(x = factor(SCHID), y = RMSNetworkEffect,
                                       color = factor(Dimension), group = Dimension)) +
  geom_line(size = 0.65) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("2" = "#24527A", "3" = "#B36B00", "4" = "#8A4F7D")) +
  labs(
    title = "Magnitude of fitted network effects by school",
    subtitle = "Root-mean-square student-level network effect",
    x = "School ID",
    y = "RMS network effect",
    color = "Dimension"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1))
ggsave(
  file.path(results_root, "school_network_effect_magnitude_sensitivity.png"),
  rms_plot,
  width = 10,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

for (d in dims) {
  x <- node_effects[[which(dims == d)]]
  x$SchoolID <- factor(x$SCHID, levels = unique(x$SCHID))
  network_plot <- ggplot(x, aes(StudentIndex, NetworkEffect, color = SchoolID)) +
    geom_point(alpha = 0.75, size = 1.1) +
    geom_hline(yintercept = 0, color = "grey45", size = 0.4) +
    labs(
      title = paste0("Student-level network effects, d = ", d),
      x = "Student index (ordered by school)",
      y = expression(hat(Z) %*% hat(beta)[z]),
      color = "School ID"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")
  ggsave(
    file.path(results_root, paste0("d", d), "student_network_effects_by_school.png"),
    network_plot,
    width = 10,
    height = 6,
    dpi = 300,
    bg = "white"
  )
}

eigen_spectrum <- read.csv(file.path(results_root, "school_eigen_spectrum_top10.csv"))
aggregate_spectrum <- eigen_spectrum %>%
  group_by(Component) %>%
  summarise(
    MedianRelativeMagnitude = median(RelativeToFirst),
    Q25 = quantile(RelativeToFirst, 0.25),
    Q75 = quantile(RelativeToFirst, 0.75),
    .groups = "drop"
  )

eigen_gaps <- bind_rows(lapply(split(eigen_spectrum, eigen_spectrum$School), function(x) {
  x <- x[order(x$Component), ]
  data.frame(
    School = x$School[1],
    GapAfter = 1:5,
    EigenvalueRatio = x$AbsEigenvalue[1:5] / x$AbsEigenvalue[2:6]
  )
}))
eigen_gap_summary <- eigen_gaps %>%
  group_by(GapAfter) %>%
  summarise(
    MedianRatio = median(EigenvalueRatio),
    Q25 = quantile(EigenvalueRatio, 0.25),
    Q75 = quantile(EigenvalueRatio, 0.75),
    .groups = "drop"
  )
largest_gap_by_school <- eigen_gaps %>%
  group_by(School) %>%
  slice_max(EigenvalueRatio, n = 1, with_ties = FALSE) %>%
  ungroup()
write.csv(eigen_gap_summary, file.path(results_root, "eigen_gap_summary.csv"), row.names = FALSE)
write.csv(largest_gap_by_school, file.path(results_root, "largest_eigen_gap_by_school.csv"), row.names = FALSE)

aggregate_plot <- ggplot(aggregate_spectrum, aes(Component, MedianRelativeMagnitude)) +
  geom_ribbon(aes(ymin = Q25, ymax = Q75), fill = "#9BB7D4", alpha = 0.35) +
  geom_line(color = "#24527A", size = 0.9) +
  geom_point(color = "#24527A", size = 2) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "#B36B00") +
  scale_x_continuous(breaks = seq_len(max(aggregate_spectrum$Component))) +
  labs(
    x = "Eigencomponent ordered by absolute eigenvalue",
    y = "Relative absolute eigenvalue"
  ) +
  theme_minimal(base_size = 11)
ggsave(file.path(results_root, "aggregate_eigen_spectrum.png"), aggregate_plot,
       width = 8, height = 5, dpi = 300, bg = "white")

school_spectrum_plot <- ggplot(eigen_spectrum, aes(Component, RelativeToFirst, group = School)) +
  geom_line(color = "#24527A", alpha = 0.35, size = 0.45) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "#B36B00") +
  scale_x_continuous(breaks = seq_len(max(eigen_spectrum$Component))) +
  labs(
    title = "Normalized adjacency spectrum for each school",
    subtitle = "Each line is one school; dashed line marks the current d = 2 cutoff",
    x = "Eigencomponent ordered by absolute eigenvalue",
    y = "Relative absolute eigenvalue"
  ) +
  theme_minimal(base_size = 11)
ggsave(file.path(results_root, "school_eigen_spectra_overlay.png"), school_spectrum_plot,
       width = 8, height = 5, dpi = 300, bg = "white")

format_num <- function(x, digits = 4) formatC(x, format = "f", digits = digits)
format_p <- function(x) {
  ifelse(x < 0.0001, "$<0.0001$", format_num(x, 4))
}

table_rows <- lapply(unname(display_names[key_coefficients]), function(v) {
  z <- coefficient_comparison %>% filter(as.character(Variable) == v) %>% arrange(Dimension)
  paste0(
    v, " & ",
    paste(sprintf("%s (%s), %s", format_num(z$Estimate), format_num(z$Bootstrap_SE), format_p(z$P_Value)),
          collapse = " & "),
    " \\\\"
  )
})

latex_lines <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\small",
  "\\caption{Sensitivity of conventional and selected school-effect estimates to the embedding dimension. Each entry reports estimate (bootstrap SE), p-value.}",
  "\\label{tab:real-data-dimension-sensitivity}",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  "Variable & $d=2$ & $d=3$ & $d=4$ \\\\",
  "\\midrule",
  unlist(table_rows),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(latex_lines, file.path(results_root, "dimension_sensitivity_table.tex"))

cat("Sensitivity summaries written to", normalizePath(results_root), "\n")
