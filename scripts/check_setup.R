args <- commandArgs(trailingOnly = TRUE)
run_tests <- "--tests" %in% args

required <- c(
  "AID", "MASS", "Matrix", "RSpectra", "Rcpp", "RcppArmadillo",
  "abind", "combinat", "doParallel", "dplyr", "foreach", "ggplot2",
  "gridExtra", "igraph", "mvtnorm", "pracma", "prettyR", "tidyr"
)

installed <- rownames(installed.packages())
missing <- setdiff(required, installed)
if (length(missing)) {
  message("Missing optional/full-repository packages: ", paste(missing, collapse = ", "))
} else {
  message("All declared R packages are installed.")
}

files <- list.files("code", pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
errors <- lapply(files, function(path) {
  tryCatch({ parse(path); NULL }, error = function(e) conditionMessage(e))
})
bad <- !vapply(errors, is.null, logical(1))
if (any(bad)) {
  for (i in which(bad)) message(files[i], ": ", errors[[i]])
  stop("One or more R files could not be parsed.")
}
message("Parsed ", length(files), " R files successfully.")

if (run_tests) {
  if (!all(c("Rcpp", "RcppArmadillo") %in% installed)) {
    stop("Rcpp and RcppArmadillo are required for --tests.")
  }
  tests <- c(
    "code/revision-experiments/disjoint-aggregation/test_disjoint.R",
    "code/revision-experiments/disjoint-aggregation/test_graph_counts.R"
  )
  for (test in tests) {
    message("Running ", test)
    status <- system2(file.path(R.home("bin"), "Rscript"), test)
    if (!identical(status, 0L)) stop("Test failed: ", test)
  }
  message("All revision unit tests passed.")
}

