# Assumption-lean inference for network-linked data

This repository contains the simulation and real-data analysis code for
**“Assumption-lean inference for network-linked data”** by Wei Li, Nilanjan
Chakraborty, and Robert Lunde (manuscript BEJ2511-037).

It includes both the code used for the originally submitted experiments and
the additional experiments conducted for the revision. Compact versions of
the reported tables and figures are included under `results/reported`.

## Repository structure

```text
code/
  original-experiments/       Original manuscript simulation and analysis code
  revision-experiments/       Experiments added during revision
data/                          Data-access and preparation instructions only
results/reported/              Reported summary tables and figures
scripts/                       Reproducibility checks
REPRODUCIBILITY.md             Commands for reproducing each analysis
```

The revision experiments are:

1. comparison of full OLS, bias-corrected OLS, and downsampled OLS using local
   two-star and triangle covariates;
2. aggregation of estimators from fixed disjoint regression subsets;
3. bootstrap inference under a continuous, indefinite GRDPG; and
4. sensitivity of the school-conflict analysis to embedding dimensions
   `d = 2, 3, 4`.

## Software

The reported runs used R 4.2.1. A C++ compiler is required for experiments
that use `Rcpp`. Install the R dependencies with:

```r
install.packages(c(
  "AID", "MASS", "Matrix", "RSpectra", "Rcpp", "RcppArmadillo",
  "abind", "combinat", "doParallel", "dplyr", "foreach", "ggplot2",
  "gridExtra", "igraph", "mvtnorm", "pracma", "prettyR", "tidyr"
))
```

The revision experiments use only a subset of these packages. Run the local
preflight check from the repository root:

```sh
Rscript scripts/check_setup.R
```

To compile and exercise the two principal revision test suites as well:

```sh
Rscript scripts/check_setup.R --tests
```

## Reproduction

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the exact commands and the
recommended order of execution. The formal Monte Carlo studies are
computationally intensive; small diagnostic commands are provided first.

## School-conflict data

The source data are not redistributed in this repository. They are available
as ICPSR study 37070:

> Paluck, Elizabeth Levy, Hana R. Shepherd, and Peter Aronow. *Changing
> Climates of Conflict: A Social Network Experiment in 56 Schools, New Jersey,
> 2012–2013*. ICPSR 37070, version 2.
> https://doi.org/10.3886/ICPSR37070.v2

See [data/README.md](data/README.md) for placement and preprocessing details.
No individual-level school data are committed here.

## Results and storage policy

`results/reported` contains the small tables, summaries, and figures used in
the revision. Large repetition-level Monte Carlo objects, checkpoints, logs,
and processed individual-level data are intentionally excluded. The scripts
record seeds and configuration information and regenerate those files.

## Citation

Please cite the manuscript and this repository. Citation metadata are provided
in `CITATION.cff`.

