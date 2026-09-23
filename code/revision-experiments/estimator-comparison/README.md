# Coverage-width tradeoff for two-star and triangle covariates

This experiment uses linear responses with true coefficients `(1,20,3,2)` and
the graphon `w(u,v)=u+v`. It treats "triangle" as normalized local triangle
frequency, not the transitivity ratio.

## Main design

- `n=6000`, 400 independent graphs per sparsity, and 999 Gaussian multiplier
  draws for the full-graph structural bootstrap.
- `delta=-.55,-.40,-.30,-.25`, with `rho=n^delta`.
- Downsample sizes `m=25,40,65,100,160,250,400,650,1000`; all network
  covariates are constructed from the full graph before regression rows are
  selected. Subsets are nested within a graph.
- Primary comparison: full uncorrected OLS, full bias-corrected OLS, and
  downsampled uncorrected OLS. Corrected moments within each downsample are
  retained as a clearly labeled diagnostic.
- Every cell records coverage, interval width, estimate, bias, RMSE, fit status,
  conditioning, mean degree, and reproducible seeds for all four coefficients.

The downsampled confidence interval uses the exact conditional quantiles of the
linearized iid Gaussian multiplier distribution. Since the multiplier is
Gaussian, these quantiles are analytic and contain no avoidable inner Monte
Carlo error. Full estimators use the order-aware first-projection multiplier
bootstrap, including re-estimation of sparsity. The leading corrected squared
moments are the nonoverlapping four-star term for local two-stars and the
nonoverlapping bowtie term for local triangles.

`count_subgraphs.cpp` generates each graph without materializing an R adjacency
matrix and uses bitsets to count triangles. Exact finite-sample identities check
all-role first projections and both order-five corrected moments on the first
replicate of every sparsity. Failed or ill-conditioned fits are recorded rather
than regularized or removed.

The four sparsities are run separately so that checkpoint files are available
every ten graph replications. After completion, run:

```sh
Rscript combine_and_plot.R results
```

The combined output includes raw results, all-coefficient summaries,
network-coefficient summaries, graph diagnostics, and separate coverage-width
figures for the two motifs.

## Completed main run

All 1,600 graphs completed, producing 256,000 coefficient-level records. Every
populated cell contains 400 replications and all reported intervals are finite.
The independent audit in `validate_results.R` exactly reproduced coverage,
average width, mean estimate, and RMSE from the raw records. Both final figures
were visually inspected.

For ordinary downsampling, the highest network-coefficient coverage over the
predeclared grid was:

| delta | Two-star: coverage (m; width) | Triangle: coverage (m; width) |
|---:|---:|---:|
| -.55 | .300 (25; 6.55) | .073 (25; 7.19) |
| -.40 | .708 (25; 4.17) | .673 (25; 4.42) |
| -.30 | .808 (25; 3.01) | .778 (40; 2.54) |
| -.25 | .858 (25; 2.54) | .855 (40; 2.09) |

Full uncorrected OLS coverage for two-star/triangle was respectively 0/0 at
delta=-.55, 0/0 at -.40, .295/.313 at -.30, and .670/.723 at -.25. Thus small
subsets materially recover coverage, but none reaches .95 in this design.
Coverage then falls as m grows because intervals narrow while attenuation
persists.

Full corrected OLS coverage was .990/1.000, .975/.965, .963/.953, and
.965/.940 for two-star/triangle across increasing density. It was also much
more efficient than the small-m intervals. At delta=-.55, however, its
structural bootstrap standard errors were conservative relative to empirical
Monte Carlo standard deviations (mean-SE/empirical-SD ratios 1.36 for two-star
and 1.68 for triangle); the overcoverage should not be presented as ideal
calibration.

The corrected-downsample diagnostic is deliberately omitted from the primary
figures. It is nearly centered but can be extremely unstable at small m; for
triangle at delta=-.55 and m=25, its mean interval width is 769 because a few
near-singular corrected moment matrices create extreme intervals. These values
remain in the raw data and summaries.
