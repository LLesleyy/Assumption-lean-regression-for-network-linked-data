# Completed disjoint-subset rerun

Completed 2026-09-19. All 400 datasets in each of four sparsity regimes
finished. Both motifs were analyzed at m=25,50,100 and K=1,5,10, using
fixed disjoint node blocks. Previous experiments and manuscript files were
not modified.

## Findings

Across all 48 comparisons of K=5 or K=10 against K=1 at the same m,
sparsity and motif, empirical SD, mean interval width and RMSE decreased.
Coverage decreased in 46 comparisons and was unchanged in two (zero in
both cases). No comparison improved coverage.

For the local two-star coefficient at delta=-0.25 and m=25:

| K | Bias | Empirical SD | RMS estimated SE | Coverage | Mean width |
|---|---:|---:|---:|---:|---:|
| 1 | -0.304 | 0.766 | 0.661 | 0.8375 | 2.493 |
| 5 | -0.337 | 0.404 | 0.293 | 0.7100 | 1.141 |
| 10 | -0.330 | 0.323 | 0.207 | 0.5725 | 0.810 |

The triangle results have the same qualitative pattern; at m=25,
delta=-0.25, coverage is 0.8425, 0.6500 and 0.5300 for K=1,5,10.
These are new fixed-block results; the old randomly selected K=1 subsets
were different and their table entries should not be copied into this study.

## What the variance diagnostic establishes

At m=25, delta=-0.25, the oracle two-star empirical variance ratio
Var(K=10)/Var(K=1) is 0.09182, close to the independent reference 0.1.
The plug-in ratio is 0.17822. The empirical cross-block contribution to
the K=10 plug-in variance is 0.04222, with total variance 0.10458.
For the corresponding oracle, these are -0.000405 and 0.013668.
Finite Monte Carlo covariance estimates need not be exactly zero even
for independent oracle estimates; no sign or significance claim follows
merely from these estimates.

Thus disjointness does not produce exact finite-sample independence of
the graph-based estimates. The common graph and estimated sparsity still
matter. The independent-block sandwich variance has an asymptotic
justification under negligible network-feature error, not an exact
finite-sample justification for the actual estimator.

Estimated SEs also understate the empirical SD. Oracle HC0 intervals at
m=25 have coverage 0.8675 at K=1 and about 0.91 at K=10, so finite-block
sandwich calibration contributes even without graph estimation error.
At m=100, oracle coverage is about 0.925--0.945. Do not attribute the
actual coverage loss solely to structural bias, or describe these
intervals as a conservative finite-sample guarantee.

The design is useful as a fixed-K variance-stabilization study, not as
a demonstration of improved coverage or exact independence. The short
oracle-equivalence derivation is in `fixed_K_justification.tex`.
No results-based alteration of the requested grid or variance estimator
was made.

## Validation

- 230400 aggregate coefficient records; 576 cells, each with 400 valid
  replications; no fit failures.
- Every fixed block allocation checked for disjointness and exact node count.
- Every estimate and variance independently reconstructed from saved
  block-level fits; every numerical summary independently recomputed.
- All source-code hashes match the recorded configuration.
- Every graph matches the original saved density and mean degree.
- Every single-block fit agrees with the original moment fitter, maximum
  estimate/SE discrepancy 4.86e-13.
- Independent toy lm/sandwich calculations and direct Gaussian simulations
  validate the linearized covariance, including off-diagonal entries.
- Both PNG figures visually inspected: axes, labels, legends and uncertainty
  bars are visible, with shared scales across sparsity panels within each row.
  Coverage bars use Wilson Monte Carlo intervals, including at zero coverage.

## Deliverables

All combined numerical results and replacement figures are in
`results/formal_n6000_r400_combined/`:

- `disjoint_aggregation_tables.tex`: both plain journal-style tables.
- `disjoint_two_star_table.tex`, `disjoint_triangle_table.tex`: individual tables.
- `disjoint_aggregation_two_star.png`, `disjoint_aggregation_triangle.png`.
- `summary_network.csv`, `summary_all_coefficients.csv`: practical and oracle results.
- `variance_diagnostics.csv`: empirical variance decomposition and K=1 ratios.
- `paired_changes_vs_K1.csv`: paired changes with Monte Carlo SEs.
- `validation_checks.csv`, `replicate_network.rds`, `full_benchmark_summary.csv`.

Per-sparsity directories preserve seeds, configuration, fixed block lists,
individual block estimates and variances, and all coefficient-level outcomes.
