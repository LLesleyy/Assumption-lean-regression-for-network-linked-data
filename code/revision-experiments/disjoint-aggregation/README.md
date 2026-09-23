# Fixed-disjoint-block aggregation study

This separate study replaces neither the old experiment nor any manuscript file.
It implements the requested grid m=25,50,100 and K=1,5,10, for both local
two-star and triangle covariates, at n=6000 and delta=-0.55,-0.40,-0.30,-0.25.
There are 400 Monte Carlo datasets per sparsity, replayed from the original
comparison's recorded data seeds. Common seeds pair sparsity regimes.

## Allocation and estimation

I_k = {(k-1)m+1,...,km}. These deterministic label blocks are fixed independently
of the data, disjoint within each m, and identical in every replication. The
first K blocks are used, so comparisons across K are paired. At most 1000 nodes
enter regression; all 6000 nodes still contribute to network-feature estimation.
Original node labels are generated iid, so selecting these fixed labels does
not select on the outcome, features, degree, or latent position.

The point estimator is the equal-weight mean of K separately fitted OLS
estimators, not OLS fitted once to the union and not the mean of CI endpoints.
The graphon, response, and covariate models are unchanged from
`../estimator-comparison/coverage_experiment.R`:

- w(u,v)=u+v, rho=n^delta, no probability truncation active on this grid.
- Clean center-rooted two-star feature: (xi+1/2)^2.
- Clean triangle feature: xi^2+7xi/6+1/3.
- Y=1+20 Z+3 X1+2 X2+epsilon, epsilon standard normal.
- (X1,X2) Gaussian with means (1,3), covariance [[1,.6],[.6,4]].
- Observed features use full-graph counts and the empirical edge density.

`disjoint_plugin` is the implementable method. `oracle_diagnostic` replaces
estimated features by the true node-wise features, on identical datasets and
blocks. It is a diagnostic, not a competing practical estimator.

## Variance: what is exact and what is asymptotic

For each block, use its own OLS fit, residuals and HC0 covariance V_k, matching
the original single-subset inference. Set V_aggregate = sum_k V_k / K^2.
This is exactly the conditional covariance of the linearized Gaussian
multiplier perturbation using independent multipliers on the disjoint nodes.
Normal quantiles are analytic; there is no inner bootstrap simulation count.
Features and empirical sparsity are held fixed in these multiplier perturbations.

Under the paper's iid-node model (conditional on the common mixing variable
when applicable), the oracle block estimators are independent and identically
distributed. If their second moments exist, the exact identity is

    Var(mean_k beta_oracle,k) = Var(beta_oracle,1) / K.

The graph-based block estimators are NOT exactly independent: they use the same
graph and sparsity estimate even when the node sets are disjoint. A valid
asymptotic justification uses the downsampling theorem's oracle-equivalence
argument. For each of finitely many blocks,

    beta_hat,k - beta* = m^(-1) sum_{i in I_k} psi_i + o_P(m^(-1/2)).

For fixed K, averaging these expansions gives

    sqrt(Km) (beta_bar - beta*)
      = (Km)^(-1/2) sum_{i in union I_k} psi_i + o_P(1).

Thus its first-order covariance is Sigma/(Km), and V_aggregate consistently
estimates it provided the block sandwich consistency assumptions also hold.
This is a finite-K extension using the oracle approximation, not a claim that
marginal CLTs alone imply independence. Fixed K changes the variance constant,
not the convergence-rate exponent. Letting K grow requires additional uniform
remainder control. For a feature bound a_n, sqrt(Km)*a_n -> 0 is a useful
sufficient aggregate-scale requirement; with fixed K it is equivalent to
sqrt(m)*a_n -> 0. The finite grid itself does not certify either asymptotic
negligibility or accurate finite-sample coverage.

Disjointness removes direct repeated-node terms, not all graph uncertainty.
Nor does averaging structurally correct bias: identically distributed block
estimates have the same expectation as their average when expectations exist.
Keep empirical SD, RMS estimated SE, bias and coverage distinct in interpretation.
No data-driven corrections or tuning are performed after viewing outcomes.

## Verification and outputs

- `test_disjoint.R`: independent lm/sandwich calculation, agreement with the
  prior shared-node implementation on disjoint blocks, simulated Gaussian
  covariance, rejection of overlapping blocks, and singular-fit propagation.
- The unchanged graph-count helper is tested against explicit toy adjacency
  matrices by `test_graph_counts.R`.
- Every dataset's edge density and mean degree match the original study.
- Every K=1 fit is checked against the original moment-fitting function on
  the SAME fixed block. Old K=1 table entries need not be identical because
  the old experiment used a random subset instead of these fixed blocks.
- `prepare_results.R` independently reconstructs ALL aggregate estimates and
  variances from block records and recomputes every summary from replicate data.
- Configuration files retain exact code hashes and seeds. Block-level and
  aggregate-level RDS files preserve the audit trail. Nothing is published.

Results under `results/formal_n6000_r400_combined` contain both motif tables,
figures, all coefficient summaries, oracle checks, empirical cross-block
covariance diagnostics, paired changes and full-sample reference summaries.
An empirical covariance decomposition is an accounting check, not evidence
that the finite-sample graph covariances must be zero.

## Reproduce (use fresh output directories)

```sh
Rscript code/revision-experiments/disjoint-aggregation/test_disjoint.R
Rscript code/revision-experiments/disjoint-aggregation/test_graph_counts.R
Rscript code/revision-experiments/disjoint-aggregation/run_disjoint.R code/revision-experiments/disjoint-aggregation/results/formal_n6000_r400_delta_m055 400 -0.55
Rscript code/revision-experiments/disjoint-aggregation/run_disjoint.R code/revision-experiments/disjoint-aggregation/results/formal_n6000_r400_delta_m040 400 -0.40
Rscript code/revision-experiments/disjoint-aggregation/run_disjoint.R code/revision-experiments/disjoint-aggregation/results/formal_n6000_r400_delta_m030 400 -0.30
Rscript code/revision-experiments/disjoint-aggregation/run_disjoint.R code/revision-experiments/disjoint-aggregation/results/formal_n6000_r400_delta_m025 400 -0.25
Rscript code/revision-experiments/disjoint-aggregation/prepare_results.R
```

The result-preparation script reads the four formal directory names used in
this run. R, Rcpp and a working C++ compiler suffice; plots use base R.
