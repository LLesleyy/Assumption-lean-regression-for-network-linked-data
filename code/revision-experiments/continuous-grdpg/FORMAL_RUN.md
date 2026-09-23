# Formal experiment launched 2026-09-11

The denser run completed at 2:28 p.m. America/Chicago on 2026-09-11. Verified
tables and meeting interpretation are in dense_results/ beside this document.
The sparse run continues independently; consult its run.log for current progress.

Requested comparison: one denser regime and one sparse regime, with downsampling
only in the sparse regime. Inferential improvement is to be assessed from the
results, not assumed. No exact reconstruction of the old SBM table is claimed.

| Setting | Denser | Sparse |
|---|---|---|
| n | 3000 | 3000 |
| rho_n | n^(-0.3) | n^(-0.55) |
| Expected degree | 129.28 | 17.47 |
| Regression sample sizes | 3000 | 3000, 50, 100, 200, 400 |
| Independent repetitions per signal | 400 | 400 |
| Bootstrap draws per inference | 2000 | 2000 |
| Signals | null, weak, strong | null, weak, strong |
| Workers | 2 | 2 |

The two cases run concurrently with two workers each. There are 1200 independent
graph repetitions within each case, 2400 total case-specific graph realizations.
Within each sparse repetition all m values share the full graph, response and
ASE. The same base seeds across sparsity settings provide shared random-number
streams; different repetitions within a setting are independent. The 400
repetitions are not 400 independent bootstrap draws: each repetition has its own
graph, covariates and response, followed by 2000 bootstrap draws per method.

Model, response, known-rho scaling, CI construction, upper-tail test and
diagnostic definitions are in README.md. The estimator
and bootstrap are identical between full and subset fits, with regression
sample size m used for the subset moments and test statistic.

## Frozen commands

Run from the project root:

```sh
Rscript code/revision-experiments/continuous-grdpg/strict_grdpg_bootstrap.R --mode=full --n=3000 --nsim=400 --nboot=2000 --delta=-0.3 --msizes=none --signals=null,weak,strong --seed=2601 --cores=2 --output=code/revision-experiments/continuous-grdpg/results/formal_n3000_b2000_r400/dense
Rscript code/revision-experiments/continuous-grdpg/strict_grdpg_bootstrap.R --mode=full --n=3000 --nsim=400 --nboot=2000 --delta=-0.55 --msizes=50,100,200,400 --signals=null,weak,strong --seed=2601 --cores=2 --output=code/revision-experiments/continuous-grdpg/results/formal_n3000_b2000_r400/sparse
```

Each output folder contains run.log and checkpoints saved after every 20
repetitions. The final repetition-level CSV, summary CSV, and RDS appear when
all three signal cases finish. Repeating the identical command resumes saved
batches; a changed statistical configuration or source fingerprint is rejected.
Do not edit the running experiment script before completion/resumption. Logs
show the process ID and checkpoint progress. No inference results are selected
or suppressed because coverage, size or power is unfavorable.

## Planned reporting

For every signal and method report conventional beta_x coverage and average CI
width. Report rejection as type-I error under the exact null and as power under
weak/strong alternatives. Show inference failure rates, valid draw counts, mean
degree and spectral recovery diagnostics. Keep conditional-on-success inference
summaries labeled if any failures occur. Retain all attempted replications.
Nominal coverage is not required to demonstrate a relative practical benefit,
but improved coverage must be assessed together with width, size and power.
Slow root-m convergence explains lost precision; remaining miscalibration may
also reflect residual network error or the finite-sample bootstrap approximation.
