# Strict-GRDPG simulation

The formal comparison was launched on 2026-09-11. See `FORMAL_RUN.md` for the
frozen commands: n=3000, B=2000, 400 repetitions per signal, delta=-0.3 full
regression only, delta=-0.55 full plus m=50,100,200,400. Read the dense/sparse
run.log files under results/formal_n3000_b2000_r400 for current progress.

This experiment extends the Section 5.2 SBM study to a continuous and genuinely
indefinite generalized random dot product graph.

## Model

For independent `xi_i ~ Uniform(0, 1)`, define

```text
Z_i1 = 0.69 + 0.05 cos(2 pi xi_i)
Z_i2 = 0.45 sin(2 pi xi_i)
Z_i3 = 0.40 cos(4 pi xi_i)
```

and generate edges independently conditional on the latent positions with

```text
P_ij = rho_n (Z_i1 Z_j1 + Z_i2 Z_j2 - Z_i3 Z_j3).
```

The three coordinate functions are mutually orthogonal in `L2[0,1]`. The
corresponding signed eigenvalues are

```text
0.47735, 0.10125, -0.08000.
```

The negative eigenvalue makes this a strict GRDPG rather than an ordinary RDPG.
The continuous latent distribution makes it non-SBM. The kernel is bounded
between `0.0471` and `0.9101`, and its mean is `0.4761`, almost exactly matching
the submitted SBM's mean connectivity `0.47625`. Therefore the same
`rho_n = n^(-0.3)` yields a controlled comparison at essentially the same
expected degree.

## Outcomes and inference

The response follows the linear model stated in the manuscript, with the
accepted correction to clean-position scaling:

- two correlated standard-normal conventional covariates with correlation 0.3;
- exact-null coefficients `(1, 2, 0, 0, 0)` for type-I error;
- strong coefficients `(1, 2, 1, 2, 1)`;
- weak coefficients `(1, 2, 0, 0.01, 0)`;
- 95% multiplier-bootstrap intervals for the first two coefficients;
- the upper-tail, rotation-invariant test of no overall network effect.

The adjacency spectral embedding selects the three eigenvalues of largest
absolute magnitude, which is necessary for an indefinite GRDPG.

The rejection rule is `T_observed > quantile(T_bootstrap, 0.95)`, as in the
manuscript. The reported p-value is the upper-tail estimate
`(1 + sum(T_bootstrap >= T_observed)) / (1 + B)`; it is never doubled.
The +1 adjustment and sample-quantile interpolation can give slightly different
threshold decisions with finite B, so rejection rates always use the stated
critical-value rule. Existing strict-GRDPG pilots already used this rule.

## Relationship to the older SBM drafts

The new experiment is not an exact reproduction of the unidentified final SBM
run. Before the corrections accepted on 2026-09-11,
`RDPG_paralleled_bootstrap_coverage_power.R` generated graphs with `n^(-0.3)`
while its declared `nu` and labels used `n^(-0.25)`. It now consistently uses
the declared setting. Several drafts also multiplied positions built from
the unscaled B matrix by `nu^(-0.5)` when generating Y; that extra factor has
now been removed from all seven older scripts. For `P = nu Z0 Z0'`, raw ASE
has scale `sqrt(nu) Z0`, so dividing
ASE by `sqrt(nu)` targets Z0, not `Z0 / sqrt(nu)`. The new experiment generates
Y from the unscaled clean positions and divides raw ASE by the same known
sparsity factor used for graph generation. Its conventional covariates have
zero means as stated in the paper; the older code uses means `(1, 2)`.

This first implementation uses known sparsity. Empirical edge density estimates
`0.4761 * rho_n` for this kernel (and `0.47625 * rho_n` for the old SBM), not
rho_n itself. An estimated-density version must explicitly adopt a compatible
kernel/position normalization. The choice `delta=-0.3` is a proposed baseline,
not a verified reconstruction of the setting behind the submitted table.

## Running

From the project root, a small diagnostic run is:

```sh
Rscript code/revision-experiments/continuous-grdpg/strict_grdpg_bootstrap.R --mode=pilot
```

The prepared full-run configuration is:

```sh
Rscript code/revision-experiments/continuous-grdpg/strict_grdpg_bootstrap.R --mode=full
```

Full mode uses `n=3000`, 400 independent Monte Carlo repetitions per signal,
2000 bootstrap draws, four workers, and `rho_n=n^(-0.3)`. It includes
`--signals=null,weak,strong`; individual cases or subsets can be selected.
Every setting can be overridden, for example with `--n=6000`, `--nsim=200`,
`--nboot=999`, `--cores=8`, or `--delta=-0.25`.

## Sparse regimes and downsampled inference

The sparse comparison uses `delta=-0.55`.
Use the same three signal settings, with full-sample inference and one
uniformly chosen subset per m per independent graph:

```sh
Rscript code/revision-experiments/continuous-grdpg/strict_grdpg_bootstrap.R --mode=full --delta=-0.55 --msizes=50,100,200,400 --output=code/revision-experiments/continuous-grdpg/results/sparse055
```

These m values are exploratory choices, not finite-sample guarantees from the
asymptotic theorem. Full-graph ASE is computed once and reused for every subset.
The subset regression uses the same Gaussian multiplier perturbation of its
own OLS moments, normalized by m. The test statistics are
`m * ||beta_hat_z||^2` and `m * ||beta_boot_z - beta_hat_z||^2`. This is inference
for a single randomly chosen subset, not intervals formed by repeatedly fitting
many subsets of a fixed graph. No induced-subgraph embedding or Procrustes
alignment is used.

Both full and downsampled methods report conventional-coefficient coverage and
width, null rejection/type-I error, and weak/strong power. The summary groups
by sparsity, signal, and regression sample size. Full and subset results within
a repetition share the graph, response, and ASE. Subset-specific seeds ensure
the chosen nodes do not change when m ordering or bootstrap count changes.

For a multi-regime comparison, pass `--delta=-0.3,-0.55`. More extreme optional
stress tests use `--delta=-0.6,-0.7`. The m grid applies to every supplied delta;
omit `--msizes` to use full regression only. At n=3000, expected degrees are
129.28, 17.47, 11.71, and 5.26 for delta=-0.3,-0.55,-0.6,-0.7 respectively.
Regime and m choices should be recorded before the final Monte Carlo run.

Diagnostics include signed ASE eigenvalues, principal-angle subspace recovery
(without using any alignment to fit the regression), isolated-node proportion,
design condition number, and valid bootstrap draw counts. A failed inference
is retained as a row with an error and missing inference values. Summary
coverage/power/width condition on successful inference and must be read with
the displayed failure rate; valid-and-cover proportions also use all attempts.
Graph generation or eigensolver failures stop the run rather than silently
dropping networks. Weak/no signal and X independent of the latent positions can
yield valid conventional inference even when the recovered network subspace is
poor. Downsampling benefits and full-sample failure are hypotheses, not assumed
outcomes; falling outside a sufficient theorem condition does not prove failure.

## Interpretation of the manuscript's rho factor

With the fixed base matrix B0=V D0 V', generation uses P=rho Theta B0 Theta'.
Clean positions are Z0=Theta V sqrt(D0). Raw ASE has scale sqrt(rho) Z0, so
dividing raw ASE by sqrt(rho) targets Z0. If instead V Dn V' decomposes the
already-scaled matrix Bn=rho B0, the expression
rho^(-1/2) Theta V sqrt(Dn) is correct and equals Z0. The submitted paragraph
lists the fixed B0 but applies this latter normalization to its eigenvalues.
This agrees with the extra clean scaling found in the old drafts; the author's
original intention cannot be inferred from these files. The original manuscript
has not been silently edited; this distinction is saved for its later revision.

The 400-repetition Monte Carlo standard error at 95% coverage (or 5% size)
is about 0.011. Formal runs have been launched with the configuration in
FORMAL_RUN.md; launch does not imply completion. Pilot results do not determine
which settings are reported.

The script writes repetition-level CSV data, a summary CSV, and a complete RDS
object under the requested output directory. Filenames include signal cases,
regimes, the m grid, the seed, and implementation version v3. The RDS includes runtime information.
Re-running an identical configuration in the same output folder replaces that
configuration's final output after resuming completed checkpoints. Checkpoints
are saved every 20 repetitions, and mismatched settings/source fingerprints are
rejected. Use a separate `--output` directory for a changed experiment.

Historical SBM code is retained separately under `../../original-experiments/grdpg-sbm`.
