# Reproducibility guide

Run all commands from the repository root. Output directories must not already
exist unless a script explicitly supports checkpoint resumption.

## 1. Original manuscript experiments

The scripts used for the submitted manuscript are retained under
`code/original-experiments`. They are organized by analysis family:

- `local-subgraph`: corrected and uncorrected local-subgraph regression,
  bootstrap inference, population targets, and misspecification experiments;
- `grdpg-sbm`: the original SBM/GRDPG bootstrap and downsampling experiments;
- `downsampling`: local-average and local-transitivity experiments; and
- `real-data`: data processing and the school-conflict analysis.

These are archival scripts and generally run from their own directory. Several
write multiple large `.rds` files. Read the script-level settings before
starting a formal run. The revision experiments below provide cleaner command
line interfaces and explicit output directories.

## 2. Estimator comparison added during revision

The formal design has `n=6000`, 400 repetitions per sparsity regime, and 999
structural-bootstrap draws. A quick smoke run is:

```sh
cmp=code/revision-experiments/estimator-comparison
Rscript "$cmp/coverage_experiment.R" "$cmp/results/smoke" 80 2 49 -0.25
```

Run the four formal regimes separately:

```sh
cmp=code/revision-experiments/estimator-comparison
Rscript "$cmp/coverage_experiment.R" "$cmp/results/main_n6000_r400_delta_m055" 6000 400 999 -0.55
Rscript "$cmp/coverage_experiment.R" "$cmp/results/main_n6000_r400_delta_m040" 6000 400 999 -0.40
Rscript "$cmp/coverage_experiment.R" "$cmp/results/main_n6000_r400_delta_m030" 6000 400 999 -0.30
Rscript "$cmp/coverage_experiment.R" "$cmp/results/main_n6000_r400_delta_m025" 6000 400 999 -0.25
Rscript "$cmp/combine_and_plot.R" "$cmp/results"
Rscript "$cmp/validate_results.R" "$cmp/results/main_n6000_r400_combined"
```

## 3. Aggregation over fixed disjoint subsets

This experiment replays the data seeds from the estimator comparison, so
Section 2 must be completed first. Unit tests can be run independently:

```sh
agg=code/revision-experiments/disjoint-aggregation
Rscript "$agg/test_disjoint.R"
Rscript "$agg/test_graph_counts.R"
```

Formal runs and result preparation:

```sh
agg=code/revision-experiments/disjoint-aggregation
Rscript "$agg/run_disjoint.R" "$agg/results/formal_n6000_r400_delta_m055" 400 -0.55
Rscript "$agg/run_disjoint.R" "$agg/results/formal_n6000_r400_delta_m040" 400 -0.40
Rscript "$agg/run_disjoint.R" "$agg/results/formal_n6000_r400_delta_m030" 400 -0.30
Rscript "$agg/run_disjoint.R" "$agg/results/formal_n6000_r400_delta_m025" 400 -0.25
Rscript "$agg/prepare_results.R"
```

For each `m`, the deterministic blocks are
`I_k = {(k-1)m+1, ..., km}`. Network features are computed from the full graph;
only the regression observations are partitioned.

## 4. Continuous GRDPG experiment

A small diagnostic run is:

```sh
grdpg=code/revision-experiments/continuous-grdpg
Rscript "$grdpg/strict_grdpg_bootstrap.R" \
  --mode=pilot --n=300 --nsim=4 --nboot=99 --cores=1 \
  --output="$grdpg/results/smoke"
```

The formal experiment reported in the revision is:

```sh
grdpg=code/revision-experiments/continuous-grdpg
Rscript "$grdpg/strict_grdpg_bootstrap.R" \
  --mode=full --n=3000 --nsim=400 --nboot=2000 --delta=-0.3 \
  --msizes=none --signals=null,weak,strong --seed=2601 --cores=2 \
  --output="$grdpg/results/formal_n3000_b2000_r400/dense"
Rscript "$grdpg/build_dense_summary.R" \
  "$grdpg/results/formal_n3000_b2000_r400/dense" \
  "$grdpg/results/formal_n3000_b2000_r400/dense-summary"
```

## 5. School-conflict analysis and dimension sensitivity

First follow `data/README.md` to create `data/ProcessedData.Rda`. Then run:

```sh
real=code/revision-experiments/real-data-dimension-sensitivity
for d in 2 3 4; do
  Rscript "$real/global_embedding_response1.R" \
    --dimension="$d" --bootstrap-reps=500 --seed=20260911 \
    --data-file=data/ProcessedData.Rda \
    --results-dir="$real/results/dimension-sensitivity"
done
Rscript "$real/summarize_dimension_sensitivity.R" \
  "$real/results/dimension-sensitivity"
```

The same students and schools are retained across dimensions. Each school is
embedded separately, using its own estimated sparsity, before the school-level
embeddings are assembled into the global block design.

## 6. Included reported results

The checked-in files under `results/reported` allow the manuscript tables and
figures to be inspected without rerunning the expensive experiments. They are
summaries rather than substitutes for the repetition-level audit files created
by the formal commands above.

