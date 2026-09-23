# Original manuscript experiments

This directory preserves the R scripts used during development of the
original manuscript. They are grouped by analysis family and retain the
original seeds and simulation settings.

The scripts predate the command-line wrappers used for the revision. Run them
from their own directory and inspect their output filenames before execution.
They may require substantial memory and may write multiple large `.rds` files.

The primary files corresponding to the manuscript analyses are:

- `local-subgraph/Subgraph_count_bootstrap_samples_Corrected_Noncorrected.R`:
  corrected versus uncorrected local-subgraph regression;
- `local-subgraph/Bootstrap_LS_misspecification.R`: projection inference under
  outcome-model misspecification;
- `grdpg-sbm/RDPG_parallel_boostrap_coverage_power.R`: SBM/GRDPG bootstrap
  coverage and network-effect testing;
- `downsampling/downsample_local_average_200_rep.R` and
  `downsampling/downsample_local_trans_200_rep.R`: downsampling examples; and
- `real-data/DataProcessing-all.R` together with the regression/embedding
  scripts: school-conflict analysis.

The revision-specific, publication-facing implementations are documented
separately under `../revision-experiments`.

