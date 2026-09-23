# School-conflict data

The individual-level data are not included in this repository.

1. Obtain version 2 of ICPSR study 37070 from
   https://doi.org/10.3886/ICPSR37070.v2.
2. Place `37070-0001-Data.rda` in this directory.
3. From the repository root, create the processed analysis object:

```sh
cd data
Rscript ../code/original-experiments/real-data/DataProcessing-all.R
cd ..
```

This creates `data/ProcessedData.Rda`, which is read by the school-conflict
analysis. Both files are ignored by Git. Users are responsible for complying
with the ICPSR terms of use and any applicable disclosure restrictions.

