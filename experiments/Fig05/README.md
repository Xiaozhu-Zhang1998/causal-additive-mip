# Fig05

The microarray analysis uses 30 bootstrap datasets (trials 0–29), each containing 118 rows and 15 genes. Rows are drawn with replacement from the standardized matrix with 118 observations. Both the input and result folders are named `subsamples`.

From the repository root, reproduce the figure using the supplied score table:

```bash
Rscript experiments/Fig05/plot.R
```

To recalculate that table from the supplied graph estimates:

```bash
Rscript scripts/prepare_microarray_results.R
```

The left panel shows `(correct-reversed)/max(edges,1)`; the right panel shows `correct/max(edges,1)`. Correct and reversed counts refer to reachable gene pairs among the 91 known pathway orderings, whereas `edges` is the number of estimated directed edges. Scores can exceed 1. Gene names, column order, pathway edges and scoring are defined in `shared/pathway.R`.

Fit a single saved dataset using the common baseline and MIP pipeline:

```bash
Rscript experiments/Fig05/run_baselines.R --trial 0
python experiments/Fig05/run_baselines.py --trial 0 --method notears
python experiments/Fig05/run_mip.py --trial 0
python experiments/Fig05/aggregate.py
```

The R command fits CAM, NPVAR, EqVar and RESIT and saves the gamboost superstructure. Python baseline methods are NOTEARS, GLOBE, TOPIC and GSGES; select one with `--method`, or omit it to run all four. MIP uses lambda 0.08, inverse gamma 5 (ridge coefficient 0.2), cubic splines with three internal knots, the estimated superstructure, no bootstrap priors, relative MIP gap 0, and a time limit of 120 seconds per optimization call. Included superstructures allow MIP to run directly.

New estimates and summaries go to `results/rerun/Expression/subsamples/`. To plot their aggregated `RS.csv`, set `JASA_RESULTS_ROOT` to `results/rerun`. `run_all.py` with no arguments runs all 30 trials; use `--dry-run --all-trials` to preview it.

`prepare_data.R` selects the shared gene order and standardizes the original data into `data/generated/microarray/`. `resample.py` generates 30 new bootstrap inputs under `data/generated/Expression/subsamples/`. To use newly standardized values:

```bash
Rscript experiments/Fig05/prepare_data.R
python experiments/Fig05/resample.py --input data/generated/microarray/Arabidopsis_MVA_standardized.csv
```

The supplied estimates correspond to the included bootstrap datasets. `resample.py` generates new draws using NumPy. To fit newly generated datasets, pass `--data-root data/generated` to each fitting command; first create their superstructures with `run_baselines.R`.
