# Fig01

Introductory simulation with four choices of h, sample sizes from 100 to 2000 in steps of 100, and 100 trials per setting. Figure 1 shows the results with unequal variances. The fitting scripts process both variance regimes.

Numerical results: `results/rawdata_teaser/*/RS.csv`. From the repository root:

```bash
Rscript experiments/Fig01/plot.R
```

Run methods on the saved datasets:

```bash
python experiments/Fig01/run_mip.py --trial 0 --n 500 --function "atan(x^2)"
```

`run_all.py` with no arguments runs all settings and trials for this result. Use `--dry-run --all-settings --all-trials` to preview the commands. `aggregate.py` combines fitted graphs and runtimes into a numerical summary.

`generate_data.py` writes new datasets under `data/generated/`.

Experiment parameters are defined in `shared/settings.py`; data formats and variable definitions are in `data/README.md`.
