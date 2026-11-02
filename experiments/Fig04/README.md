# Fig04

Causal Chamber analysis with sample sizes from 200 to 800 in steps of 100 and 10 trials per sample size.

Numerical results: `results/causal_chamber/RS.csv`. From the repository root:

```bash
Rscript experiments/Fig04/plot.R
```

Run methods on the saved datasets:

```bash
python experiments/Fig04/run_mip.py --trial 0 --n 800
```

`run_all.py` with no arguments runs all settings and trials for this result. Use `--dry-run --all-settings --all-trials` to preview the commands. `aggregate.py` combines fitted graphs and runtimes into a numerical summary. Add `--from-saved` to summarize the supplied estimates in `results/reaggregated/`.

`prepare_data.py` creates the analysis subsamples from `data/real/causal_chamber/X.csv`.

Experiment parameters are defined in `shared/settings.py`; data formats and variable definitions are in `data/README.md`.
