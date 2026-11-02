# Fig10

Sensitivity to Student-t noise with n=500, p=10, tau=0.2, edge probability 0.3, and 50 trials. The plotted degrees of freedom are 10, 15, 20, and 25.

Numerical results: `results/rawdata_tdisn/*/RS.csv`. From the repository root:

```bash
Rscript experiments/Fig10/plot.R
```

Run methods on the saved datasets:

```bash
python experiments/Fig10/run_mip.py --trial 0 --df 10
```

`run_all.py` with no arguments runs all settings and trials for this result. Use `--dry-run --all-settings --all-trials` to preview the commands. `aggregate.py` combines fitted graphs and runtimes into a numerical summary. Add `--from-saved` to summarize the supplied estimates in `results/reaggregated/`.

`generate_data.py` writes new datasets under `data/generated/`.

Experiment parameters are defined in `shared/settings.py`; data formats and variable definitions are in `data/README.md`.
