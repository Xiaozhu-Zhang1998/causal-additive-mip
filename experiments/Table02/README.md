# Table02

Comparison of the cutting plane and vanilla MIP formulations on five fixed graphs, using n=500 and 50 trials. Inverse gamma is 1 and lambda is 0.005.

Numerical results: `results/rawdata_vanillacomp/*/RS.csv`. From the repository root:

```bash
Rscript experiments/Table02/make_table.R
```

Run methods on the saved datasets:

```bash
python experiments/Table02/run_comparison.py --trial 0
```

`run_all.py` with no arguments runs all settings and trials for this result. Use `--dry-run --all-settings --all-trials` to preview the commands. `aggregate.py` combines fitted graphs and runtimes into a numerical summary.

`generate_data.py` writes new datasets under `data/generated/`.

Experiment parameters are defined in `shared/settings.py`; data formats and variable definitions are in `data/README.md`.
