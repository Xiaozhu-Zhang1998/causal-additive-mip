# Fig02

Erdős–Rényi graphs with p=15 and edge probability 0.3. Each combination of n=300 or 500 and tau=0.1 or 0.2 has 50 trials.

Numerical results: `results/rawdata_c*/rawdata_p_15_030_c*/RS.csv`. From the repository root:

```bash
Rscript experiments/Fig02/plot.R
```

Run methods on the saved datasets:

```bash
python experiments/Fig02/run_mip.py --trial 0 --n 500 --p 15 --tau 0.2 --edge_prob 0.3
```

`run_all.py` with no arguments runs all settings and trials for this result. Use `--dry-run --all-settings --all-trials` to preview the commands. `aggregate.py` combines fitted graphs and runtimes into a numerical summary. Add `--from-saved` to summarize the supplied estimates in `results/reaggregated/`.

The ER experiments use the routines in `shared/`. This folder selects the settings shown in this result. `scripts/run_er_all.py` runs all 32 ER settings.

`generate_data.py` writes new datasets under `data/generated/`.

Experiment parameters are defined in `shared/settings.py`; data formats and variable definitions are in `data/README.md`.
