# Table07

Full results for ten fixed benchmark graphs with n=500, tau=0.1, and 50 trials per graph.

Numerical results: `results/rawdata_graphs/rawdata_graph_*/RS.csv`. From the repository root:

```bash
Rscript experiments/Table07/make_table.R
```

Run methods on the saved datasets:

```bash
python experiments/Table07/run_mip.py --trial 0
```

`run_all.py` with no arguments runs all settings and trials for this result. Use `--dry-run --all-settings --all-trials` to preview the commands. `aggregate.py` combines fitted graphs and runtimes into a numerical summary. Add `--from-saved` to summarize the supplied estimates in `results/reaggregated/`.

`generate_data.py` writes new datasets under `data/generated/`.

Experiment parameters are defined in `shared/settings.py`; data formats and variable definitions are in `data/README.md`.
