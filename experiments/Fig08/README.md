# Fig08

Sensitivity to inverse gamma values from 3 to 11, with n=300, p=10, tau=0.2, edge probability 0.3, and 50 trials.

Numerical results: `results/rawdata_gamma/RS.csv`. From the repository root:

```bash
Rscript experiments/Fig08/plot.R
```

Run methods on the saved datasets:

```bash
python experiments/Fig08/run_mip.py --trial 0
```

`run_all.py` with no arguments runs all settings and trials for this result. Use `--dry-run --all-settings --all-trials` to preview the commands. `aggregate.py` combines fitted graphs and runtimes into a numerical summary. Add `--from-saved` to summarize the supplied estimates in `results/reaggregated/`.

For newly generated inputs, run `run_baselines.R --method pns` to estimate the superstructure. `run_all.py` includes this step. The supplied superstructures can also be used directly.

`generate_data.py` writes new datasets under `data/generated/`.

Experiment parameters are defined in `shared/settings.py`; data formats and variable definitions are in `data/README.md`.
