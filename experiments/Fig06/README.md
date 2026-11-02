# Fig06

Early stopping on the InsuranceSmall graph with p=15, n=500, and 50 trials. Relative gap thresholds range from 0 to 0.1 in steps of 0.01.

Numerical results: `results/rawdata_earlystop/RS.csv`. From the repository root:

```bash
Rscript experiments/Fig06/plot.R
```

Run methods on the saved datasets:

```bash
python experiments/Fig06/run_mip.py --trial 0
```

`run_all.py` with no arguments runs all settings and trials for this result. Use `--dry-run --all-settings --all-trials` to preview the commands. `aggregate.py` combines fitted graphs and runtimes into a numerical summary. Add `--from-saved` to summarize the supplied estimates in `results/reaggregated/`.

The set of candidate parents is the true moral graph. The displayed stopping thresholds range from 0 through 0.1, with 50 trials.

`generate_data.py` writes new datasets under `data/generated/`.

Experiment parameters are defined in `shared/settings.py`; data formats and variable definitions are in `data/README.md`.
