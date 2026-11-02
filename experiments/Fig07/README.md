# Fig07

Comparison of spline and radial bases with n=500, tau=0.1, and 50 trials for each of four graph settings. The radial specification uses six Gaussian basis functions with bandwidth 0.5.

Numerical results: `results/rawdata_newbasis/*/RS.csv`. From the repository root:

```bash
Rscript experiments/Fig07/plot.R
```

Run methods on the saved datasets:

```bash
python experiments/Fig07/run_mip_radial.py --trial 0
```

`run_all.py` with no arguments runs all settings and trials for this result. Use `--dry-run --all-settings --all-trials` to preview the commands. `aggregate.py` combines fitted graphs and runtimes into a numerical summary. Add `--from-saved` to summarize the supplied estimates in `results/reaggregated/`.

Experiment parameters are defined in `shared/settings.py`; data formats and variable definitions are in `data/README.md`.
