# Convex Mixed-Integer Programming for Causal Additive Models

Code, data and numerical results for **Convex Mixed-Integer Programming for Causal Additive Models with Optimization and Statistical Guarantees**. Experiment folders follow the manuscript's figure and table numbers; `docs/experiment_map.csv` identifies their inputs and rendering scripts.

## Reproduce figures and tables from saved results

With R and the packages listed in `requirements-plotting-R.txt`, run from the repository root:

```bash
Rscript --vanilla scripts/reproduce_saved.R
```

This creates Figures 1–10 and Tables 2–7 under `outputs/`. It reads the supplied numerical results using R. Individual results can be rendered directly:

```bash
Rscript experiments/Fig02/plot.R
Rscript experiments/Fig05/plot.R
Rscript experiments/Table03/make_table.R
```

Figure 5 reads `results/Expression/influence_scores.csv` directly. The optional `scripts/prepare_microarray_results.R` rebuilds this intermediate table from the supplied graph estimates.

## Rerun experiments on saved datasets

Use the Python dependencies in `requirements.txt`, the R dependencies in `docs/software.md`, and a valid Gurobi license. The recorded Python and Gurobi versions are 3.9.6 and 12.0.3, respectively.

Inputs and true graphs are included under `data/real/` and `data/simulated/`. Supplied estimates, superstructures, prior matrices, raw timings and numerical summaries are under `results/`.

For one ER trial, the execution order is:

```bash
# R baseline methods and preliminary neighborhood selection.
Rscript experiments/Fig02/run_baselines.R --trial 15 --n 500 --p 15 --tau 0.2 --edge_prob 0.3

# Bootstrap priors for settings that use them.
Rscript experiments/Fig02/run_bootstrap.R --trial 15 --n 500 --p 15 --tau 0.2 --edge_prob 0.3

# Proposed MIP method and an independently selected Python baseline.
python experiments/Fig02/run_mip.py --trial 15 --n 500 --p 15 --tau 0.2 --edge_prob 0.3
python experiments/Fig02/run_baselines.py --method notears --trial 15 --n 500 --p 15 --tau 0.2 --edge_prob 0.3

# Aggregate completed estimates and their recorded runtimes.
python experiments/Fig02/aggregate.py --n 500 --p 15 --tau 0.2 --edge_prob 0.3
```

Trial indices start at zero. MIP can use the included superstructures and priors directly; newly computed intermediates take precedence. ER settings with p=10 and edge probability 0.3 or 0.4 do not use bootstrap priors; the remaining ER settings do.

The supplied `RS.csv` files provide the numerical values for the manuscript figures and tables. Re-estimation parameters are defined in `shared/settings.py`.

New fits write to `results/rerun/`. Existing fitted outputs require `--overwrite` to replace. `--data-root` selects another input tree and `--output-root` selects another output directory. Preview commands without fitting:

```bash
python experiments/Fig02/run_mip.py --trial 15 --dry-run
python experiments/Fig02/run_all.py --dry-run --all-settings --all-trials
```

A result folder's `run_all.py` with no arguments runs all settings and trials for that figure or table and may take days. Figure 2 selects four sample/noise settings, Figure 3 selects four graph settings, and each ER table selects eight graph settings. `scripts/run_er_all.py` covers all 32 ER settings. The sequential runners call the same functions as the commands for individual trials.

To plot rerun summaries, set `JASA_RESULTS_ROOT` to their result tree. See `docs/execution.md` for execution details and runtime normalization, and `data/README.md` for data formats and variable definitions.

## Repository layout

- `experiments/`: execution and rendering scripts by figure or table.
- `shared/`: DGPs, MIP algorithms, baseline wrappers, evaluation and plotting functions.
- `third_party/`: external implementations and their import dependencies.
- `data/`: supplied real/synthetic datasets and true graphs.
- `results/`: supplied numerical results, estimates, timings and intermediate results; new fits use `results/rerun/`.
- `docs/`: experiment map, software information and execution guide.

Provided data and results are tracked in Git. Generated `outputs/`, `data/generated/`, `results/rerun/` and `results/reaggregated/` are ignored. Run the scripts directly from the repository root.

Software credits and component licenses are listed in `third_party/README.md` and `LICENSE.md`.
