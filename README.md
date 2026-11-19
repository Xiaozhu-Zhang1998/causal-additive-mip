# Convex Mixed-Integer Programming for Causal Additive Models

Code, data and numerical results for **Convex Mixed-Integer Programming for Causal Additive Models with Optimization and Statistical Guarantees**. Experiment folders follow the manuscript's figure and table numbers; `docs/experiment_map.csv` identifies their inputs and rendering scripts.

## Reproduce figures and tables from saved results

With R and the packages listed in `requirements-plotting-R.txt`, run from the repository root:

```bash
Rscript --vanilla scripts/reproduce_saved.R
```

The script uses the supplied numerical results to create Figures 1–10 and Tables 1–9 under `outputs/`. To reproduce an individual figure or table, run its script directly.

```bash
Rscript experiments/Fig02/plot.R
Rscript experiments/Fig05/plot.R
Rscript experiments/Table01/make_table.R
```

Figure 5 reads `results/Expression/influence_scores.csv` directly. The optional `scripts/prepare_microarray_results.R` rebuilds this intermediate table from the supplied graph estimates.

## Rerun experiments on saved datasets

Use the Python dependencies in `requirements.txt`, the R dependencies in `docs/software.md`, and a valid Gurobi license. The recorded Python and Gurobi versions are 3.9.6 and 12.0.3, respectively.

Inputs and true graphs are included under `data/real/` and `data/simulated/`. Supplied estimates, superstructures, prior matrices, raw timings and numerical summaries are under `results/`.

For one ER (Erdos-Renyi) trial, the execution order is:

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

The supplied `RS.csv` files provide the numerical values for most figures and tables. The [interaction experiment](experiments/Table01/README.md) has its own saved summaries, graph files, and fitting script for Tables 1, 3, and 4. Re-estimation parameters for the shared fitting scripts are defined in `shared/settings.py`.

New fits write to `results/rerun/`. Existing fitted outputs require `--overwrite` to replace. `--data-root` selects another input tree and `--output-root` selects another output directory. Preview commands without fitting:

```bash
python experiments/Fig02/run_mip.py --trial 15 --dry-run
python experiments/Fig02/run_all.py --dry-run --all-settings --all-trials
```

A result folder's `run_all.py` with no arguments runs all settings and trials for that figure or table and may take days. Figure 2 selects four sample/noise settings, Figure 3 selects four graph settings, and each ER table selects eight graph settings. `scripts/run_er_all.py` covers all 32 ER settings. The sequential runners call the same functions as the commands for individual trials.

To plot rerun summaries, set `JASA_RESULTS_ROOT` to their result tree. See `docs/execution.md` for execution details and runtime normalization, and `data/README.md` for data formats and variable definitions.

## Generate new inputs and repeat the analyses

Use the included DGP code to generate new synthetic datasets, or prepare new analysis inputs from the original real data. Install the Python and R dependencies listed in `requirements.txt`, `requirements-plotting-R.txt` and `docs/software.md`, and use a valid Gurobi license for MIP fitting. Run the commands from the repository root. New inputs are written under `data/generated/`, and fitted estimates and summaries under `results/rerun/`.

### Synthetic simulations

For the complete Figure 2 experiment, generate the four sample-size/noise settings and their 50 trials, then fit the methods, aggregate the estimates, and plot the results:

```bash
python experiments/Fig02/generate_data.py --all-settings --all-trials
python experiments/Fig02/run_all.py --data-root data/generated --all-settings --all-trials
JASA_RESULTS_ROOT=results/rerun Rscript --vanilla experiments/Fig02/plot.R
```

The fitting runner performs the R baseline methods and neighborhood selection, constructs bootstrap priors where used, fits the Python baselines and MIP, and aggregates the results. For one setting and trial, replace `--all-settings --all-trials` in both Python commands with `--trial 0 --n 500 --p 15 --tau 0.2 --edge_prob 0.3`; plotting the complete figure requires completing its full grid. Full experiment runs may take hours or days.

Other simulation folders provide their own `generate_data.py` scripts and experiment settings. Figure 7 reuses ER datasets generated through Figure 3. The [interaction experiment](experiments/Table01/README.md) has a separate runner that generates observations and fits the two MIP estimators.

### Causal Chamber

Prepare the analysis subsamples from the included original `data/real/causal_chamber/X.csv`, run all sample sizes and trials, and plot the new results:

```bash
python experiments/Fig04/prepare_data.py
python experiments/Fig04/run_all.py --data-root data/generated --all-settings --all-trials
JASA_RESULTS_ROOT=results/rerun Rscript --vanilla experiments/Fig04/plot.R
```

### Microarray

Preprocess the included original gene expression data, generate 30 new samples with 118 observations drawn with replacement, fit the methods, calculate the influence scores, and plot:

```bash
Rscript --vanilla experiments/Fig05/prepare_data.R
python experiments/Fig05/resample.py --input data/generated/microarray/Arabidopsis_MVA_standardized.csv
python experiments/Fig05/run_all.py --data-root data/generated --all-trials
JASA_RESULTS_ROOT=results/rerun Rscript --vanilla scripts/prepare_microarray_results.R
JASA_RESULTS_ROOT=results/rerun Rscript --vanilla experiments/Fig05/plot.R
```

The new samples are saved under `data/generated/Expression/subsamples/`. The fitting runner first estimates their superstructures, then fits the baseline methods and MIP. Re-estimation settings, including the microarray MIP penalty of 0.08, are defined in `shared/settings.py`.

These workflows follow the paper's experimental settings using newly generated inputs.

## Repository layout

- `experiments/`: execution and rendering scripts by figure or table.
- `shared/`: DGPs, MIP algorithms, baseline wrappers, evaluation and plotting functions.
- `third_party/`: external implementations and their import dependencies.
- `data/`: supplied real/synthetic datasets and true graphs.
- `results/`: supplied numerical results, estimates, timings and intermediate results; new fits use `results/rerun/`.
- `docs/`: experiment map, software information and execution guide.

Provided data and results are tracked in Git. Generated `outputs/`, `data/generated/`, `results/rerun/` and `results/reaggregated/` are ignored. Run the scripts directly from the repository root.

Software credits and component licenses are listed in `third_party/README.md` and `LICENSE.md`.
