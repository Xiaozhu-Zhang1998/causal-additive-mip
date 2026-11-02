# Experiment execution

`main('teaser','all',experiment='Fig01')` enters the Figure 1 workflow through `shared/cli.py`. The `teaser` family selects the introductory experiment:

```python
if family == 'teaser':
    return teaser_main(a, action)
```

The call runs `teaser_main(a,'all')`. For ER experiments, the `experiment` argument selects the settings associated with a particular figure or table.

## Figure 1 execution sequence

For a concrete example, select h=atan(x²), n=500, trial 15. The code does the following, in order:

| Step | File and function | Work performed |
|---|---|---|
| 1 | `shared/cli.py`: `main` → `teaser_main` | Reads command-line arguments and enters the requested loops over functions, sample sizes, and trials. |
| 2 | `shared/special_fits.py`: `fit_teaser` | Reads the two saved unequal/equal-noise datasets. For each, fits nonlinear MIP and linear MIP, then saves graph estimates, exact recovery, SHD and elapsed time in `python_500_15.csv`. |
| 2a | `shared/mip/basis.py`: `gen_basis_feat` | Constructs the spline basis for nonlinear MIP or centered linear features for linear MIP. |
| 2b | `shared/mip/cutting_plane.py`: `solve_sem_DAG_OA_withinitial1` | Runs the Gurobi cutting-plane MIP solver. This is called four times for the example: two basis choices × two noise regimes. |
| 3 | `Rscript shared/baselines/teaser.R`, with mode=`fit` | Reads the same datasets. Calls `CAM::CAM(X,pruning=TRUE)` and `prune(X,NPVAR(X))`, using `third_party/NPVAR/NPVAR.R` and `utils.R`. Saves graphs and recovery/SHD/time rows in `R_500_15.csv`. |
| 4 | `shared/aggregation.py`: `aggregate_teaser` | After all requested loops finish, joins the `python_*.csv` and `R_*.csv` batches into an `RS.csv` for each h function.  |

Inputs for this example are:

```text
 data/simulated/rawdata_teaser/atan(x^2)/unequal_500_15.csv
 data/simulated/rawdata_teaser/atan(x^2)/equal_500_15.csv
```

Outputs are under `results/rerun/rawdata_teaser/atan(x^2)/`. Use `generate_data.py` to generate new data and `plot.R` to draw Figure 1.

## What “all” means

`'all'` selects the Python MIP fits, R baseline fits and CSV aggregation. The `--all-settings` and `--all-trials` flags select the full experiment grid.

When `experiments/Fig01/run_all.py` is run with no arguments, it adds `--all-settings --all-trials`. Those flags select four h functions, n=100,200,…,2000, and trials 0–99. If `main` were called without those flags, the parser defaults instead select atan(x²), n=500, trial 0.

This is an equivalent outline of the teaser branch:

```python
for function in selected_functions:
    for n in selected_sample_sizes:
        for trial in selected_trials:
            fit_teaser(trial, n, function, ...)
            subprocess.run(['Rscript', '--vanilla', '.../teaser.R', ...])
aggregate_teaser(args)
```

Preview the execution sequence from the repository root:

```bash
python experiments/Fig01/run_all.py --n 500 --trial 15 --dry-run
```

The dry run prints the names of the fitting functions, solvers, R scripts and aggregation functions.

## Other experiments

For ER experiments, `shared/cli.py` instead calls `select_specs`, then the R wrapper for neighborhood selection and baseline methods, optional bootstrap priors, Python baseline wrappers, `shared/fitting.py:fit_mip`, and aggregation. Figure 3's full batch is restricted to its four graph settings ×50 trials. The complete ER grid of 32 settings has an explicit runner at `scripts/run_er_all.py`.

## Results and inputs

All provided numerical results and intermediates are in one `results/` tree. `results/Expression/influence_scores.csv` is the Figure 5 metric table, calculated from all 30 trials in `results/Expression/subsamples/`. The input files in `data/real/microarray/subsamples/` are bootstrap samples with 118 rows drawn with replacement. Default plotting reuses the influence table; the optional `scripts/prepare_microarray_results.R` can rebuild it if requested. New estimation runs write under `results/rerun/` to preserve provided numbers.

`data/real/causal_chamber/X.csv` is the original full Chamber input dataset; `X_l.csv` are its analysis subsamples.

The included datasets and true graphs are under `data/simulated/` or `data/real/`; the fitting functions resolve them through `shared/paths.py:source_folder`. Provided estimated graphs, superstructures, priors and timings are separate files under `results/`. `auxiliary_file` prefers newly fitted intermediates in `results/rerun/`, then the supplied intermediates associated with the input tree. For an alternative data root, intermediate files are read from that input tree or its selected output directory. For example, a default Figure 2 MIP trial 15 reads `data/simulated/rawdata_c3/rawdata_p_15_030_c3/X_15.csv`, its saved PNS/prior matrices from `results/rawdata_c3/rawdata_p_15_030_c3/`, and writes new estimates into the corresponding `results/rerun/` folder.

`aggregate.py --from-saved` reads provided estimates/timings from `results/` and truth from `data/`, writing a new summary under `results/reaggregated/`. `--results-root` can select an alternative estimate tree for this aggregation mode.


## Numerical summaries and runtimes

Summary columns `n`, `p`, `tau` and `edge_density` identify sample size, variable count, noise heterogeneity and ER edge probability. In early-stopping summaries, `tau` identifies the stopping threshold.

Aggregation records `time_total`, `tuning_configurations` and `time_per_configuration`. Table entries use elapsed time divided by the number of tuning configurations: MIP uses 6 for ER p=10/rho=0.3 or 0.4, and 12 for the remaining ER and fixed-graph settings; NOTEARS/GSGES use 3. Prior generation uses 20 resamples, and gamma sensitivity uses 8 lambda values per gamma. Microarray MIP uses one lambda value. Figure 6 and Table 2 report unscaled elapsed times.

The Table 2 comparison adds p*(1+log(2)) to the projected cutting-plane objective before comparing it with the vanilla objective. Its true edges are defined in the Python comparison driver. In solver calls, gamma is inverse ridge strength (5 corresponds to coefficient 0.2), and `lam2` is the edge penalty passed directly to the optimization.

The synthetic data and saved trial inputs can be reused directly. Stochastic fitting procedures use their existing randomization; rerun graph estimates can vary. Each experiment's README identifies its inputs and terminal commands.
