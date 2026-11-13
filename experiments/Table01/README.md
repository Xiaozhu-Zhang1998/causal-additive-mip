# Table 1: interaction effects

Table 1 reports the results for noise heterogeneity `tau = 0.1`. Supplementary Table 3 reports `tau = 0.2` from the same experiment. There are 13 settings with 20 trials each, eight variables, and 300 observations per trial. See Section 8 and Appendix J.2–J.3 for the model and simulation design.

From the repository root, create the tables as CSV and LaTeX files under `outputs/tables/`:

```bash
Rscript --vanilla experiments/Table01/make_table.R
Rscript --vanilla experiments/Table03/make_table.R
```

The tables use the EDF BIC results for MIP-main and MIP-interaction and calculate directed SHD from the saved baseline DAGs. Directed SHD counts one for each differing adjacency entry, so an edge reversal counts as two. Each cell reports the mean and standard error across 20 trials.

## Data and saved results

- `data/simulated/interaction/tau_*/setting_*/X_<trial>.csv` contains the observations. Settings (0–12), trials (0–19), and node indices in the graph files start at zero. Matrix columns follow the manuscript order X1 through X8.
- `data/simulated/interaction/ground_truth.csv` lists true edges and interaction triples with their function assignments. In the pure interaction setting, an edge can have a blank main effect. For the ratio interaction, `reverse_inputs` is true when `parent2` is the numerator.
- `results/interaction/mip_results.csv` records the two MIP estimators under both EDF BIC and nominal OLS BIC for each dataset, including SHD, SID, solver status, and method time in seconds. `mip_summary.csv` reports means and standard errors by setting, noise level, method, and BIC criterion.
- `results/interaction/mip_dags.npz` stores `true_adjacency` and `estimated_adjacency`, each with shape (2080, 8, 8). The first axis follows `graph_number` in the CSV. Use `numpy.load(..., allow_pickle=False)` to read it. Baseline adjacency matrices and timings are under `results/interaction/tau_*/setting_*/`.

To recompute the MIP graph metrics and summary from the saved estimates, run:

```bash
python experiments/Table01/run_interactions.py --summarize-results results/interaction
```

This writes the recalculated CSVs under `results/reaggregated/interaction/`.

## Rerun the experiment

Use the Python packages in the root `requirements.txt`. The paper settings use Gurobi. The script also supports `--backend scipy`.

```bash
# Regenerate one dataset and its ground truth using the saved seed.
python experiments/Table01/run_interactions.py --seed 20261002 --settings 0 --noise-taus 0.1 --trial 0 --save-data-only

# Fit both MIP estimators for that trial, using the full penalty grid.
python experiments/Table01/run_interactions.py --seed 20261002 --settings 0 --noise-taus 0.1 --trial 0

# Fit all four baselines on the corresponding supplied observations.
python experiments/Table01/run_baselines.py --setting 0 --tau 0.1 --trial 0
```

To fit both MIP estimators on all 520 datasets, run the script with `--seed 20261002` and omit `--settings`, `--noise-taus`, and `--trial`. The script uses eight threads by default and saves checkpoints, selected models, and summaries under `results/rerun/interaction/mip/`. Use `--output-dir` to select a separate run directory. The full batch can take days. Generated inputs go to `data/generated/interaction/`. The baseline command uses the shared Python fitting code and writes under `results/rerun/interaction/tau_*/setting_*/`.

Both MIP estimators start from all directed edges between distinct nodes and screen candidate supports using the continuous relaxation. Both ridge penalties equal one. The interaction penalty is the edge penalty multiplied by `rho`.
