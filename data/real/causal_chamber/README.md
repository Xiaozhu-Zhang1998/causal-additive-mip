# Causal Chamber inputs

`X.csv` contains the original 10000-by-8 observational matrix. `subsamples/X_l.csv` and `subsamples/dag_oracle_l.csv` contain the analysis datasets and true graphs. Estimates, summaries and timings are in `results/causal_chamber/`.

Columns are R, G, B, Ce, Ie1, Ie2, Ie3, and theta1-minus-theta2. `experiments/Fig04/prepare_data.py` standardizes the full matrix and selects observations at fixed positions within blocks to form the subsamples. Variable definitions and file formats are in `data/README.md`.
