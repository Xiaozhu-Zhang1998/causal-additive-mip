# Data inputs and formats

The datasets, true graphs and numerical results are included with this repository. Trial indices start at zero.

| Location | Contents |
|---|---|
| `simulated/rawdata_*/` | Synthetic observations `X_l.csv` and true graphs `dag_oracle_l.csv`; collection names match the result folders. |
| `simulated/interaction/` | 520 datasets with 300 observations and eight variables, plus true edges and interaction functions in `ground_truth.csv`. See [Tables 1 and 3](../experiments/Table01/README.md). |
| `real/causal_chamber/X.csv` | Original 10000-by-8 observational matrix. |
| `real/causal_chamber/subsamples/` | Standardized analysis datasets and true graphs for 70 combinations of sample size and trial. |
| `real/microarray/` | Original expression source, raw and standardized 118-by-15 matrices. |
| `real/microarray/subsamples/` | Thirty bootstrap datasets, each with 118 rows sampled with replacement. |

Each observation matrix has n rows and p numeric columns. Synthetic columns identify variables 1 through p. Adjacency entry `[i,j]=1` represents the directed edge i to j, using the same variable order as the observation matrix. Files are headerless: Python matrices use commas; R matrices use whitespace. The shared readers handle both formats.

Causal Chamber columns are R, G, B, Ce, Ie1, Ie2, Ie3, and theta1-minus-theta2: the red/green/blue light inputs, the C and three I measurements, and the difference between polarization angles. The analysis values are standardized and dimensionless. `experiments/Fig04/prepare_data.py` standardizes the full input matrix and selects observations at fixed positions in consecutive blocks of ten. The data source is the `uniform_reference` environment of `lt_interventions_standard_v1` in [Causal Chamber](https://github.com/juangamella/causal-chamber).

Microarray columns are AACT1, AACT2, HMGS, HMGR1, HMGR2, MK, MPDC1, MPDC2, IPPI2, FPPS1, FPPS2, DPPS1, DPPS3, GGPPS3, and GGPPS4. Raw values are gene expression measurements; standardized values are dimensionless. `experiments/Fig05/prepare_data.R` selects this order and applies R `scale` over the 118 observations. The original [expression source file](https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fgb-2004-5-11-r92/MediaObjects/13059_2004_896_MOESM1_ESM.txt) is included locally. `resample.py` creates 30 bootstrap datasets.

Estimated graphs, matrices of candidate parents, prior matrices and runtimes are under `results/`. `dag_<method>_l.csv` is an estimated graph; `dag_superstructure_l.csv` permits candidate parents; `dag_CAMpartial_l.csv` and `dag_NPVARpartial_l.csv` contain prior edge selection probabilities. `time_<method>_l.txt` stores elapsed seconds. `RS.csv` contains numerical summaries. Figure 5 reads `results/Expression/influence_scores.csv`, calculated from `results/Expression/subsamples/`.

SHD counts binary adjacency differences, with two units for a reversed edge. SID uses `SID::structIntervDist`; directed F1, Chamber TDR and FDR use directed edge comparisons. Microarray score definitions are in `experiments/Fig05/README.md`.

New inputs go to `data/generated/`; new estimates and summaries go to `results/rerun/` or `results/reaggregated/`. The optional `scripts/build_data_archive.py` creates an archive of the supplied inputs.
