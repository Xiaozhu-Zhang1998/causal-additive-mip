# Table 4: structural intervention distance

Table 4 reports mean structural intervention distance (SID) and its standard error for both noise levels in the [interaction study](../Table01/README.md). It reads the saved EDF BIC summaries for the two MIP estimators and computes SID for the four baseline methods from their saved DAGs.

With the R package SID installed, run from the repository root:

```bash
Rscript --vanilla experiments/Table04/make_table.R
```

The script writes `outputs/tables/Table04.csv` and `outputs/tables/Table04.tex`. The Table 1 README gives commands for regenerating datasets and fitting methods.
