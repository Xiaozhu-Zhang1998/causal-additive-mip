# Table 3: interaction effects at tau = 0.2

This table uses the same experiment as [Table 1](../Table01/README.md), at noise heterogeneity `tau = 0.2`. From the repository root:

```bash
Rscript --vanilla experiments/Table03/make_table.R
```

For the corresponding MIP runs, pass `--noise-taus 0.2` to `experiments/Table01/run_interactions.py`. For the baselines, pass `--tau 0.2` to `experiments/Table01/run_baselines.py`.
