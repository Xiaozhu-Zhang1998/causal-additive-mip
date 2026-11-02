# Software and computation

The verified environment uses Python 3.9.6, Gurobi 12.0.3 and R 4.5.0. Activate the Python environment before running the terminal scripts. Install the recorded Python dependencies with:

```bash
python -m pip install -r requirements.txt
```

Gurobi estimation requires a valid license. Rendering figures and tables from saved results uses R and the packages listed in `requirements-plotting-R.txt`.

| R packages | Verified versions, in the same order |
|---|---|
| ggplot2, dplyr, tidyr, purrr | 4.0.3, 1.2.1, 1.3.2, 1.2.2 |
| readr, stringr, patchwork, ggh4x | 2.1.5, 1.6.0, 1.3.2, 0.3.1 |
| igraph, SID | 2.3.4, 1.1 |
| CAM, mgcv, mboost | 1.0, 1.9.3, 2.9.11 |
| np, glmnet, earth | 0.60.18, 5.1, 5.3.5 |

In R, install the plotting packages with:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "purrr", "readr",
                   "stringr", "patchwork", "ggh4x", "igraph"))
```

Estimation also uses CAM, mgcv, mboost, np, glmnet and earth; graph aggregation uses SID. Archived releases are available through the [CAM archive](https://cran.r-project.org/src/contrib/Archive/CAM/) and [SID package page](https://cran.r-project.org/package=SID). External baseline sources and notices are listed in `third_party/README.md`.

GLOBE uses rpy2, R and earth. Its Python interpreter and R shared library must use matching CPU architectures, such as arm64 on macOS. R terminal adapters call `python3`; set `JASA_PYTHON` to select another interpreter. Python entry points include `third_party/` on the import path.

MIP runners use three joblib workers by default; `--n-jobs` changes this. Gurobi and BLAS thread defaults are retained. Time limits apply per optimization call. Figure 1 has no solver time limit. Full experiment batches can take hours or days. You can run individual trials or render figures and tables from saved results separately.
