# External methods

The following external implementations are used by the experiment scripts. Shared wrappers are in `shared/baselines/`.

| Component | Original source | Implementation and license |
|---|---|---|
| NOTEARS | https://github.com/xunzheng/notears | Nonlinear implementation; LocallyConnected/LBFGSB/trace helpers. Apache-2.0. |
| NPVAR | https://github.com/MingGao97/NPVAR | Algorithm plus teaser pruning helper and internal ordering helper. Apache2.0. |
| EqVarDAG | https://github.com/WY-Chen/EqVarDAG | Top-down ordering and helpers for constructing a DAG from an ordering and for variable selection. GPL-3.0; see `EqVarDAG/COPYING`. |
| RESIT | Peters' nonlineCodeANM archive cited in supplement K.4 | GAM fitter, fit/independence interface and HSIC helpers. Copyright and licensing information are identified in the source headers and upstream archive. |
| GLOBE | https://github.com/osman-mian/Globe | GLOBE modules with adjustments to rpy2 data conversion. Upstream README identifies MIT licensing. |
| TOPIC | https://github.com/srhmm/topic | Topic/GAM fitting and the required GP, Fourier and embedding helpers. MIT. |
| TAM | https://github.com/MingGao97/TAM | `compute_caus_order` in `shared/utils/graphs.py`. Apache2.0; see `TAM_LICENSE`. |
| CAM | https://cran.r-project.org/package=CAM | Installed R dependency. Local wrappers perform neighborhood selection, ordering and pruning. |
| GSGES | https://github.com/py-why/causal-learn | Installed dependency. Local code handles fitting with a generalized score and parameter tuning, and includes a wrapper for extending a CPDAG to a DAG. |

Each component retains its own licensing terms; the MIT/Apache notices do not apply to every component. No license covering the entire repository is assigned to the authors’ code or data.
