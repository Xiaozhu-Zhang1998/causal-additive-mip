"""Figure 1: for every requested function/sample size/trial:

1. shared.special_fits.fit_teaser runs spline and linear MIP.
2. shared/baselines/teaser.R runs CAM and NPVAR.
3. shared.aggregation.aggregate_teaser joins the per-trial result CSVs.
Use --dry-run to print the calls without fitting.
"""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p/'shared/settings.py').exists())
sys.path[:0] = [str(ROOT),str(ROOT/'third_party')]
from shared.cli import main

if len(sys.argv)==1: sys.argv.extend(['--all-settings','--all-trials'])

if __name__ == '__main__': main('teaser','all', experiment='Fig01')
