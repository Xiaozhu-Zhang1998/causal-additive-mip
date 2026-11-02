"""Run fitting and aggregation for Table 7."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p/'shared/settings.py').exists())
sys.path[:0] = [str(ROOT),str(ROOT/'third_party')]
from shared.cli import main

if len(sys.argv)==1: sys.argv.extend(['--all-settings','--all-trials'])

if __name__ == '__main__': main('fixed','all', experiment='Table07')
