"""Fit the Python baseline methods for Figure 7."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p/'shared/settings.py').exists())
sys.path[:0] = [str(ROOT),str(ROOT/'third_party')]
from shared.cli import main


if __name__ == '__main__': main('radial','python_baselines', experiment='Fig07')
