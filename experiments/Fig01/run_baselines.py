"""Figure 1 CAM/NPVAR baseline launcher; invokes the shared R implementation."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p/'shared/settings.py').exists())
sys.path[:0] = [str(ROOT),str(ROOT/'third_party')]
from shared.cli import main


if __name__ == '__main__': main('teaser','r_baselines', experiment='Fig01')
