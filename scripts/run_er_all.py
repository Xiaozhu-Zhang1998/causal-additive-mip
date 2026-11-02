"""Run all32 ER settings, using the same per-trial functions as the figure/table runners."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'third_party')]
from shared.cli import main
if len(sys.argv)==1:sys.argv.extend(['--all-settings','--all-trials'])
if __name__=='__main__':main('er','all')
