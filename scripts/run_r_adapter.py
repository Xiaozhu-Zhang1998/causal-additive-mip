"""Resolve R trial/settings arguments without importing numerical Python packages."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from shared.cli import main
family,action,experiment=sys.argv[1:4]
sys.argv=sys.argv[:1]+sys.argv[4:]
main(family,action,experiment=experiment)
