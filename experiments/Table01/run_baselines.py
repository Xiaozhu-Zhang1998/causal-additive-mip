"""Fit one or more Python baselines on a supplied interaction dataset."""

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / 'third_party')]

from shared.fitting import fit_python_baseline

METHODS = ('globe', 'topic', 'notears', 'gsges')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--setting', type=int, choices=range(13), metavar='0..12', required=True)
    parser.add_argument('--tau', type=float, choices=(0.1, 0.2), required=True)
    parser.add_argument('--trial', type=int, choices=range(20), metavar='0..19', required=True)
    parser.add_argument('--method', choices=(*METHODS, 'all'), default='all')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    spec = {
        'family': 'interaction',
        'folder': f'interaction/tau_{args.tau:g}/setting_{args.setting:02d}',
        'trials': 20,
    }
    methods = METHODS if args.method == 'all' else (args.method,)
    for method in methods:
        fit_python_baseline(spec, args.trial, method, overwrite=args.overwrite)


if __name__ == '__main__':
    main()
