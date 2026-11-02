"""Generate 30 bootstrap datasets, each with 118 rows sampled with replacement."""
from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd

if __name__=='__main__':
    root=next(p for p in Path(__file__).resolve().parents if (p/'shared/settings.py').exists())
    sys.path.insert(0,str(root))
    from shared.paths import assert_output_separate
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=root/'data/real/microarray/Arabidopsis_MVA_standardized.csv')
    parser.add_argument('--output-root',type=Path,default=root/'data/generated')
    args=parser.parse_args()
    assert_output_separate(args.output_root)
    X=pd.read_csv(args.input,header=None)
    if X.shape!=(118,15) or not np.isfinite(X.to_numpy()).all():
        raise ValueError('Expected a finite 118-by-15 standardized matrix in the documented gene order.')
    B=30
    dest=args.output_root/'Expression/subsamples'
    targets=[dest/f'X_{l}.csv' for l in range(B)]
    for target in targets:
        if target.exists():raise FileExistsError(target)
    dest.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng()
    for target in targets:
        idx=rng.choice(len(X),size=len(X),replace=True)
        X.iloc[idx].to_csv(target,index=False,header=False)
    print(f'Generated {B} bootstrap samples of size {len(X)} in {dest}.')
