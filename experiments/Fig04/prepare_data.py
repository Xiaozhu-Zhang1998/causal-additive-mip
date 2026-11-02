"""Read the original Chamber data and recreate the 70 standardized block-position subsamples."""
from pathlib import Path
import sys,argparse
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'shared/settings.py').exists())
sys.path.insert(0,str(ROOT))
from shared.paths import data_root,assert_output_separate
from shared.utils.io import load_matrix,save_matrix
import numpy as np

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',help='Alternative data root (curated or legacy layout)')
    parser.add_argument('--input',type=Path,help='Original full Chamber matrix; defaults to data/real/causal_chamber/X.csv')
    parser.add_argument('--output-root',default=str(ROOT/'data/generated'))
    args=parser.parse_args()
    base=data_root(args.data_root)
    raw=args.input or base/'real/causal_chamber/X.csv'
    if args.input is None and not raw.exists():raw=base/'causal_chamber/X.csv'
    X=load_matrix(raw)
    assert_output_separate(args.output_root,args.data_root)
    if X.shape!=(10000,8):raise ValueError('Expected the saved 10000 x 8 Chamber X.csv')
    X=(X-X.mean(axis=0))/X.std(axis=0)
    truth=np.array([[0,0,0,1,1,1,1,0],[0,0,0,1,1,1,1,0],[0,0,0,1,1,1,1,0],
                    [0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0],[0,0,0,0,0,0,1,0]])
    dest=Path(args.output_root)/'causal_chamber'
    for ni,n in enumerate(range(200,801,100)):
        for trial in range(10):
            l=ni*10+trial
            if (dest/f'X_{l}.csv').exists():raise FileExistsError(dest/f'X_{l}.csv')
            save_matrix(dest/f'X_{l}.csv',X[np.arange(trial,10000,10)[:n]])
            save_matrix(dest/f'dag_oracle_{l}.csv',truth)
