"""Read headerless numeric matrices stored with commas or whitespace."""
from pathlib import Path
import numpy as np

def load_matrix(path):
    path=Path(path)
    with path.open() as f:
        first=next((line.strip() for line in f if line.strip()),'')
    if not first:raise ValueError(f'Empty matrix: {path}')
    result=np.loadtxt(path,delimiter=',' if ',' in first else None,ndmin=2)
    if not np.isfinite(result).all():raise ValueError(f'Non-finite matrix: {path}')
    return result

def save_matrix(path,matrix):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    np.savetxt(path,matrix,delimiter=',')
