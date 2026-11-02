"""Create a self-contained reproducibility ZIP including provided data and results."""
from pathlib import Path
import argparse,zipfile

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();root=Path(__file__).resolve().parents[1];output=args.output.expanduser().resolve()
    if output.exists():raise FileExistsError(output)
    excluded={'.git','__pycache__','.venv','.DS_Store','.Rhistory','__MACOSX','.ipynb_checkpoints'}
    output.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(output,'x',compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(root.rglob('*')):
            rel=path.relative_to(root)
            if rel.parts[0]=='outputs' or excluded.intersection(rel.parts) or rel.parts[:2] in [('data','generated'),('results','rerun'),('results','reaggregated')]:continue
            if path==output or path.is_symlink() or not path.is_file():continue
            archive.write(path,arcname=str(rel))
    print(output)
