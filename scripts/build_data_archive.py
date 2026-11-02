"""Optional archive of provided data/real and data/simulated inputs; not required for Git setup."""
import argparse,csv,hashlib,tarfile
from pathlib import Path

def data_files(source):
    """Discover actual input files in the self-contained data layout."""
    for name in ['real','simulated']:
        for path in sorted((source/name).rglob('*')):
            if path.is_file() and not path.is_symlink() and path.name not in ['.DS_Store']:
                yield path

def build_archive(source,output):
    source=Path(source).resolve();output=Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    if output==source or source in output.parents:raise ValueError('Write the archive outside the input data folder')
    hashes=output.with_suffix(output.suffix+'.sha256.csv')
    if hashes.exists():raise FileExistsError(hashes)
    output.parent.mkdir(parents=True,exist_ok=True)
    count=0
    with tarfile.open(output,'w:gz') as archive,hashes.open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['path','sha256'])
        for path in data_files(source):
            rel=Path('data')/path.relative_to(source)
            archive.add(path,arcname=str(rel),recursive=False)
            h=hashlib.sha256()
            with path.open('rb') as data:
                for block in iter(lambda:data.read(1024*1024),b''):h.update(block)
            writer.writerow([str(rel),h.hexdigest()]);count+=1
    return count

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=Path(__file__).resolve().parents[1]/'data')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    count=build_archive(args.source,args.output)
    print(args.output,f'({count} files)')
