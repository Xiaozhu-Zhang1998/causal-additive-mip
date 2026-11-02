"""Resolve preserved inputs, provided estimates, and separate rerun outputs."""
from pathlib import Path
import os

ROOT = Path(__file__).resolve().parents[1]

def data_root(value=None):
    """Input root: repository data/ by default, or an explicit alternative tree."""
    return Path(value or os.environ.get('JASA_DATA_ROOT', ROOT/'data')).expanduser().resolve()

def source_folder(spec, root=None):
    """Locate observations and true graphs in provided or generated data folders."""
    base=data_root(root)
    if spec['family']=='chamber':folder=base/'real/causal_chamber/subsamples'
    elif spec['family']=='microarray':folder=base/'real/microarray/subsamples'
    else:folder=base/'simulated'/spec['folder']
    if folder.is_dir():return folder
    legacy=base/spec['folder']
    if legacy.is_dir():return legacy
    if spec['family']=='radial':
        from shared.settings import make_spec
        return source_folder(make_spec('er',n=500,p=spec['p'],tau=.1,edge_prob=spec['edge_density']),root)
    return folder

def provided_result_folder(spec, root=None):
    """Original estimates, priors, timings and summaries, copied without recomputation."""
    return Path(root or ROOT/'results').expanduser().resolve()/spec['folder']

def output_folder(spec, root=None):
    return Path(root or ROOT/'results/rerun')/spec['folder']

def auxiliary_file(spec, name, input_root=None, output_root=None):
    """Prefer newly computed PNS/prior files, then provided intermediate estimates."""
    fitted=output_folder(spec,output_root)/name
    if fitted.exists():return fitted
    base=data_root(input_root)
    if (base/'real').is_dir() or (base/'simulated').is_dir():
        saved=provided_result_folder(spec,base.parent/'results')/name
        if saved.exists():return saved
    # Alternative input folders can include their own intermediate matrices.
    return source_folder(spec,input_root)/name

def assert_output_separate(destination, input_root=None):
    """Protect supplied datasets while permitting data/generated/ as a new-output tree."""
    base=data_root(input_root);dest=Path(destination).expanduser().resolve()
    protected=[(base/'real').resolve(),(base/'simulated').resolve()] if (base/'real').exists() or (base/'simulated').exists() else [base]
    if dest==base or any(dest==p or p in dest.parents or dest in p.parents for p in protected):
        raise ValueError('Output must be outside preserved input-data trees.')

def actual_trial(spec, trial):
    if not 0 <= trial < spec['trials']:
        raise ValueError(f"--trial must be between 0 and {spec['trials']-1} for this setting")
    return trial+spec.get('trial_offset',0)
