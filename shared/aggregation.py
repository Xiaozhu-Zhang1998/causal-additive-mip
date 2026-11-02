"""Describe graph/time files for the shared R SID/SHD/F1 aggregator."""
import csv,json,subprocess
from pathlib import Path
from shared.paths import ROOT,source_folder,provided_result_folder,output_folder,actual_trial
from shared.settings import METHODS,ZETAS,EARLY_GAPS

def jobs(spec,args):
    family=spec['family'];src=source_folder(spec,args.data_root)
    fitted=provided_result_folder(spec,getattr(args,'results_root',None)) if args.from_saved else output_folder(spec,args.output_root)
    result=[]
    methods=['mip'] if family in ['gamma','basis_size','earlystop'] else METHODS+['mip_rad'] if family=='radial' else METHODS
    if args.method!='all':methods=[args.method]
    for trial in range(spec['trials']):
        l=actual_trial(spec,trial)
        variants=range(3,12) if family=='gamma' else ZETAS if family=='basis_size' else EARLY_GAPS if family=='earlystop' else [None]
        for value in variants:
            for method in methods:
                stem='mipRad' if method=='mip_rad' else method
                suffix='' if value is None else f'_{value}'
                # Accept both 0 and 0.0 in numeric filename tags.
                graphfile=fitted/f'dag_{stem}_{l}{suffix}.csv'
                if not graphfile.exists() and value==0:
                    alternate=fitted/f'dag_{stem}_{l}_0.0.csv'
                    if alternate.exists():suffix='_0.0';graphfile=alternate
                timefile=fitted/f'time_{stem}_{l}{suffix}.txt'
                count=8 if family=='gamma' else 1 if family in ['basis_size','earlystop','microarray'] else (12 if spec['priors'] else 6) if method in ['mip','mip_rad'] else 3 if method in ['notears','gsges'] else 1
                if family=='microarray':count=3 if method in ['notears','gsges'] else 1
                meta=fitted/f'metadata_{stem}_{l}{suffix}.json'
                if meta.exists():count=json.loads(meta.read_text())['tuning_configurations']
                row={k:spec.get(k,'') for k in ['n','p','tau','edge_density','graph','df']}
                row.update(family=family,method=method,l=l,trial=l,
                    estimate=str(graphfile),truth=str(src/f'dag_oracle_{l}.csv'),runtime=str(timefile),
                    tuning_configurations=count,zeta=value if family=='basis_size' else '',
                    gamma=value if family=='gamma' else '')
                if family=='earlystop':row['tau']=value
                result.append(row)
        if args.method=='all' and family in ['er','fixed','radial']:
            for method,count in [('pns',1),('bootstrap',20)]:
                row={k:spec.get(k,'') for k in ['n','p','tau','edge_density','graph','df']}
                row.update(family=family,method=method,l=l,trial=l,estimate='',truth='',
                           runtime=str(fitted/f'time_{method}_{l}.txt'),tuning_configurations=count,zeta='',gamma='')
                result.append(row)
    return result

def aggregate(specs,args):
    for spec in specs:
        dest=output_folder(spec,args.output_root or (ROOT/'results/reaggregated' if args.from_saved else None))
        if spec['family']=='vanilla':
            if args.from_saved:raise ValueError('Table2 historical RS.csv is already provided; aggregate.py combines new comparison_*.csv batches.')
            if args.dry_run:
                print(f'{spec["folder"]}: join comparison_*.csv -> {dest}/RS.csv')
                continue
            rows=[]
            for f in sorted(dest.glob('comparison_*.csv')):
                with f.open() as stream:rows.extend(csv.DictReader(stream))
            if not rows:print('No completed comparison files:',dest);continue
            with (dest/'RS.csv').open('w',newline='') as stream:
                w=csv.DictWriter(stream,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
            continue
        rows=jobs(spec,args)
        if args.dry_run:print(f'{spec["folder"]}: {len(rows)} candidate graph files -> {dest}/RS.csv');continue
        dest.mkdir(parents=True,exist_ok=True)
        jobfile=dest/'aggregation_inputs.csv'
        with jobfile.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        subprocess.run(['Rscript','--vanilla',str(ROOT/'shared/aggregate.R'),str(ROOT),str(jobfile),str(dest)],check=True)

def aggregate_teaser(args):
    if args.from_saved:raise ValueError('Historical teaser RS.csv files are already provided; this command joins new Python/R per-trial batches.')
    base=Path(args.output_root or ROOT/'results/rerun')/'rawdata_teaser'
    if args.dry_run:
        print(f'CALL aggregate_teaser: join python_*.csv and R_*.csv under {base} into each RS.csv (preview only)')
        return
    for folder in sorted(base.glob('*')):
        if not folder.is_dir():continue
        rows=[]
        for pattern in ['python_*.csv','R_*.csv']:
            for f in sorted(folder.glob(pattern)):
                with f.open() as stream:rows.extend(csv.DictReader(stream))
        if rows:
            with (folder/'RS.csv').open('w',newline='') as stream:
                w=csv.DictWriter(stream,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
            print(folder/'RS.csv',len(rows),'rows')
