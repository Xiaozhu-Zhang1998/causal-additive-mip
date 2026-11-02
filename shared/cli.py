"""Small command-line adapters; scientific routines are imported only when run."""
import argparse,csv,subprocess,sys
from pathlib import Path
from shared.settings import make_spec,all_specs,experiment_specs,METHODS,TEASER_FUNCTIONS
from shared.paths import ROOT,source_folder,output_folder,auxiliary_file,assert_output_separate,actual_trial

def parser(family,action):
    p=argparse.ArgumentParser(description=f'{family}: {action}. Saved inputs are read-only; new outputs go to results/rerun.')
    p.add_argument('--trial',type=int,default=0,help='Zero-based trial within this setting (Chamber 0..9).')
    p.add_argument('--n',type=int);p.add_argument('--p',type=int);p.add_argument('--tau',type=float)
    p.add_argument('--edge_prob',type=float);p.add_argument('--graph');p.add_argument('--df',type=int)
    p.add_argument('--data-root');p.add_argument('--results-root',help='Provided estimates for --from-saved aggregation; defaults to results/.');p.add_argument('--output-root');p.add_argument('--overwrite',action='store_true')
    method_choices=(['all','mip'] if action=='mip' else ['all']) if family=='teaser' else METHODS+['pns','all']
    p.add_argument('--method',choices=method_choices,default='all')
    p.add_argument('--n-jobs',type=int,default=3)
    p.add_argument('--all-settings',action='store_true');p.add_argument('--all-trials',action='store_true')
    p.add_argument('--dry-run',action='store_true',help='Print commands/paths without loading a solver or fitting.')
    p.add_argument('--from-saved',action='store_true',help='Aggregation only: read provided results/ graphs and write results/reaggregated.')
    if family=='teaser':p.add_argument('--function',choices=TEASER_FUNCTIONS,default='atan(x^2)')
    return p

def select_specs(family,a,experiment=None):
    """Resolve result-specific defaults and reject settings outside that result."""
    if family!='er' or experiment is None:
        if not a.all_settings:
            return [make_spec(family,**{k:getattr(a,k) for k in ['n','p','tau','edge_prob','graph','df']})]
        fields={'n':'n','p':'p','tau':'tau','edge_prob':'edge_density','graph':'graph','df':'df'}
        selected=[s for s in all_specs(family) if all(getattr(a,arg) is None or
                    s.get(field)==getattr(a,arg) for arg,field in fields.items())]
        if not selected:raise ValueError('Specified filters select no documented settings.')
        return selected
    candidates=experiment_specs(family,experiment)
    fields={'n':'n','p':'p','tau':'tau','edge_prob':'edge_density'}
    candidates=[s for s in candidates if all(getattr(a,arg) is None or
                getattr(a,arg)==s[field] for arg,field in fields.items())]
    if not candidates:raise ValueError(f'Specified setting is not part of {experiment}.')
    if a.graph is not None or a.df is not None:raise ValueError('ER experiments do not take --graph or --df.')
    if a.all_settings:return candidates
    # Prefer the shared n500/p15/tau.2/rho.3 representative where available.
    preferred={'n':500,'p':15,'tau':.2,'edge_density':.3}
    return [max(candidates,key=lambda s:sum(s[k]==v for k,v in preferred.items()))]

def run_r(spec,trial,mode,a):
    l=actual_trial(spec,trial);source=source_folder(spec,a.data_root);dest=output_folder(spec,a.output_root)
    sup=auxiliary_file(spec,f'dag_superstructure_{l}.csv',a.data_root,a.output_root)
    command=['Rscript','--vanilla',str(ROOT/'shared/baselines/run_r.R'),str(ROOT),mode,
             str(source/f'X_{l}.csv'),str(dest),str(l),str(sup),a.method]
    if a.dry_run:print(' '.join(command))
    else:
        protected=dest/f'dag_CAMpartial_{l}.csv' if mode=='bootstrap' else dest/f'dag_superstructure_{l}.csv'
        if protected.exists() and not a.overwrite:raise FileExistsError(f'{protected}; pass --overwrite for rerun files')
        subprocess.run(command,check=True)

def main(family,action,experiment=None):
    a=parser(family,action).parse_args()
    allowed={'r_baselines':{'all','pns','cam','npvar','eqvar','resit'},
             'python_baselines':{'all','notears','globe','topic','gsges'},
             'mip':{'all','mip'},'all':{'all'}}
    if action in allowed and a.method not in allowed[action]:
        raise ValueError(f'--method for {action} must be one of {sorted(allowed[action])}')
    if a.output_root:
        assert_output_separate(a.output_root,a.data_root)
    if family=='teaser':
        return teaser_main(a,action)
    specs=select_specs(family,a,experiment)
    if action=='aggregate':
        from shared.aggregation import aggregate
        return aggregate(specs,a)
    for spec in specs:
        trials=range(spec['trials']) if a.all_trials else [a.trial]
        for trial in trials:
            if action=='all':
                if family not in ['vanilla','earlystop']:
                    run_r(spec,trial,'pns' if family in ['gamma','basis_size'] else 'methods',a)
                    if spec['priors']:run_r(spec,trial,'bootstrap',a)
                actions=['python_baselines','mip'] if family not in ['gamma','basis_size','earlystop','vanilla'] else ['comparison' if family=='vanilla' else 'mip']
            else:actions=[action]
            for item in actions:
                if item in ['r_baselines','bootstrap']:
                    run_r(spec,trial,'bootstrap' if item=='bootstrap' else 'methods',a);continue
                if a.dry_run:
                    if item=='generate':
                        print(f'generate: {spec} output={output_folder(spec,a.output_root or ROOT/"data/generated")}')
                    else:print(f'{item}: input={source_folder(spec,a.data_root)}/X_{actual_trial(spec,trial)}.csv output={output_folder(spec,a.output_root)}')
                    continue
                if item=='mip':
                    from shared.fitting import fit_mip
                    fit_mip(spec,trial,a.data_root,a.output_root,a.n_jobs,a.overwrite)
                    if action=='all' and family=='radial':
                        fit_mip(dict(spec,family='er'),trial,a.data_root,a.output_root,a.n_jobs,a.overwrite)
                elif item=='python_baselines':
                    from shared.fitting import fit_python_baseline
                    methods=['notears','globe','topic','gsges'] if a.method=='all' else [a.method]
                    for method in methods:fit_python_baseline(spec,trial,method,a.data_root,a.output_root,a.overwrite)
                elif item=='generate':
                    from shared.generation import generate
                    generate(spec,trial,a.output_root,a.overwrite)
                elif item=='comparison':
                    from shared.special_fits import fit_vanilla
                    fit_vanilla(spec,trial,a.data_root,a.output_root,a.overwrite)
                else:raise ValueError(item)
    if action=='all' and not a.dry_run:
        from shared.aggregation import aggregate
        aggregate(specs,a)

def teaser_main(a,action):
    if action not in {'mip','r_baselines','all','generate','aggregate'}:
        raise ValueError(f'Unsupported Figure 1 action: {action}')
    if action=='r_baselines' and a.method!='all':
        raise ValueError('Figure 1 fits CAM and NPVAR together; omit --method or use --method all.')
    functions=TEASER_FUNCTIONS if a.all_settings else [a.function]
    ns=range(100,2001,100) if a.all_settings else [500 if a.n is None else a.n]
    trials=range(100) if a.all_trials else [a.trial]
    if any(n not in range(100,2001,100) for n in ns) or any(not 0<=trial<100 for trial in trials):
        raise ValueError('Figure 1 uses n=100,200,...,2000 and trials 0..99.')
    if action=='aggregate':
        from shared.aggregation import aggregate_teaser
        return aggregate_teaser(a)
    for function in functions:
        for n in ns:
            for trial in trials:
                if a.dry_run:
                    if action in ['mip','all']:
                        print(f"CALL shared/special_fits.py: fit_teaser(trial={trial}, n={n}, function={function!r})")
                        print('  -> shared/mip/basis.py: gen_basis_feat; shared/mip/cutting_plane.py: solve_sem_DAG_OA_withinitial1')
                        print('  -> spline and linear fits, each on unequal/equal-noise saved datasets')
                    if action in ['r_baselines','generate','all']:
                        mode='generate' if action=='generate' else 'fit'
                        print(f"RUN Rscript shared/baselines/teaser.R: mode={mode!r}, n={n}, trial={trial}, function={function!r}")
                        if mode=='fit':print('  -> CAM::CAM; NPVAR; prune (unequal/equal-noise saved datasets)')
                    continue
                if action in ['mip','all']:
                    from shared.special_fits import fit_teaser
                    fit_teaser(trial,n,function,a.data_root,a.output_root,a.overwrite)
                if action in ['r_baselines','generate','all']:
                    dest=Path(a.output_root or ROOT/('data/generated' if action=='generate' else 'results/rerun'))/'rawdata_teaser'/function
                    protected=[dest/f'{noise}_{n}_{trial}.csv' for noise in ['equal','unequal']] if action=='generate' else [dest/f'R_{n}_{trial}.csv']
                    if not a.overwrite:
                        for path in protected:
                            if path.exists():raise FileExistsError(f'{path}; pass --overwrite for rerun files')
                    subprocess.run(['Rscript','--vanilla',str(ROOT/'shared/baselines/teaser.R'),str(ROOT),
                        'generate' if action=='generate' else 'fit',str(source_folder(dict(family='teaser',folder='rawdata_teaser/'+function),a.data_root)),
                        str(dest),str(n),str(trial),function],check=True)
    if action=='all':
        if a.dry_run:print('CALL shared/aggregation.py: aggregate_teaser -> join per-trial Python/R CSV files into RS.csv (preview only)')
        else:
            from shared.aggregation import aggregate_teaser
            aggregate_teaser(a)
