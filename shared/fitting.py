"""One-dataset fits shared by terminal and sequential experiment runners."""
import json,time
from pathlib import Path
import numpy as np
from shared.settings import MAIN_LAMBDAS, ZETAS, EARLY_GAPS
from shared.paths import ROOT,source_folder,output_folder,auxiliary_file,actual_trial
from shared.utils.io import load_matrix,save_matrix

def basis_setting(kind='spline'):
    return {'type':kind,**({'n_knots':3,'degree':3,'include_intercept':False} if kind=='spline'
                         else {'n_centers':6,'bandwidth':.5} if kind=='rbf' else {})}

def fit_mip(spec,trial,input_root=None,output_root=None,n_jobs=3,overwrite=False):
    from shared.mip.basis import gen_basis_feat,loss_from_edges
    from shared.mip.cutting_plane import solve_sem_DAG_OA_withinitial1
    from shared.mip.priors import prior_edges
    from shared.utils.graphs import edge2adj
    l=actual_trial(spec,trial);dest=output_folder(spec,output_root);dest.mkdir(parents=True,exist_ok=True)
    X=load_matrix(source_folder(spec,input_root)/f'X_{l}.csv');n,p=X.shape
    if (n,p)!=(spec['n'],spec['p']):raise ValueError(f'Dataset shape {X.shape} disagrees with setting {(spec["n"],spec["p"])}')
    family=spec['family']
    if family=='earlystop':
        from shared.dgp.fixed_graphs import GRAPHS
        moral=GRAPHS['inssmall']['moral'];edges=[(i-1,j-1) for i,j in moral]+[(j-1,i-1) for i,j in moral]
    else:
        sup=load_matrix(auxiliary_file(spec,f'dag_superstructure_{l}.csv',input_root,output_root))
        edges=[tuple(map(int,x)) for x in np.argwhere(sup)]
    priors=[(None,None)]
    if spec['priors']:
        priors=[prior_edges(load_matrix(auxiliary_file(spec,f'dag_{name}partial_{l}.csv',input_root,output_root)),edges,p) for name in ['CAM','NPVAR']]
    configurations=[(5,None,None)]
    if family=='gamma':configurations=[(g,None,None) for g in range(3,12)]
    elif family=='basis_size':configurations=[(5,z,None) for z in ZETAS]
    elif family=='earlystop':configurations=[(5,None,t) for t in EARLY_GAPS]
    for gamma,zeta,gap in configurations:
        tag=f'_{gamma}' if family=='gamma' else f'_{zeta}' if zeta is not None else f'_{gap}' if gap is not None else ''
        method='mipRad' if family=='radial' else 'mip'
        dagfile=dest/f'dag_{method}_{l}{tag}.csv'
        if dagfile.exists() and not overwrite:raise FileExistsError(f'{dagfile}; pass --overwrite for rerun outputs only')
        start=time.perf_counter()
        setting=basis_setting('rbf' if family=='radial' else 'linear' if zeta==0 else 'spline')
        if zeta is not None and zeta!=0:setting['n_knots']=int(n**zeta)
        basis,R=gen_basis_feat(X,setting)
        lambdas=[.01] if family in ['basis_size','earlystop'] else [.08] if family=='microarray' else [.01,.015,.02,.025,.03,.035,.04,.045] if family=='gamma' else MAIN_LAMBDAS
        candidates=[]
        for order,stable in priors:
            for lam in lambdas:
                g,_,obj=solve_sem_DAG_OA_withinitial1(X,basis,edges,gamma=gamma,lam2=lam,
                    mip_gap=gap or 0,partial_order=order,known_edges=stable,n_jobs=n_jobs,
                    time_limit=500 if family=='earlystop' else 120 if family=='microarray' else 180,verbose=False)
                selected=[e for e,v in g.items() if v==1]
                bic=float(loss_from_edges(selected,X,setting)+len(selected)*R*np.log(n)/n)
                candidates.append((bic,selected,lam,float(obj)))
        best=min(candidates,key=lambda c:c[0]);elapsed=time.perf_counter()-start
        save_matrix(dagfile,edge2adj(best[1],p))
        (dest/f'time_{method}_{l}{tag}.txt').write_text(str(elapsed)+'\n')
        metadata=dict(trial=l,method=method,n=n,p=p,basis=setting,inverse_gamma=gamma,
                      lambda_grid=lambdas,selected_lambda=best[2],bic=best[0],objective=best[3],
                      tuning_configurations=len(candidates),time_total=elapsed,
                      time_per_configuration=elapsed/len(candidates),relative_gap=gap or 0)
        (dest/f'metadata_{method}_{l}{tag}.json').write_text(json.dumps(metadata,indent=2)+'\n')
        print(dagfile)

def fit_python_baseline(spec,trial,method,input_root=None,output_root=None,overwrite=False):
    from shared.mip.basis import gen_basis_feat,loss_from_adj
    l=actual_trial(spec,trial);dest=output_folder(spec,output_root);dest.mkdir(parents=True,exist_ok=True)
    target=dest/f'dag_{method}_{l}.csv'
    if target.exists() and not overwrite:raise FileExistsError(f'{target}; pass --overwrite')
    X=load_matrix(source_folder(spec,input_root)/f'X_{l}.csv');n,p=X.shape
    start=time.perf_counter();count=1
    if method=='notears':
        import torch
        from notears.nonlinear import NotearsMLP,notears_nonlinear
        from shared.baselines.notears import notears_bic
        torch.set_default_dtype(torch.double);candidates=[]
        for lam in [.001,.01,.1]:
            model=NotearsMLP(dims=[p,5,1],bias=True)
            w=notears_nonlinear(model,X,lambda1=lam,lambda2=.01)
            candidates.append((notears_bic(model,X,w)['bic']/n,w!=0))
        adj=min(candidates,key=lambda x:x[0])[1];count=3
    elif method=='gsges':
        from shared.baselines.gsges import fit_gsges
        setting=basis_setting();_,R=gen_basis_feat(X,setting);candidates=[]
        for lam in [.001,.01,.1]:
            a=fit_gsges(X,lam=lam,maxP=5,kfold=5)
            candidates.append((float(loss_from_adj(a,X,setting)+np.sum(a)*R*np.log(n)/n),a))
        adj=min(candidates,key=lambda x:x[0])[1];count=3
    elif method=='globe':
        from shared.baselines.globe import fit_globe
        adj=fit_globe(X,max_interactions=2,log_results=False,verbose=False)
    elif method=='topic':
        from shared.baselines.topic import fit_topic
        adj=fit_topic(X,extra_refinement=True,verbose=0)['adjacency']
    else:raise ValueError('Python methods: notears,globe,topic,gsges')
    elapsed=time.perf_counter()-start;save_matrix(target,adj)
    (dest/f'time_{method}_{l}.txt').write_text(str(elapsed)+'\n')
    (dest/f'metadata_{method}_{l}.json').write_text(json.dumps(dict(tuning_configurations=count,time_total=elapsed,time_per_configuration=elapsed/count),indent=2)+'\n')
    print(target)
