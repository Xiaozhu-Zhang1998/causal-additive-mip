"""Introductory example and comparison of the two MIP formulations."""
import csv,time,json
from pathlib import Path
import numpy as np
from shared.paths import output_folder,source_folder,actual_trial,ROOT
from shared.utils.io import load_matrix,save_matrix
from shared.utils.graphs import edge2adj,edge_diff

def fit_teaser(trial,n,function,input_root=None,output_root=None,overwrite=False):
    from shared.mip.basis import gen_basis_feat
    from shared.mip.cutting_plane import solve_sem_DAG_OA_withinitial1
    if not 0<=trial<100 or n not in range(100,2001,100):raise ValueError('Teaser: trials 0..99; n=100,200,...,2000')
    folder=source_folder(dict(family='teaser',folder='rawdata_teaser/'+function),input_root)
    dest=Path(output_root or ROOT/'results/rerun')/'rawdata_teaser'/function
    dest.mkdir(parents=True,exist_ok=True);target=dest/f'python_{n}_{trial}.csv'
    if target.exists() and not overwrite:raise FileExistsError(target)
    truth=edge2adj([(0,1),(0,2),(1,2)],5);rows=[]
    for var in ['unequal','equal']:
        X=load_matrix(folder/f'{var}_{n}_{trial}.csv')
        for kind,method in [('spline','Mips'),('linear','linear')]:
            setting={'type':kind,**({'n_knots':3,'degree':3,'include_intercept':False} if kind=='spline' else {})}
            basis,_=gen_basis_feat(X,setting);start=time.perf_counter()
            g,_,_=solve_sem_DAG_OA_withinitial1(X,basis,[(i,j) for i in range(5) for j in range(5)],
                gamma=10,lam2=.01,mip_gap=0,n_jobs=3,time_limit=None,verbose=False)
            a=edge2adj([e for e,v in g.items() if v==1],5);elapsed=time.perf_counter()-start
            save_matrix(dest/f'dag_{method}_{var}_{n}_{trial}.csv',a)
            rows.append(dict(n=n,m=trial,var=var,method=method,rec=bool(np.all(a==truth)),
                             SHD=float(np.abs(a-truth).sum()),func=function,time=elapsed))
    with target.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(target)

def fit_vanilla(spec,trial,input_root=None,output_root=None,overwrite=False):
    from shared.dgp.fixed_graphs import GRAPHS
    from shared.mip.basis import gen_basis_feat,gen_suff_stat
    from shared.mip.cutting_plane import solve_sem_DAG_OA_withinitial1
    from shared.mip.vanilla import Mips_DAG_diffvar
    l=actual_trial(spec,trial);dest=output_folder(spec,output_root);dest.mkdir(parents=True,exist_ok=True)
    target=dest/f'comparison_{l}.csv'
    if target.exists() and not overwrite:raise FileExistsError(target)
    X=load_matrix(source_folder(spec,input_root)/f'X_{l}.csv');g=GRAPHS[spec['graph']];p=g['p']
    truth=[(i-1,j-1) for i,j in g['edges']]
    sup=[(i-1,j-1) for i,j in g['moral']]+[(j-1,i-1) for i,j in g['moral']]
    setting=dict(type='spline',n_knots=3,degree=3,include_intercept=False);rows=[]
    start=time.perf_counter();basis,_=gen_basis_feat(X,setting)
    a,_,obj=solve_sem_DAG_OA_withinitial1(X,basis,sup,gamma=1,lam2=.005,mip_gap=0,n_jobs=3,verbose=False)
    edges=[e for e,v in a.items() if v==1]
    # Restore the constant omitted by the projected objective.
    rows.append(dict(l=l,method='MIP_OA',objval=float(obj+p*(1+np.log(2))),
                     time=time.perf_counter()-start,shd=float(edge_diff(truth,edges,p))))
    save_matrix(dest/f'dag_MIP_OA_{l}.csv',edge2adj(edges,p))
    start=time.perf_counter();suff=gen_suff_stat(X,setting)
    edges,obj=Mips_DAG_diffvar(suff['Sigma'],p,.005,1,suff['R'],sup,tau=0,verbose=0,time_limit=300)
    rows.append(dict(l=l,method='MIP_vanilla',objval=float(obj),time=time.perf_counter()-start,
                     shd=float(edge_diff(truth,edges,p))))
    save_matrix(dest/f'dag_MIP_vanilla_{l}.csv',edge2adj(edges,p))
    with target.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(target)

