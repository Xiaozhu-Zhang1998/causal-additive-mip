"""Generate synthetic datasets from the experiment settings."""
import numpy as np
from shared.paths import ROOT,output_folder,actual_trial
from shared.utils.io import save_matrix
from shared.utils.graphs import edge2adj

def generate(spec,trial,output_root=None,overwrite=False):
    from shared.dgp.er import generate_er_dag
    from shared.dgp.sem import gen_data_from_adj,gen_data_from_adj_tdisn,function_pool
    from shared.dgp.fixed_graphs import GRAPHS
    l=actual_trial(spec,trial)
    dest=output_folder(spec,output_root or ROOT/'data/generated');dest.mkdir(parents=True,exist_ok=True)
    if (dest/f'X_{l}.csv').exists() and not overwrite:raise FileExistsError(dest/f'X_{l}.csv')
    family=spec['family'];p=spec['p'];n=spec['n'];rng=np.random.default_rng()
    if 'graph' in spec:
        g=GRAPHS[spec['graph']];edges=[(i-1,j-1) for i,j in g['edges']];a=edge2adj(edges,p)
    else:a,_=generate_er_dag(p,spec['edge_density'])
    if family=='vanilla':
        from shared.dgp.vanilla import gen_data
        sd=.5*rng.beta(1,2,size=p)+.5
        X=gen_data(edges,n,p,lambda x:(np.sin(x)+np.cos(x))/2,sd)
    else:
        sd=np.exp(spec['tau']*rng.normal(size=p));sd/=np.sqrt(np.mean(sd**2))
        if family=='student_t':
            X=gen_data_from_adj_tdisn(a,n,p,function_pool,sd=sd,df=spec['df'])
        else:X=gen_data_from_adj(a,n,p,function_pool,sd=sd)
    save_matrix(dest/f'X_{l}.csv',X);save_matrix(dest/f'dag_oracle_{l}.csv',a)
    print(dest/f'X_{l}.csv')
