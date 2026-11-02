"""Experiment settings and dataset collections. Trial indices start at zero.
"""
from itertools import product

ER_SETTINGS = [(10,.3),(10,.4),(10,.5),(15,.15),(15,.3),(20,.15),(20,.2),(25,.2)]
CONDITIONS = {(500,.1):'c1', (300,.1):'c2', (500,.2):'c3', (300,.2):'c4'}
GRAPH_SIZES = dict(dsep=6, asia=8, bowling=9, rain=14, inssmall=15,
                   cloud=16, funnel=18, galaxy=20, insurance=27, factors=27)
MAIN_LAMBDAS = [.015,.02,.025,.03,.035,.04]  # Penalty grid for the main simulations.
METHODS = ['gsges','resit','notears','globe','topic','eqvar','npvar','cam','mip']
ZETAS = [0,.13,.2,.27,.33,.4,.47,.53]
EARLY_GAPS = [0,.01,.02,.03,.04,.05,.06,.07,.08,.09,.1]
TEASER_FUNCTIONS = ['atan(x^2)','exp(x)','0.5x^3','abs(x)']


def make_spec(family, n=None, p=None, tau=None, edge_prob=None, graph=None, df=None):
    """Identify a saved dataset collection and its experiment-specific prior policy."""
    if family not in ['er','radial']:
        defaults={'fixed':(500,GRAPH_SIZES.get(graph or 'inssmall'),.1,None),
          'vanilla':(500,GRAPH_SIZES.get(graph or 'insurance'),.1,None),
          'student_t':(500,10,.2,.3),'chamber':(n or 800,8,None,None),
          'microarray':(118,15,None,None),'gamma':(300,10,.2,.3),
          'basis_size':(500,15,.2,.3),'earlystop':(500,15,.1,.3)}
        if family in defaults:
            for name,given,expected in zip(['n','p','tau','edge_prob'],[n,p,tau,edge_prob],defaults[family]):
                if given is not None and given!=expected:
                    raise ValueError(f'{family} fixes {name}={expected}; this CLI does not silently change the design.')
    if family in ['er','radial']:
        n = n if n is not None else 500
        tau = tau if tau is not None else (.1 if family=='radial' else .2)
        p = p if p is not None else 15
        rho = edge_prob if edge_prob is not None else .3
        if (n,tau) not in CONDITIONS or (p,rho) not in ER_SETTINGS:
            raise ValueError('Choose one of the documented ER sample/noise/graph settings.')
        ck = CONDITIONS[n,tau]
        if family=='radial':
            if (n,tau)!=(500,.1) or (p,rho) not in [(10,.3),(10,.4),(15,.15),(15,.3)]:
                raise ValueError('Figure 7 uses n=500,tau=.1 and its four documented graph settings.')
            folder=f'rawdata_newbasis/rawdata_p_{p}_{round(rho*100):03d}'
        else:
            folder=f'rawdata_{ck}/rawdata_p_{p}_{round(rho*100):03d}_{ck}'
        return dict(family=family, folder=folder, n=n,p=p,tau=tau,edge_density=rho,
                    trials=50, priors=not(p==10 and rho in [.3,.4]))
    if family in ['fixed','vanilla']:
        graph=graph or ('inssmall' if family=='fixed' else 'insurance')
        if graph not in GRAPH_SIZES or (family=='vanilla' and graph not in ['dsep','asia','bowling','inssmall','insurance']):
            raise ValueError('Graph is not part of this experiment.')
        folder=f'rawdata_graphs/rawdata_graph_{graph}' if family=='fixed' else f'rawdata_vanillacomp/{graph}'
        return dict(family=family,folder=folder,n=500,p=GRAPH_SIZES[graph],graph=graph,
                    tau=.1,trials=50,priors=family=='fixed')
    if family=='student_t':
        df=25 if df is None else df
        if df not in [5,10,15,20,25,30]:raise ValueError('df must be 5,10,15,20,25,30.')
        return dict(family=family,folder=f'rawdata_tdisn/rawdata_p_10_{df:02d}',
                    n=500,p=10,tau=.2,edge_density=.3,df=df,trials=50,priors=False)
    if family=='chamber':
        n=800 if n is None else n
        if n not in range(200,801,100):raise ValueError('Chamber n must be 200,300,...,800.')
        return dict(family=family,folder='causal_chamber',n=n,p=8,trials=10,
                    trial_offset=(n//100-2)*10,priors=False)
    if family=='microarray':
        return dict(family=family,folder='Expression/subsamples',n=118,p=15,trials=30,priors=False)
    if family in ['gamma','basis_size','earlystop']:
        return dict(family=family,folder={'gamma':'rawdata_gamma','basis_size':'rawdata_rn',
                    'earlystop':'rawdata_earlystop'}[family],n=300 if family=='gamma' else 500,
                    p=10 if family=='gamma' else 15,tau=.1 if family=='earlystop' else .2,
                    edge_density=.3,trials=50,priors=False,
                    **({'graph':'inssmall'} if family=='earlystop' else {}))
    raise ValueError(f'Unknown family {family}')


def all_specs(family):
    if family=='er':
        return [make_spec(family,n=n,tau=t,p=p,edge_prob=r) for (n,t),(p,r) in product(CONDITIONS,ER_SETTINGS)]
    if family=='radial':return [make_spec(family,p=p,edge_prob=r) for p,r in [(10,.3),(10,.4),(15,.15),(15,.3)]]
    if family in ['fixed','vanilla']:
        names=GRAPH_SIZES if family=='fixed' else ['dsep','asia','bowling','inssmall','insurance']
        return [make_spec(family,graph=g) for g in names]
    if family=='student_t':return [make_spec(family,df=d) for d in [5,10,15,20,25,30]]
    if family=='chamber':return [make_spec(family,n=n) for n in range(200,801,100)]
    return [make_spec(family)]


def experiment_specs(family, experiment=None):
    """Settings displayed by one result; other families keep their documented grid."""
    if family!='er' or experiment is None:
        return all_specs(family)
    if experiment=='Fig02':
        return [make_spec('er',n=n,p=15,tau=t,edge_prob=.3) for n,t in CONDITIONS]
    if experiment=='Fig03':
        return [make_spec('er',n=500,p=p,tau=.2,edge_prob=r)
                for p,r in [(10,.5),(15,.3),(20,.2),(25,.2)]]
    conditions={'Table03':(500,.1),'Table04':(300,.1),
                'Table05':(500,.2),'Table06':(300,.2)}
    if experiment in conditions:
        n,t=conditions[experiment]
        return [make_spec('er',n=n,p=p,tau=t,edge_prob=r) for p,r in ER_SETTINGS]
    raise ValueError(f'Unknown ER experiment: {experiment}')
