import csv,json,re
from collections import Counter
import numpy as np
from scipy.optimize import minimize
from analyze import design,ref,fit,huber,dump_csv


def read(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))


def verify(root,reps):
    s=json.loads((root/'results/analysis_results.json').read_text());cells=s['cells']
    a=np.array([r['age'] for r in cells]);m=np.array([r['mileage'] for r in cells]);checks=[]
    def check(name,condition):
        checks.append(dict(check=name,passed=bool(condition)))
        if not condition:raise AssertionError(name)
    check('18 unique grouped inputs',len(cells)==len({r['cell'] for r in cells})==18)
    check('count discrepancy documented',s['counts_sum']==1000 and len(s['count_mismatch_cells'])==17)
    for p in ['CO','NMHC','NOx']:
        y=np.array([r[p+'_mean'] for r in cells]);sd=np.array([r[p+'_sd'] for r in cells])
        for model in ['A','AM','CAT']:
            X=design(a,m,model);w=1/sd**2;b,e,h,V,H,B=fit(X,y,w)
            independent=np.linalg.solve(X.T@(w[:,None]*X),X.T@(w*y))
            check(p+' '+model+' normal equations',np.allclose(b,independent,atol=1e-8))
            check(p+' '+model+' scale invariance',np.allclose(b,fit(X,y,w*100)[0],atol=1e-8))
            check(p+' '+model+' leverage trace',abs(h.sum()-X.shape[1])<1e-9)
            check(p+' '+model+' symmetric covariance',all(np.allclose(C,C.T) for C in [V,H,B]))
            errors=[]
            for g in range(18):
                k=np.arange(18)!=g;deleted=y[g]-X[g]@fit(X[k],y[k],w[k])[0]
                errors.append(abs(deleted-e[g]/(1-h[g])))
            check(p+' '+model+' PRESS',max(errors)<1e-8)
        X=design(a,m,'AM');w=1/sd**2;bh,scale,_=huber(X,y,w);q=w/w.sum()
        def obj(b):
            u=(y-X@b)/scale;z=abs(u);return np.sum(q*np.where(z<=1.345,u*u/2,1.345*z-1.345**2/2))
        def grad(b):return -X.T@(q*np.clip((y-X@b)/scale,-1.345,1.345))/scale
        opt=minimize(obj,bh,jac=grad,method='BFGS',options={'gtol':1e-9})
        check(p+' independent Huber objective',abs(obj(opt.x)-obj(bh))<1e-8)
        check(p+' Huber gradient',max(abs(grad(bh)))<1e-6)
        baseline=float(ref('AM')@fit(X,y,w)[0])
        check(p+' sensitivity baseline',all(r['low']-1e-8<=baseline<=r['high']+1e-8 for r in s['sensitivity'][p]))
        check(p+' sensitivity feasibility',all(v['constraint_violation']<1e-7 for r in s['sensitivity'][p] for v in r.get('runs',[]) if v['success']))
    check('648 successful linear deletion fits',sum(r['failed_folds'] for r in s['fits'])==0 and len(read(root/'results/cell_diagnostics.csv'))==648)
    check('54 successful Huber deletion fits',all(r['failed_folds']==0 for r in s['huber'].values()))
    check('30000 original simulation fits at full run',sum(r['replications'] for r in s['simulations'])==30*reps and all(r['failures']==0 for r in s['simulations']))
    ext=read(root/'results/extended_replications.csv')
    check('44000 extension fits at full run',len(ext)==44*reps)
    check('finite extension errors and standard errors',all(all(np.isfinite(float(r[k])) for k in ['error','plugin_se','working_se','hc3_se']) for r in ext))
    assigned=read(root/'data/resampled_assignments.csv');counts=Counter(r['Cell'] for r in assigned)
    check('resampled frequencies differ in 17 cells',sum(counts[c['cell']]!=c['n'] for c in cells)==17)
    comparison=[]
    def compare(table,key,metric,reported,calculated):
        text=str(reported);v=float(text);dec=len(text.split('.')[1]) if '.' in text else 0;tol=.50001*10**(-dec)
        comparison.append(dict(table=table,key=key,metric=metric,reported=v,recalculated=float(calculated),absolute_difference=float(abs(v-calculated)),rounding_tolerance=tol,matches_reported_rounding=bool(abs(v-calculated)<=tol)))
    def tab(n):
        with (root/'data'/f'reference_table_{n:02}.csv').open(newline='') as f:return list(csv.reader(f))[1:]
    def f(p,model,scheme='inverse_variance'):return next(r for r in s['fits'] if r['pollutant']==p and r['model']==model and r['scheme']==scheme)
    for r,c in zip(tab(1),cells):
        check('Input '+r[0]+' cell identifier',r[0]==c['cell'])
        for j,k in enumerate(['age','mileage','n','CO_mean','CO_sd','NMHC_mean','NMHC_sd','NOx_mean','NOx_sd'],1):compare('1',r[0],k,r[j],c[k])
    for r in tab(3):
        rr=[f(p,r[0]) for p in ['CO','NMHC','NOx']]
        compare('3',r[0],'coefficients',r[2],rr[0]['coefficients'])
        low,high=r[3].split('–')
        compare('3',r[0],'R2_min',low,min(v['R2'] for v in rr))
        compare('3',r[0],'R2_max',high,max(v['R2'] for v in rr))
        for value,p,v in zip(r[4].split(' / '),['CO','NMHC','NOx'],rr):compare('3',r[0]+' '+p,'LOCO_MSE',value,v['LOCO_MSE'])
    for r in tab(2):
        actual=f(r[0],r[1])
        for j,k in enumerate(['reference','working_se','hc3_se','R2','LOCO_MSE'],2):compare('2',r[0]+' '+r[1],k,r[j],actual[k])
    schemes={'Equal':'equal','Inverse SD':'inverse_sd','Inverse dispersion':'inverse_variance','Count informed':'count_informed'}
    for r in tab(4):
        actual=f(r[0],'AM',schemes[r[1]])
        for j,k in enumerate(['reference','LOCO_MSE','N_eff'],2):compare('4',r[0]+' '+r[1],k,r[j],actual[k])
    for r in tab(5):
        actual=f(r[0],'AM');hu=s['huber'][r[0]]
        for j,k,val in [(1,'max_deletion',actual['max_deletion']),(3,'max_leverage',actual['max_leverage']),(4,'Huber_reference',hu['reference']),(5,'Huber_LOCO_MSE',hu['LOCO_MSE'])]:compare('5',r[0],k,r[j],val)
    for r in tab(6):
        for p,text in zip(['CO','NMHC','NOx'],r[1:]):
            actual=next(v for v in s['sensitivity'][p] if v['rho']==float(r[0]))
            for k,v in zip(['low','high'],text.split(' to ')):compare('6',p+' rho '+r[0],k,v,actual[k])
    scenarios=list(dict.fromkeys(r['scenario'] for r in s['simulations']))
    for row,sc in zip(tab(7),scenarios):
        for text,method in zip(row[1:],['equal_AM','dispersion_AM','oracle_AM','dispersion_CAT','Huber_AM']):
            actual=next(r for r in s['simulations'] if r['scenario']==sc and r['method']==method);numbers=re.findall(r'[0-9.]+',text)
            for k,v in zip(['RMSE','RMSE_MCSE'],numbers):compare('7',sc+' '+method,k,v,actual[k])
    folds=read(root/'results/age_level_folds.csv')
    for r in tab(8):
        for j,model in enumerate(['A','AM','CAT'],1):compare('8',r[0]+' '+model,'MSE',r[j],np.mean([float(v['squared_error']) for v in folds if v['pollutant']==r[0] and v['model']==model]))
        for j,age in [(4,.5),(5,5)]:compare('8',r[0]+' CAT '+str(age),'fold_MSE',r[j],np.mean([float(v['squared_error']) for v in folds if v['pollutant']==r[0] and v['model']=='CAT' and float(v['held_age'])==age]))
    for r in tab(9):
        actual=f(r[0],r[1],'count_informed')
        for j,k in enumerate(['reference','working_se','hc3_se','plugin_se'],2):compare('9',r[0]+' '+r[1],k,r[j],actual[k])
    extension=read(root/'results/extended_simulations.csv')
    for table,mode in [(10,'estimated'),(11,'known')]:
        for r in tab(table):
            actual=next(v for v in extension if v['scenario']==r[0] and v['model']==r[1] and v['dispersion']==mode)
            numbers=re.findall(r'[0-9.]+',r[2])
            for k,v in zip(['RMSE','RMSE_MCSE'],numbers):compare(str(table),r[0]+' '+r[1],k,v,float(actual[k]))
            for j,k in enumerate(['plugin_coverage','working_coverage','hc3_coverage'],3):compare(str(table),r[0]+' '+r[1],k,r[j],float(actual[k]))
    dump_csv(root/'results/manuscript_comparison.csv',comparison)
    report=dict(checks=checks,passed=sum(c['passed'] for c in checks),failed=sum(not c['passed'] for c in checks),manuscript_values_compared=len(comparison),manuscript_values_matching=sum(r['matches_reported_rounding'] for r in comparison),manuscript_values_differing=sum(not r['matches_reported_rounding'] for r in comparison))
    (root/'validation/verification.json').write_text(json.dumps(report,indent=2))
    print('Verification:',report['passed'],'checks passed;',report['manuscript_values_differing'],'manuscript rounding differences.',flush=True)
