import csv,json
from pathlib import Path
import numpy as np
from scipy.stats import t
from analyze import design,ref,fit,dump_csv


def run(root,reps=1000):
    out=root/'results'
    cells=json.loads((out/'analysis_results.json').read_text())['cells']
    a=np.array([r['age'] for r in cells]);m=np.array([r['mileage'] for r in cells])
    age_rows=[];curve_rows=[];cov_rows=[]
    for p in ['CO','NMHC','NOx']:
        y=np.array([r[p+'_mean'] for r in cells]);s=np.array([r[p+'_sd'] for r in cells]);n=np.array([r['n'] for r in cells])
        for model in ['A','AM','CAT']:
            X=design(a,m,model);c=ref(model)
            for age in np.unique(a):
                keep=a!=age;b=fit(X[keep],y[keep],1/s[keep]**2)[0]
                for g in np.flatnonzero(~keep):
                    age_rows.append(dict(pollutant=p,model=model,held_age=float(age),cell=cells[g]['cell'],observed=float(y[g]),predicted=float(X[g]@b),squared_error=float((y[g]-X[g]@b)**2)))
            b,e,h,V,H,B=fit(X,y,n/s**2)
            cov_rows.append(dict(pollutant=p,model=model,reference=float(c@b),working_se=float(np.sqrt(c@V@c)),hc3_se=float(np.sqrt(c@H@c)),plugin_se=float(np.sqrt(c@B@c))))
        b=fit(design(a,m,'CAT'),y,1/s**2)[0]
        for cl in [1,2,3]:
            aa=np.linspace(.5,5,100);pred=design(aa,np.full(100,cl),'CAT')@b
            curve_rows.extend(dict(pollutant=p,mileage=cl,age=float(ag),fitted_mean=float(v)) for ag,v in zip(aa,pred))
    dump_csv(out/'age_level_folds.csv',age_rows);dump_csv(out/'covariance_comparison.csv',cov_rows);dump_csv(out/'fitted_curves.csv',curve_rows)
    counts=np.array([20,80,30,25,100,20,40,90,25,30,80,20,40,100,30,20,80,25.])
    sd=np.linspace(1,4,18);se=sd/np.sqrt(counts);base=design(a,m,'AM')@np.array([2,.3,.2,.05])
    scenarios=[('Ordinal mean',base,None,0,False),('Monotone category mean',base+.4*(m==2)+1.5*(m==3),None,0,False),('Middle class dip',base-1.5*(m==2),None,0,False),('Middle class rise',base+1.5*(m==2),None,0,False)]
    scenarios += [(f'g{g+1:02} plus {k} SD',base,g,k,False) for g in [0,8,16] for k in [3,6]]
    scenarios += [('Ordinal mean t5 errors',base,None,0,True)]
    rng=np.random.default_rng(20261030);rows=[];raw=[]
    for label,mu,g,k,heavy in scenarios:
        stats={(mode,model):[] for mode in ['estimated','known'] for model in ['AM','CAT']}
        for rep in range(reps):
            z=rng.standard_t(5,18)*np.sqrt(3/5) if heavy else rng.normal(size=18)
            y=mu+se*z
            if g is not None:y[g]+=k*se[g]
            estimated=sd**2*rng.chisquare(counts-1)/(counts-1)
            for mode,var in [('estimated',estimated),('known',sd**2)]:
                for model in ['AM','CAT']:
                    X=design(a,m,model);c=ref(model);b,e,h,V,H,B=fit(X,y,counts/var)
                    err=float(c@b-mu[16]);ses=[np.sqrt(c@cov@c) for cov in [B,V,H]]
                    crit=[1.959964,t.ppf(.975,18-X.shape[1]),t.ppf(.975,18-X.shape[1])]
                    covers=[int(abs(err)<=q*s) for q,s in zip(crit,ses)]
                    stats[mode,model].append((err,*covers))
                    raw.append(dict(scenario=label,replication=rep+1,dispersion=mode,model=model,error=err,plugin_se=float(ses[0]),working_se=float(ses[1]),hc3_se=float(ses[2]),plugin_cover=covers[0],working_cover=covers[1],hc3_cover=covers[2]))
        for (mode,model),values in stats.items():
            ar=np.array(values);err=ar[:,0];rmse=np.sqrt(np.mean(err**2))
            rows.append(dict(scenario=label,dispersion=mode,model=model,replications=len(err),failures=0,bias=float(err.mean()),RMSE=float(rmse),RMSE_MCSE=float(np.std(err**2,ddof=1)/(2*rmse*np.sqrt(len(err)))),plugin_coverage=float(100*ar[:,1].mean()),working_coverage=float(100*ar[:,2].mean()),hc3_coverage=float(100*ar[:,3].mean())))
        print('Extension:',label,flush=True)
    dump_csv(out/'extended_simulations.csv',rows);dump_csv(out/'extended_replications.csv',raw)
