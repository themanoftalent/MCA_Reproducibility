import sys, json, csv, hashlib, argparse, platform
from pathlib import Path
import numpy as np
import openpyxl
from scipy.optimize import minimize, linprog

SEED=20261017
POL=['CO','NMHC','NOx']

def design(a,m,kind):
 if kind=='A':return np.column_stack([np.ones(len(a)),a])
 if kind=='AM':return np.column_stack([np.ones(len(a)),a,m,a*m])
 if kind=='CAT':return np.column_stack([np.ones(len(a)),a,m==2,m==3,a*(m==2),a*(m==3)])
 raise ValueError(kind)
def ref(kind):return design(np.array([5.]),np.array([2.]),kind)[0]
def fit(X,y,w):
 Z=X*np.sqrt(w[:,None]);z=y*np.sqrt(w)
 b,_,rank,_=np.linalg.lstsq(Z,z,rcond=None)
 if rank<X.shape[1]:raise ValueError('rank deficient')
 resid=y-X@b;inv=np.linalg.inv(Z.T@Z);h=w*np.einsum('ij,jk,ik->i',X,inv,X)
 cov=(np.sum(w*resid**2)/(len(y)-X.shape[1]))*inv
 hc3=inv@(X.T@((w*resid/(1-h))[:,None]**2*X))@inv
 return b,resid,h,cov,hc3,inv

def huber(X,y,w):
 w=w/w.sum();b,res,_,_,_,_=fit(X,y,w)
 scale=1.4826*np.median(abs(res-np.median(res)))
 if not np.isfinite(scale) or scale<1e-12:raise ValueError('invalid Huber scale')
 for it in range(200):
  u=(y-X@b)/scale;adj=np.minimum(1,1.345/np.maximum(abs(u),1e-15));bn=fit(X,y,w*adj)[0]
  if np.max(abs(X@(bn-b)))<1e-9*(1+np.max(abs(X@b))):return bn,scale,it+1
  b=bn
 raise ValueError('Huber failed to converge')

def sensitivity(X,y,w,c,rng):
 q0=w/w.sum();out=[]
 for rho in [0.,.25,.5,1.]:
  if rho==0:
   v=float(c@fit(X,y,q0)[0]);out.append(dict(rho=rho,low=v,high=v,successes=2,attempts=2,max_violation=0.));continue
  lo=q0*np.exp(-rho);hi=q0*np.exp(rho)
  starts=[q0]
  for k in range(7):
   lp=linprog(rng.normal(size=len(y)),A_eq=np.ones((1,len(y))),b_eq=[1],bounds=list(zip(lo,hi)),method='highs')
   if lp.success:starts.append(.95*lp.x+.05*q0)
  vals=[];successful=0;violation=0;logs=[]
  for sign in [1.,-1.]:
   def fun(q):
    b=fit(X,y,q)[0];return sign*float(c@b)
   def jac(q):
    b,res,*_=fit(X,y,q);return sign*(X@np.linalg.solve(X.T@(q[:,None]*X),c))*res
   for start in starts:
    opt=minimize(fun,start,jac=jac,bounds=list(zip(lo,hi)),constraints={'type':'eq','fun':lambda q:q.sum()-1,'jac':lambda q:np.ones(len(q))},method='SLSQP',options={'maxiter':1000,'ftol':1e-10})
    v=max(abs(opt.x.sum()-1),float(np.max(lo-opt.x)),float(np.max(opt.x-hi)),0)
    violation=max(violation,v)
    logs.append(dict(direction="minimum" if sign==1 else "maximum",success=bool(opt.success),status=int(opt.status),message=str(opt.message),iterations=int(opt.nit),constraint_violation=float(v),attained_target=float(c@fit(X,y,opt.x)[0])))
    if opt.success and v<1e-7:successful+=1;vals.append(float(c@fit(X,y,opt.x)[0]))
  if not vals:raise ValueError('all sensitivity optimizations failed')
                                                                                          
  vals.append(float(c@fit(X,y,q0)[0]))
  out.append(dict(rho=rho,low=min(vals),high=max(vals),successes=successful,attempts=2*len(starts),max_violation=violation,runs=logs))
 return out

def dump_csv(path,rows):
 with open(path,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main(datafile,assignfile,outdir,reps):
 out=Path(outdir);out.mkdir(exist_ok=True,parents=True)
 rr=list(csv.reader(open(datafile, newline='', encoding='utf-8')));headers=rr.pop(0);rr=[[r[0]]+[float(v) for v in r[1:]] for r in rr]
 cells=[dict(zip(['cell','age','mileage','n','CO_mean','CO_sd','NMHC_mean','NMHC_sd','NOx_mean','NOx_sd'],r)) for r in rr if r[0]]
 assert len(cells)==18 and len(set(r['cell'] for r in cells))==18
 a=np.array([r['age'] for r in cells],float);m=np.array([r['mileage'] for r in cells],float);n=np.array([r['n'] for r in cells],float)
 assert np.all(n>1) and np.all(n==np.floor(n))
 assign=list(csv.reader(open(assignfile, newline='', encoding='utf-8')))[1:]
 assignments={c['cell']:sum(r[1]==c['cell'] for r in assign) for c in cells}
 summary={'seed':SEED,'reps':reps,'counts_sum':int(n.sum()),'assignments':assignments,'count_mismatch_cells':[r['cell'] for r in cells if assignments[r['cell']]!=r['n']],'hashes':{},'cells':cells,'fits':[],'sensitivity':{},'simulations':[]}
 for file in [datafile,assignfile]:summary['hashes'][Path(file).name]=hashlib.sha256(Path(file).read_bytes()).hexdigest()
 rng=np.random.default_rng(SEED)
 foldrows=[];coefrows=[];huberfoldrows=[]
 for pollutant in POL:
  y=np.array([r[pollutant+'_mean'] for r in cells]);sd=np.array([r[pollutant+'_sd'] for r in cells]);assert np.all(sd>0) and np.all(np.isfinite(y))
  for scheme in ['equal','inverse_sd','inverse_variance','count_informed']:
   w={'equal':np.ones(18),'inverse_sd':1/sd,'inverse_variance':1/sd**2,'count_informed':n/sd**2}[scheme]
   for kind in ['A','AM','CAT']:
    X=design(a,m,kind);c=ref(kind);b,res,h,cov,hc3,inv=fit(X,y,w);sse=float(np.sum(w*res**2));sst=float(np.sum(w*(y-np.average(y,weights=w))**2));pred=[];changes=[];fail=[]
    for g in range(18):
     mask=np.arange(18)!=g
     try:
      bg=fit(X[mask],y[mask],w[mask])[0];pred.append(float(X[g]@bg));changes.append(float(c@(bg-b)))
      foldrows.append(dict(pollutant=pollutant,scheme=scheme,model=kind,omitted=cells[g]['cell'],observed=float(y[g]),prediction=pred[-1],reference_change=changes[-1],leverage=float(h[g]),status='ok'))
     except (ValueError,np.linalg.LinAlgError) as e:fail.append(cells[g]['cell']);pred.append(float('nan'));changes.append(float('nan'))
    rr=dict(pollutant=pollutant,scheme=scheme,model=kind,reference=float(c@b),working_se=float(np.sqrt(c@cov@c)),hc3_se=float(np.sqrt(c@hc3@c)),SSE=sse,SST=sst,R2=1-sse/sst,adjusted_R2=1-sse/sst*17/(18-X.shape[1]),LOCO_MSE=float(np.mean((y-np.array(pred))**2)),max_deletion=float(max(abs(np.array(changes)))),max_deletion_cell=cells[int(np.argmax(abs(np.array(changes))))]['cell'],max_leverage=float(h.max()),max_leverage_cell=cells[int(h.argmax())]['cell'],N_eff=float(w.sum()**2/np.sum(w**2)),rank=int(np.linalg.matrix_rank(X)),coefficients=X.shape[1],failed_folds=len(fail))
    if scheme=='count_informed':rr['plugin_se']=float(np.sqrt(c@inv@c))
    summary['fits'].append(rr)
    for j,v in enumerate(b):coefrows.append(dict(pollutant=pollutant,scheme=scheme,model=kind,index=j,coefficient=float(v)))
  X=design(a,m,'AM');c=ref('AM');w=1/sd**2;b,scale,it=huber(X,y,w);pred=[];failed=[]
  for g in range(18):
   mask=np.arange(18)!=g
   try:
    pred.append(float(X[g]@huber(X[mask],y[mask],w[mask])[0]));huberfoldrows.append(dict(pollutant=pollutant,omitted=cells[g]['cell'],observed=float(y[g]),prediction=pred[-1]))
   except ValueError:failed.append(g);pred.append(float('nan'))
  summary.setdefault('huber',{})[pollutant]=dict(reference=float(c@b),scale=float(scale),iterations=it,LOCO_MSE=float(np.mean((y-pred)**2)),failed_folds=len(failed))
  summary['sensitivity'][pollutant]=sensitivity(X,y,w,c,rng)
 print('Empirical fits complete',flush=True)
                                                                                         
 X=design(a,m,'AM');C=design(a,m,'CAT');c=ref('AM');cc=ref('CAT');base=X@np.array([2.,.3,.2,.05]);n_unequal=np.array([20,80,30,25,100,20,40,90,25,30,80,20,40,100,30,20,80,25.]);hetero=np.linspace(1.,4.,18)
 methods=['equal_AM','dispersion_AM','oracle_AM','dispersion_CAT','Huber_AM']
 scenarios=['equal_counts_equal_sd','unequal_counts','unequal_sd','counts_and_sd','contaminated_cell','categorical_misspecification']
 simulation_replications=[]
 for scenario in scenarios:
  counts=np.full(18,50.) if scenario in ['equal_counts_equal_sd','unequal_sd'] else n_unequal
  sd=np.ones(18)*2 if scenario in ['equal_counts_equal_sd','unequal_counts'] else hetero
  truth=base+ (1.5*(m==2) if scenario=='categorical_misspecification' else 0)
  target=float(c@np.array([2.,.3,.2,.05]))+(1.5 if scenario=='categorical_misspecification' else 0)
  errors={method:[] for method in methods};failures={method:0 for method in methods}
  for rep in range(reps):
   y=truth+rng.normal(size=18)*sd/np.sqrt(counts)
   if scenario=='contaminated_cell':y[16]+=6*sd[16]/np.sqrt(counts[16])
   for method in methods:
    try:
     if method=='Huber_AM':bh=huber(X,y,1/sd**2)[0];val=c@bh
     else:
      mat=C if method=='dispersion_CAT' else X;cr=cc if method=='dispersion_CAT' else c
      w=np.ones(18) if method=='equal_AM' else counts/sd**2 if method=='oracle_AM' else 1/sd**2
      val=cr@fit(mat,y,w)[0]
     errors[method].append(float(val-target));simulation_replications.append(dict(scenario=scenario,replication=rep+1,method=method,prediction=float(val),target=float(target),error=float(val-target)))
    except (ValueError,np.linalg.LinAlgError):failures[method]+=1
  for method,e in errors.items():
   e=np.array(e);sq=e**2;rmse=np.sqrt(sq.mean());summary['simulations'].append(dict(scenario=scenario,method=method,replications=len(e),failures=failures[method],bias=float(e.mean()),RMSE=float(rmse),RMSE_MCSE=float(sq.std(ddof=1)/(2*rmse*np.sqrt(len(e))))))
  print('Simulation',scenario,'complete',flush=True)
 dump_csv(out/'simulation_replications.csv',simulation_replications)
 import scipy, matplotlib
 summary['software']={'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'openpyxl':openpyxl.__version__,'matplotlib':matplotlib.__version__}
 (out/'analysis_results.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
 dump_csv(out/'huber_folds.csv',huberfoldrows)
 dump_csv(out/'grouped_inputs.csv',cells);dump_csv(out/'cell_diagnostics.csv',foldrows);dump_csv(out/'coefficients.csv',coefrows);dump_csv(out/'simulation_results.csv',summary['simulations'])
                                                     
 for pi in POL:
  y=np.array([r[pi+'_mean'] for r in cells]);sd=np.array([r[pi+'_sd'] for r in cells]);X=design(a,m,'AM');w=1/sd**2;b,res,h,*_=fit(X,y,w)
  assert np.allclose(b,fit(X,y,100*w)[0],rtol=1e-10,atol=1e-10)
  assert abs(h.sum()-4)<1e-10
  for g in range(18):
   mask=np.arange(18)!=g;deleted=y[g]-X[g]@fit(X[mask],y[mask],w[mask])[0]
   assert abs(deleted-res[g]/(1-h[g]))<1e-8
 print('Checks passed: weight scale invariance, leverage trace, direct LOCO versus PRESS identity.',flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--data',default=str(Path(__file__).parent.parent/'data'/'grouped_inputs.csv'));ap.add_argument('--assignments',default=str(Path(__file__).parent.parent/'data'/'resampled_assignments.csv'));ap.add_argument('--out',default=str(Path(__file__).parent/'results'));ap.add_argument('--reps',type=int,default=1000);args=ap.parse_args();main(args.data,args.assignments,args.out,args.reps)
