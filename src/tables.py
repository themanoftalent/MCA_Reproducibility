import csv,json
import numpy as np
from analyze import dump_csv


def tables(root):
    out=root/'results';s=json.loads((out/'analysis_results.json').read_text())
    def read(name):
        with (out/name).open() as f:return list(csv.DictReader(f))
    def f(p,model,scheme='inverse_variance'):return next(r for r in s['fits'] if r['pollutant']==p and r['model']==model and r['scheme']==scheme)
    dump_csv(out/'table01_grouped_inputs.csv',s['cells'])
    dump_csv(out/'table02_model_comparison.csv',[{k:r[k] for k in ['pollutant','model','reference','working_se','hc3_se','R2','LOCO_MSE']} for r in s['fits'] if r['scheme']=='inverse_variance'])
    specs={'A':'Age only','AM':'Age, ordinal mileage, age by mileage','CAT':'Separate age lines for three mileage classes'}
    rows=[]
    for model in specs:
        rr=[f(p,model) for p in ['CO','NMHC','NOx']]
        rows.append(dict(model=model,specification=specs[model],coefficients=rr[0]['coefficients'],R2_min=min(r['R2'] for r in rr),R2_max=max(r['R2'] for r in rr),CO_LOCO_MSE=rr[0]['LOCO_MSE'],NMHC_LOCO_MSE=rr[1]['LOCO_MSE'],NOx_LOCO_MSE=rr[2]['LOCO_MSE']))
    dump_csv(out/'table03_specifications.csv',rows)
    dump_csv(out/'table04_weighting.csv',[{k:r[k] for k in ['pollutant','scheme','reference','LOCO_MSE','N_eff']} for r in s['fits'] if r['model']=='AM'])
    rows=[]
    for p in ['CO','NMHC','NOx']:
        r=f(p,'AM');h=s['huber'][p]
        rows.append(dict(pollutant=p,max_deletion=r['max_deletion'],omitted_cell=r['max_deletion_cell'],max_leverage=r['max_leverage'],Huber_reference=h['reference'],Huber_LOCO_MSE=h['LOCO_MSE']))
    dump_csv(out/'table05_influence.csv',rows)
    dump_csv(out/'table06_sensitivity.csv',[dict(pollutant=p,**{k:v for k,v in r.items() if k!='runs'}) for p,rs in s['sensitivity'].items() for r in rs])
    dump_csv(out/'sensitivity_solver_runs.csv',[dict(pollutant=p,rho=r['rho'],start=i+1,**v) for p,rs in s['sensitivity'].items() for r in rs for i,v in enumerate(r.get('runs',[]))])
    dump_csv(out/'table07_original_simulation.csv',s['simulations'])
    folds=read('age_level_folds.csv');rows=[]
    for p in ['CO','NMHC','NOx']:
        row={'pollutant':p}
        for model in specs:row[model+'_MSE']=float(np.mean([float(r['squared_error']) for r in folds if r['pollutant']==p and r['model']==model]))
        for age in [.5,5]:row['CAT_age_'+str(age)+'_MSE']=float(np.mean([float(r['squared_error']) for r in folds if r['pollutant']==p and r['model']=='CAT' and float(r['held_age'])==age]))
        rows.append(row)
    dump_csv(out/'table08_age_validation.csv',rows)
    dump_csv(out/'table09_covariance.csv',read('covariance_comparison.csv'))
    ext=read('extended_simulations.csv')
    for num,mode in [(10,'estimated'),(11,'known')]:dump_csv(out/f'table{num:02}_simulation_{mode}.csv',[r for r in ext if r['dispersion']==mode])
