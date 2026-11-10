import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager


def figures(root):
    available={f.name for f in font_manager.fontManager.ttflist}
    family=next((name for name in ['Times New Roman','Nimbus Roman','Liberation Serif'] if name in available),'DejaVu Serif')
    font=font_manager.FontProperties(family=family)
    plt.rcParams.update({'font.family':font.get_family(),'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300,'pdf.fonttype':42,'ps.fonttype':42})
    s=json.loads((root/'results/analysis_results.json').read_text());out=root/'figures'
    def save(fig,name):
        fig.savefig(out/(name+'.png'),dpi=300,metadata={'Software':'Matplotlib'})
        fig.savefig(out/(name+'.pdf'),metadata={'Creator':'Matplotlib','CreationDate':None,'ModDate':None})
        plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(11,3.7),constrained_layout=True)
    for ax,p in zip(axes,['CO','NMHC','NOx']):
        for model,color in [('A','#666666'),('AM','#2677A7'),('CAT','#B75425')]:
            r=next(r for r in s['fits'] if r['pollutant']==p and r['scheme']=='inverse_variance' and r['model']==model)
            ax.bar(model,np.sqrt(r['LOCO_MSE']),color=color)
        ax.set_title(p);ax.set_ylabel('Held-out cell RMSE (mg/km)');ax.set_ylim(bottom=0)
    save(fig,'figure1')
    fig,axes=plt.subplots(2,3,figsize=(12,7.5),constrained_layout=True)
    labels=['Equal counts and SD','Unequal counts','Unequal SD','Unequal counts and SD','Contaminated reference cell','Category-2 shift']
    for ax,sc,label in zip(axes.flat,dict.fromkeys(r['scenario'] for r in s['simulations']),labels):
        rows=[r for r in s['simulations'] if r['scenario']==sc]
        ax.bar(['Equal','Disp.','Oracle','CAT','Huber'],[r['RMSE'] for r in rows],yerr=[1.96*r['RMSE_MCSE'] for r in rows],color=['#666666','#2677A7','#498A56','#B75425','#8E66A1'],capsize=2)
        ax.set_title(label);ax.set_ylabel('Reference RMSE (simulation units)');ax.tick_params(axis='x',rotation=30);ax.set_ylim(bottom=0)
    save(fig,'figure2')
    with (root/'results/fitted_curves.csv').open() as f:curves=list(csv.DictReader(f))
    fig,axes=plt.subplots(1,3,figsize=(12,4),constrained_layout=True)
    for ax,p in zip(axes,['CO','NMHC','NOx']):
        for cl,color in [(1,'#2677A7'),(2,'#B75425'),(3,'#498A56')]:
            rr=[r for r in curves if r['pollutant']==p and int(r['mileage'])==cl]
            ax.plot([float(r['age']) for r in rr],[float(r['fitted_mean']) for r in rr],color=color,label='Class '+str(cl))
            cc=[r for r in s['cells'] if r['mileage']==cl]
            ax.scatter([r['age'] for r in cc],[r[p+'_mean'] for r in cc],color=color,s=22)
        ax.set_title(p);ax.set_xlabel('Supplied age (years)');ax.set_ylabel('Supplied mean (mg/km)');ax.legend(frameon=False)
    save(fig,'figure3')
