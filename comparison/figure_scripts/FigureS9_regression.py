#!/usr/bin/env python
"""Figure 4S2 - Confound-REGRESSION control (audit 2026-07-06).
Answers directly: "if FC vs BG differences exist, are they just contamination / library-size / depth?"
For each dataset, the astrocyte and bergmann module-score delta (FC - BG) is recomputed after
regressing OUT the 5 contamination scores (Purkinje/MLI/granule/oligoWM/vascular) via OLS, and additionally
   library size (Hao, log10 nCount) or laminar depth (Zhuang, depth_um). The ML-indicator coefficient is the
contamination-adjusted delta. Suyama uses the canonical pre-specified 56-cell purity-gated set.
Outputs: comparison/tables/fig4S2_confound_regression.csv + FigureS9_regression.png"""
import os, numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from project_paths import ROOT as PROJECT_ROOT
from figstyle import use_style, DSL, panel_letter, save
use_style()
BASE=str(PROJECT_ROOT)
TAB=BASE+"/comparison/tables"; FIG=BASE+"/comparison/figures"
CONTAM=['contam_Purkinje','contam_MLI','contam_granule','contam_oligoWM','contam_vascular']
DSETS=[(DSL[k],k,e) for k,e in [('zhuang','depth_um'),('hao','nCount'),('suyama',None)]]

def ols_ml(y,X):
    beta,_,_,_=np.linalg.lstsq(X,y,rcond=None)
    resid=y-X@beta; n,k=X.shape; dof=max(n-k,1); s2=(resid@resid)/dof
    se=np.sqrt(s2*np.diag(np.linalg.inv(X.T@X)))
    t=beta[1]/se[1]; return beta[1],2*stats.t.sf(abs(t),dof)

def deltas(m,score,extra=None):
    d=m[m.group.isin(['FC','BG'])].copy()
    y=d[score].values.astype(float); ml=(d.group=='FC').astype(float).values
    covs=[d[c].values.astype(float) for c in CONTAM]+([d[e].values.astype(float) for e in extra] if extra else [])
    Z=[(c-c.mean())/(c.std() if c.std()>0 else 1) for c in covs]
    X=np.column_stack([np.ones(len(y)),ml]+Z)
    raw=y[ml==1].mean()-y[ml==0].mean(); b,p=ols_ml(y,X); return raw,b,p

rows=[]
for name,key,extra in DSETS:
    m=pd.read_parquet(f"{TAB}/{key}_astro_cells.parquet")
    m['group']=m['group'].replace({'FC':'FC','BG':'BG'})
    if extra=='nCount': m=m.copy(); m['nCount']=np.log10(m['nCount'].clip(lower=1))
    if extra=='depth_um': m=m.dropna(subset=['depth_um']).copy()
    nml=int((m.group=='FC').sum())
    for score in ['astrocyte','bergmann']:
        raw,adj,padj=deltas(m,score)
        rec=dict(dataset=name,score=score,n_FC=nml,raw_delta=round(raw,3),
                 contam_adj_delta=round(adj,3),contam_adj_p=padj)
        if extra:
            _,fadj,fp=deltas(m,score,extra=[extra])
            rec['full_adj_delta']=round(fadj,3); rec['full_adj_p']=fp; rec['extra_covar']=extra
        rows.append(rec)
R=pd.DataFrame(rows); R.to_csv(f"{TAB}/fig4S2_confound_regression.csv",index=False)
print(R.to_string(index=False))

# ---------------- FIGURE ----------------
def star(p): return '***' if p<1e-3 else '**' if p<1e-2 else '*' if p<0.05 else 'ns'
fig,axes=plt.subplots(1,2,figsize=(13,5.2),sharey=True)
names=[n for n,_,_ in DSETS]; xp=np.arange(len(names)); w=0.26
cols={'raw':'#bcbddc','contam':'#807dba','full':'#3f007d'}
for ax,score in zip(axes,['astrocyte','bergmann']):
    sub=R[R.score==score].set_index('dataset').reindex(names)
    raw=sub['raw_delta'].values
    cadj=sub['contam_adj_delta'].values
    fadj=sub['full_adj_delta'].values  # NaN for Suyama (no extra)
    ax.bar(xp-w,raw,w,color=cols['raw'],label='raw',edgecolor='k',lw=.4)
    ax.bar(xp,cadj,w,color=cols['contam'],label='contamination-adjusted (5 scores)',edgecolor='k',lw=.4)
    fvals=np.where(np.isnan(fadj),cadj,fadj)  # for Suyama full=contam
    ax.bar(xp+w,fvals,w,color=cols['full'],label='+ library-size (Stereo-seq) / depth (MERFISH)',edgecolor='k',lw=.4)
    for i in range(len(names)):
        p=sub['contam_adj_p'].values[i]
        ax.text(xp[i],cadj[i]+(0.02 if cadj[i]>=0 else -0.05),star(p),ha='center',fontsize=12.2,fontweight='bold')
    ax.axhline(0,c='k',lw=.8); ax.set_xticks(xp); ax.set_xticklabels(names,fontsize=12.2)
    ax.set_title(f'{score.capitalize()} — Δ score (FC − BG)',fontsize=14.9,fontweight='bold')
    ax.grid(axis='y',ls=':',alpha=.4)
axes[0].set_ylabel('Δ module score (z)',fontsize=13.5)
axes[0].legend(fontsize=10.1,loc='upper left',framealpha=.9)
panel_letter(axes[0],'a'); panel_letter(axes[1],'b')
fig.tight_layout()
save(fig, "FigureS9_regression")
plt.close(fig)
