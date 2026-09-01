#!/usr/bin/env python
"""Figure 4S - Confound controls (from adversarial review):
A. Directional contamination delta (Purkinje & granule, FC-BG) per dataset -> sign-flip.
B. Hao library-size confound: nCount by group + bergmann-score vs log10(nCount).
C. Distance-ratio bootstrap CI from the canonical identity-metrics table; >1 = closer to BG.
Also adds a spatial-covariance Purkinje flag column to Figure6_signature_table.csv."""
import os, sys, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, pearsonr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
sys.path.insert(0,os.path.dirname(__file__)); from panels import PANELS
from figstyle import use_style, GROUP, DSL, panel_letter, save
use_style()
from project_paths import ROOT as PROJECT_ROOT
BASE=str(PROJECT_ROOT); TAB=BASE+"/comparison/tables"; FIG=BASE+"/comparison/figures"
DSETS=[(DSL[k],k) for k in ('zhuang','hao','suyama')]
GCOL=GROUP
DATA={}
for name,key in DSETS:
    m=pd.read_parquet(f"{TAB}/{key}_astro_cells.parquet")
    m['group']=m['group'].replace({'FC':'FC','BG':'BG'})
    if 'cell' not in m.columns: m=m.rename(columns={'index':'cell'})
    E=pd.read_parquet(f"{TAB}/{key}_astro_expr.parquet").reindex(m['cell'].values)
    DATA[key]=dict(name=name,m=m,E=E)

# ---- A: directional contamination delta ----
contA=[]
for name,key in DSETS:
    m=DATA[key]['m']
    for sc in ['contam_Purkinje','contam_granule','contam_MLI','contam_vascular']:
        a=m[m.group=='FC'][sc]; b=m[m.group=='BG'][sc]
        try:_,p=mannwhitneyu(a,b)
        except Exception:p=1.0
        contA.append(dict(dataset=name,score=sc,delta=float(a.mean()-b.mean()),p=p))
contA=pd.DataFrame(contA); contA.to_csv(f"{TAB}/fig4S_contam_delta.csv",index=False)

# ---- C: bootstrap distance ratio (single reproducible builder) ----
boot_table=pd.read_csv(f"{TAB}/fig_bootstrap_ratios.csv").set_index('dataset')
boot={key:(boot_table.loc[key,'med'],boot_table.loc[key,'lo'],boot_table.loc[key,'hi']) for _,key in DSETS}

# ---- spatial-covariance Purkinje flag for signature table ----
sig=pd.read_csv(f"{TAB}/Figure6_signature_table.csv")
covar={}
for name,key in DSETS:
    m=DATA[key]['m']; E=DATA[key]['E']; berg=m.group.isin(['FC','BG']).values
    cp=m['contam_Purkinje'].values[berg]
    for gname in E.columns:
        v=E[gname].values[berg]
        if v.std()>0 and cp.std()>0:
            r=pearsonr(v,cp)[0]; covar.setdefault(gname,[]).append(abs(r))
def pflag(g):
    rs=covar.get(g,[]);
    return f'r={max(rs):.2f}' if rs and max(rs)>0.4 else ''
sig['purkinje_spatial_covar']=sig['Gene'].map(pflag)
sig.to_csv(f"{TAB}/Figure6_signature_table.csv",index=False)
print("genes newly covar-flagged:",sorted({g for g in sig[sig.purkinje_spatial_covar!='']['Gene']}))

# ---- D: metabolic pathway activity FC vs BG (Hao whole-tx + Suyama 5K; Zhuang panel lacks metabolic genes) ----
SU=BASE+"/datasets/suyama_2026_xenium"
METAB={
 'Glutamate–glutamine cycle':['Slc1a3','Slc1a2','Glul','Gls','Slc38a1','Slc38a3','Got1','Glud1'],
 'Glycolysis':['Hk1','Hk2','Gpi1','Pfkl','Aldoa','Aldoc','Gapdh','Pgk1','Eno1','Pkm','Ldha'],
 'Lactate shuttle (MCT)':['Ldha','Ldhb','Slc16a1','Slc16a3','Slc16a4'],
 'Glycogen metabolism':['Gys1','Pygb','Pygl','Ugp2','Gbe1','Agl','Ppp1r3c'],
 'OXPHOS / mitochondria':['Ndufa4','Ndufs1','Sdha','Sdhb','Uqcrc1','Uqcrc2','Cox4i1','Cox5a','Cox6c','Atp5f1a','Atp5f1b'],
 'Fatty-acid / lipid':['Fabp7','Apoe','Cd36','Acsbg1','Elovl2','Acsl3','Cpt1a','Fasn','Scd2','Dbi'],
 'Cholesterol synthesis':['Hmgcr','Hmgcs1','Sqle','Dhcr7','Fdps','Idi1','Msmo1','Insig1','Lss'],
 'Glutathione / antioxidant':['Gclc','Gclm','Gss','Gpx1','Gpx4','Gsr','Sod1','Sod2','Cat','Prdx6'],
 'Iron / metal handling':['Fth1','Ftl1','Mt1','Mt2','Cp','Slc40a1','Tfrc'],
 'TCA / ketone support':['Bdh1','Oxct1','Acat1','Idh1','Idh2','Aco2','Cs','Mdh1'],
}
hstat=pd.read_csv(f"{TAB}/hao_bg_gene_stats.csv").pivot_table(index='gene',columns='group',values='mean')
hstat=hstat.rename(columns={"FC":"FC","BG":"BG"})
sexp=pd.read_csv(SU+"/results/tables/bergmann_full_expr.csv")
sexp['cell']=sexp['cell'].astype(str)
sexp["group"]=sexp["group"].replace({"BG_PCL":"BG"})   # renamed by 09_adult_laminar.R (2026-07-10)
# CANONICAL FC SET (étape 3, 2026-08-10): a-priori purity gate 73->56
import json as _json
_gate56=set(str(c) for c in _json.load(open(SU+"/results/tables/suyama_fc_gate.json")))
sexp=sexp[(sexp['group']=='BG')|((sexp['group']=='FC')&(sexp['cell'].isin(_gate56)))].reset_index(drop=True)
sgm=sexp.groupby('group').mean(numeric_only=True).T          # genes x {FC,BG}
GM={'Hao':hstat,'Suyama':sgm}
PATHS=list(METAB.keys()); MCOLS=[('Hao','FC'),('Hao','BG'),('Suyama','FC'),('Suyama','BG')]
Mmet=np.full((len(PATHS),4),np.nan); ncov=np.zeros((len(PATHS),2),int); mstar=['']*len(PATHS)
for pi,p in enumerate(PATHS):
    for ci,(dk,grp) in enumerate(MCOLS):
        gm=GM[dk]; g=[x for x in METAB[p] if x in gm.index]
        if g: Mmet[pi,ci]=gm.loc[g,grp].mean(); ncov[pi,0 if dk=='Hao' else 1]=len(g)
    gs2=[x for x in METAB[p] if x in sexp.columns]            # Suyama per-cell significance (has cells)
    if gs2:
        a=sexp[sexp.group=='FC'][gs2].mean(1); b=sexp[sexp.group=='BG'][gs2].mean(1)
        try:
            _,pv=mannwhitneyu(a,b); mstar[pi]='*' if pv<0.05 else 'ns'
        except Exception: pass
# scale each dataset block to its own max -> relative activity 0..1 (presence; ML vs PCL comparable within dataset)
Mrel=Mmet.copy()
for a,b in [(0,2),(2,4)]:
    sub=Mmet[:,a:b]; mx=np.nanmax(sub); Mrel[:,a:b]=sub/(mx+1e-9)
pd.DataFrame(Mmet,index=PATHS,columns=[f'{d}_{g}' for d,g in MCOLS]).to_csv(f"{TAB}/fig4S_metabolic.csv")

# ================= FIGURE =================
fig=plt.figure(figsize=(15,4.6)); gs=GridSpec(1,3,figure=fig,wspace=0.34,top=0.9,bottom=0.14)
# A
ax=fig.add_subplot(gs[0,0]); W=0.2; scs=['contam_Purkinje','contam_granule','contam_MLI','contam_vascular']
for i,sc in enumerate(scs):
    vals=[contA[(contA.dataset==nm)&(contA.score==sc)].delta.values[0] for nm,_ in DSETS]
    ax.bar(np.arange(3)+i*W,vals,W,label=sc.replace('contam_',''))
ax.axhline(0,c='k',lw=.8); ax.set_xticks(np.arange(3)+1.5*W); ax.set_xticklabels(['MERFISH','Stereo-seq','Xenium'],fontsize=9)
ax.set_ylabel('contamination Δ (FC − BG)',fontsize=12.2); ax.legend(fontsize=9.5,ncol=2)
ax.set_title('Directional contamination',fontsize=14.9,loc='left',fontweight='bold'); panel_letter(ax,'a')
# B Hao library size
axb=fig.add_subplot(gs[0,1]); mh=DATA['hao']['m']
if 'nCount' in mh.columns:
    data=[mh[mh.group==g]['nCount'].values for g in GCOL]
    bp=axb.boxplot([d[~np.isnan(d)] for d in data],labels=['FC','BG','Parenchymal'],showfliers=False,patch_artist=True)
    for patch,g in zip(bp['boxes'],GCOL): patch.set_facecolor(GCOL[g]); patch.set_alpha(0.5)
    axb.set_ylabel('nCount (transcripts/cell)',fontsize=12.2)
    r=np.corrcoef(np.log10(mh['nCount'].clip(lower=1)),mh['bergmann'])[0,1]
    axb.set_title('Stereo-seq library-size confound',fontsize=14.9,loc='left',fontweight='bold')
    axb.tick_params(labelsize=10.8)
panel_letter(axb,'b')
# C bootstrap ratio
axc=fig.add_subplot(gs[0,2])
for i,(name,key) in enumerate(DSETS):
    med,lo,hi=boot[key]
    axc.errorbar(med,i,xerr=[[med-lo],[hi-med]],fmt='o',color='#333',capsize=4,ms=7)
    axc.text(med,i+0.18,f'{med:.2f} [{lo:.2f},{hi:.2f}]',ha='center',fontsize=10.8)
axc.axvline(1.0,c='r',ls='--',lw=1); axc.set_yticks(range(3))
# noms courts : l'étiquette longue débordait sur le panneau voisin
axc.set_yticklabels([n.split(' ')[0] for n,_ in DSETS],fontsize=12.2)
axc.set_xlabel('d(FC, other) / d(FC, BG)',fontsize=12.2); axc.set_ylim(-0.5,2.5)
axc.set_title('Distance-ratio bootstrap',fontsize=14.9,loc='left',fontweight='bold'); panel_letter(axc,'c')
save(fig, "FigureS3_confounds")
plt.close(fig)
print("\ncontam_Purkinje delta (ML-PCL):",{nm:round(contA[(contA.dataset==nm)&(contA.score=='contam_Purkinje')].delta.values[0],2) for nm,_ in DSETS})
print("bootstrap distance ratio (med [95CI]):",{k:tuple(round(x,2) for x in boot[k]) for k in boot})
