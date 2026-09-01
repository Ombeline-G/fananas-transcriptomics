#!/usr/bin/env python
"""Build the canonical Hao astroglia table from Mouse 1 sections only.
FC = swipe-validated (section+label); BG = segmentation (capped); OtherCBastro = parenchymal astro (GCL/WM).
Expression = log2(1+CP10k), panel genes only (base homogenized across datasets 2026-08-17; was ln).
Stores depth (dist_to_PCL) + nCount."""
import os, re, sys, json, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
sys.path.insert(0,os.path.dirname(__file__)); from panels import module_scores, PANELS
from project_paths import HAO as HAO_PATH, TABLES
HAO=str(HAO_PATH); TAB=str(TABLES)
PG=sorted({g for gs in PANELS.values() for g in gs})
rng=np.random.default_rng(0)
def strip(s): return re.sub(r'-\d+$','',str(s))
dec=pd.read_csv(HAO+"/swipe_validation/fc_decisions_hao.csv")
dec=dec[dec.animal.eq('Mouse1')].copy()
ml=dec[dec.decision=='yes'].copy()
ml['cl']=ml.cell_label.astype(str).map(strip)
seg=pd.read_csv(HAO+"/layer_segmentation/BG_assignments_all.csv")
seg=seg[seg.section.astype(str).str.startswith('Mouse1_')].copy()
seg['cl']=seg.cell_label.astype(str).map(strip)
pclseg=seg[seg.region=='BG']
def secpath(sec):
    m=sec.split('_')[0]; sub='stereo_seq_mouse1' if m=='Mouse1' else 'stereo_seq_mouse2'
    return f"{HAO}/data_raw/{sub}/{sec}.h5ad"
secs=sorted(set(ml.section)|set(pclseg.section))
print("sections to scan:",len(secs))
ANN={'molecular layer':'ML','granular layer':'GCL','purkinje layer':'PCL','white matter':'WM'}
frames=[]; metas=[]
for k,sec in enumerate(secs):
    p=secpath(sec)
    if not os.path.exists(p): continue
    a=ad.read_h5ad(p); a.var_names_make_unique()
    idx=a.obs.index.astype(str); idxset=set(idx); strip2pos={}
    # map stripped label -> positions
    cl=np.array([strip(x) for x in idx])
    ann=a.obs['annotation'].astype(str).map(ANN).values if 'annotation' in a.obs else np.array(['?']*len(idx))
    tot=np.asarray(a.X.sum(1)).ravel().astype(float); tot[tot==0]=1
    pg=[g for g in PG if g in set(map(str,a.var_names))]
    Xp=a[:,pg].X; Xp=Xp.toarray() if sp.issparse(Xp) else np.asarray(Xp)
    norm=np.log1p(Xp/tot[:,None]*1e4)/np.log(2)  # log2(1+CP10k): identical selection (z-invariant), log2 units
    nd=pd.DataFrame(norm,columns=pg)
    # per-section astro/berg z for astro selection
    def zmean(genes):
        gg=[g for g in genes if g in pg]
        if not gg: return np.zeros(len(nd))
        Z=(nd[gg]-nd[gg].mean())/(nd[gg].std()+1e-9); return Z.mean(1).values
    az=zmean(PANELS['astrocyte']); bz=zmean(PANELS['bergmann'])
    mlset=set(ml[ml.section==sec]['cl']); pclset=set(pclseg[pclseg.section==sec]['cl'])
    is_ml=np.isin(cl,list(mlset)); is_pcl=np.isin(cl,list(pclset))&~is_ml
    pcl_pos=np.where(is_pcl)[0]
    if len(pcl_pos)>50: pcl_pos=rng.choice(pcl_pos,50,replace=False)
    astro_ok=(~is_ml)&(~is_pcl)&np.isin(ann,['GCL','WM'])&(az>0.5)&(az>bz)
    astro_pos=np.where(astro_ok)[0]
    if len(astro_pos)>30: astro_pos=rng.choice(astro_pos,30,replace=False)
    sel=np.concatenate([np.where(is_ml)[0],pcl_pos,astro_pos]).astype(int)
    if len(sel)==0: continue
    grp=np.array(['OtherCBastro']*len(sel),dtype=object)
    grp[np.isin(sel,np.where(is_ml)[0])]='FC'; grp[np.isin(sel,pcl_pos)]='BG'
    # depth
    dep=np.full(len(sel),np.nan)
    segsec=seg[seg.section==sec].set_index('cl')['dist_to_PCL_um']
    for j,ppos in enumerate(sel):
        lab=cl[ppos]
        if lab in segsec.index:
            v=segsec.loc[lab]; dep[j]=v if np.isscalar(v) else np.asarray(v)[0]
    decsec=ml[ml.section==sec].set_index('cl')['dist_to_PCL_um']
    for j,ppos in enumerate(sel):
        if grp[j]=='FC' and cl[ppos] in decsec.index: dep[j]=decsec.loc[cl[ppos]] if np.isscalar(decsec.loc[cl[ppos]]) else np.asarray(decsec.loc[cl[ppos]])[0]
    cid=[f"{sec}|{idx[pp]}" for pp in sel]
    frames.append(pd.DataFrame(norm[sel],columns=pg,index=cid))
    metas.append(pd.DataFrame({'cell':cid,'group':grp,'section':sec,'animal':sec.split('_')[0],
                               'x':a.obs['coor_x'].values[sel].astype(float),'y':a.obs['coor_y'].values[sel].astype(float),
                               'depth_um':dep,'nCount':tot[sel]}))
    if k%10==0: print(f"  [{k}/{len(secs)}] {sec} sel={len(sel)}")
E=pd.concat(frames).fillna(0.0); meta=pd.concat(metas,ignore_index=True)
E=E.loc[~E.index.duplicated()]; meta=meta.drop_duplicates('cell').set_index('cell').reindex(E.index).reset_index()
sco,used=module_scores(E)
meta=pd.concat([meta.reset_index(drop=True),sco.reset_index(drop=True)],axis=1)
meta.to_parquet(TAB+"/hao_astro_cells.parquet"); E.to_parquet(TAB+"/hao_astro_expr.parquet")
print("groups:",pd.Series(meta.group).value_counts().to_dict())
print("by mouse:",meta.groupby(['animal','group']).size().to_dict())
print("genes used:",{k:len(v) for k,v in used.items()})
print("depth available:",meta['depth_um'].notna().sum(),"/",len(meta))
print("saved hao_astro_cells/expr")
