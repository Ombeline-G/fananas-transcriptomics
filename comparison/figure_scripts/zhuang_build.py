#!/usr/bin/env python
"""Build Zhuang unified astroglia table: FC / BG / OtherCBastro (317 Astro-CB NN).
Four animals since 2026-08-18 (ABCA-1 integrated after swipe validation, review point 35):
cell_sets.json carries the validated FC (362) and confident BG (6,084) of ABCA-1/2/3/4.
depth_um: robust_layers pipeline (ABCA-2/3/4) or BG_assignments_ABCA1 distance (ABCA-1) —
same quantity (µm above the local Purkinje line), estimator recorded in Methods/CANONICAL.
Outputs comparison/tables/zhuang_astro_{cells.parquet, expr.parquet}."""
import os, json, sys, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
sys.path.insert(0,os.path.dirname(__file__)); from panels import module_scores
from project_paths import ROOT as PROJECT_ROOT
BASE=str(PROJECT_ROOT); ZH=BASE+"/datasets/zhang_zhuang_2023_merfish"; TAB=BASE+"/comparison/tables"
cs=json.load(open(ZH+"/validated_analysis/cell_sets.json"))
ml=set(map(str,cs['FC_validated'])); pcl=set(map(str,cs['BG']))
mv=pd.read_parquet(ZH+"/cb_cells_multivote.parquet"); mv.index=mv.index.astype(str)
rl=pd.read_parquet(BASE+"/comparison/tables/robust_layers_zhuang.parquet"); rl.index=rl.index.astype(str)
# sections that carry the validated BG cells -> restrict OtherCBastro to those
secs=set(mv.loc[mv.index.isin(ml|pcl),'brain_section_label'])
astro=set(mv[(mv.subclass=='317 Astro-CB NN')&(mv.brain_section_label.isin(secs))].index)
astro-=(ml|pcl)
print("sections:",len(secs)," FC",len(ml),"BG",len(pcl),"OtherCBastro",len(astro))
want=ml|pcl|astro
# load expression from clean animals
frames=[]
for an in ['Zhuang-ABCA-1','Zhuang-ABCA-2','Zhuang-ABCA-3','Zhuang-ABCA-4']:
    f=f"{ZH}/abc_cache/expression_matrices/{an}/20230830/{an}-log2.h5ad"
    a=ad.read_h5ad(f,backed='r'); m=a.obs.index.astype(str).isin(want)
    if m.sum()==0: continue
    sub=a[m].to_memory(); X=sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
    frames.append(pd.DataFrame(X,columns=sub.var['gene_symbol'].astype(str).values,index=sub.obs.index.astype(str)))
    print(" ",an,m.sum())
E=pd.concat(frames); E=E.loc[~E.index.duplicated()]
grp=np.where(E.index.isin(ml),'FC',np.where(E.index.isin(pcl),'BG','OtherCBastro'))
sco,used=module_scores(E)
meta=pd.DataFrame({'cell':E.index,'group':grp},index=E.index)
meta['section']=mv.reindex(E.index)['brain_section_label'].values
meta['animal']=mv.reindex(E.index)['animal'].values
meta['x']=mv.reindex(E.index)['x'].values; meta['y']=mv.reindex(E.index)['y'].values
meta['depth_um']=rl.reindex(E.index)['dP_um'].values
a1=pd.read_csv(ZH+"/layer_segmentation/BG_assignments_ABCA1.csv",dtype={'cell_label':str}).set_index('cell_label')
gap=meta['depth_um'].isna()&meta['cell'].isin(a1.index)
meta.loc[gap,'depth_um']=a1.reindex(meta.loc[gap,'cell'])['dist_to_PCL_um'].values
print("depth from ABCA-1 assignments:",int(gap.sum()))
meta['perp_norm']=mv.reindex(E.index)['perp_norm'].values
meta=pd.concat([meta.reset_index(drop=True),sco.reset_index(drop=True)],axis=1)
meta.to_parquet(TAB+"/zhuang_astro_cells.parquet")
E.to_parquet(TAB+"/zhuang_astro_expr.parquet")
print("groups:",pd.Series(grp).value_counts().to_dict())
print("genes used per panel:",{k:len(v) for k,v in used.items()})
print("depth available:",meta['depth_um'].notna().sum(),"/",len(meta))
print("saved zhuang_astro_cells/expr")
