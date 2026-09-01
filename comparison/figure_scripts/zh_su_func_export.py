#!/usr/bin/env python
"""Export Zhuang + Suyama FUNCTIONAL-gene expression for FC/BG cells -> {zhuang,suyama}_func_expr.parquet"""
import os,sys,json,numpy as np,pandas as pd,anndata as ad,scipy.sparse as sp
sys.path.insert(0,os.path.dirname(__file__)); from func_panels import ALL_FUNC
from project_paths import ROOT as PROJECT_ROOT
BASE=str(PROJECT_ROOT); ZH=BASE+"/datasets/zhang_zhuang_2023_merfish"; SU=BASE+"/datasets/suyama_2026_xenium"; TAB=BASE+"/comparison/tables"

# ---- Zhuang ----
cs=json.load(open(ZH+"/validated_analysis/cell_sets.json")); ml=set(map(str,cs['FC_validated'])); pcl=set(map(str,cs['BG']))
want=ml|pcl; frames=[]
for an in ['Zhuang-ABCA-1','Zhuang-ABCA-2','Zhuang-ABCA-3','Zhuang-ABCA-4']:
    a=ad.read_h5ad(f"{ZH}/abc_cache/expression_matrices/{an}/20230830/{an}-log2.h5ad",backed='r')
    m=a.obs.index.astype(str).isin(want)
    if not m.sum(): continue
    sub=a[m].to_memory(); X=sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
    df=pd.DataFrame(X,columns=sub.var['gene_symbol'].astype(str).values,index=sub.obs.index.astype(str))
    fg=[g for g in ALL_FUNC if g in df.columns]; frames.append(df[fg])
E=pd.concat(frames); E=E.loc[~E.index.duplicated()]
grp=np.where(E.index.isin(ml),'FC','BG')
pd.concat([pd.DataFrame({'cell':E.index,'group':grp}).reset_index(drop=True),E.reset_index(drop=True)],axis=1).to_parquet(TAB+"/zhuang_func_expr.parquet")
print("Zhuang func:",E.shape,pd.Series(grp).value_counts().to_dict(),"genes",len([g for g in ALL_FUNC if g in E.columns]),"/",len(ALL_FUNC))

# ---- Suyama (from bergmann_full_expr 5006 ; FC = set canonique gate) ----
e=pd.read_csv(SU+"/results/tables/bergmann_full_expr.csv")
gate=set(map(str,json.load(open(SU+"/results/tables/suyama_fc_gate.json"))))
e['cell']=e['cell'].astype(str)
e['group']=e['group'].replace({'BG_PCL':'BG'})
e=e[(e.group=='BG')|((e.group=='FC')&(e.cell.isin(gate)))]
fg=[g for g in ALL_FUNC if g in e.columns]
out=e[['cell','group']+fg].copy()
out[fg]=out[fg]/np.log(2)  # Seurat ln(1+CP10k) -> log2(1+CP10k), exact — base homogenized 2026-08-17
out.to_parquet(TAB+"/suyama_func_expr.parquet")
print("Suyama func:",out.shape,out.group.value_counts().to_dict(),"genes",len(fg),"/",len(ALL_FUNC))
print("MISSING Suyama:",[g for g in ALL_FUNC if g not in e.columns])
