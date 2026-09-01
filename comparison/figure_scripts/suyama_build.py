#!/usr/bin/env python
"""Build the canonical Suyama astroglia tables: FC / BG / parenchymal astrocytes.

The active FC set contains 56 cells retained by the symmetric purity gate
documented in HARMONIZATION_protocol.md.
"""
import os, json, sys, numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(__file__)); from panels import module_scores
from project_paths import ROOT as PROJECT_ROOT
BASE=str(PROJECT_ROOT); SU=BASE+"/datasets/suyama_2026_xenium"; TAB=BASE+"/comparison/tables"
dec=pd.read_csv(SU+"/results/fc_decisions_suyama.csv"); ml=set(dec[dec.decision=='yes']['cell_label'].astype(str))
lam=pd.read_csv(SU+"/results/tables/adult_bergmann_laminar_purkinje.csv"); lam['cell']=lam['cell'].astype(str)
pcl=set(lam[lam.bgclass.isin(['BG','BG_PCL'])]['cell'])-ml  # 09_adult_laminar.R renamed 'BG'->'BG_PCL'; accept both
cl2t=json.load(open(SU+"/results/tables/adult_cl2type.json"))
sa=pd.read_csv(SU+"/results/tables/adult_cells_scores.csv"); sa['cell']=sa['cell'].astype(str)
sa['ct']=sa['cl'].map(cl2t)
astro=set(sa[sa.ct=='Astrocyte']['cell'])-ml-pcl-set(lam.cell)   # parenchymal astro, not Bergmann
print("FC",len(ml),"BG",len(pcl),"OtherCBastro",len(astro))
ex=pd.read_csv(SU+"/results/tables/adult_panel_expr.csv"); ex['cell']=ex['cell'].astype(str); ex=ex.set_index('cell')
ex=ex/np.log(2)  # Seurat 'data' layer is ln(1+CP10k); exact conversion to log2(1+CP10k) — base homogenized 2026-08-17
want=ml|pcl|astro
E=ex.reindex([c for c in ex.index if c in want])
grp=np.where(E.index.isin(ml),'FC',np.where(E.index.isin(pcl),'BG','OtherCBastro'))
sco0,_=module_scores(E)
# --- CANONICAL FC SET: symmetric purity gate ---
# contam_Purkinje < 1 z and contam_granule < 1 z, computed by
# build_suyama_purity_gate.py. Stable list: results/tables/suyama_fc_gate.json.
gate56=set(str(c) for c in json.load(open(SU+"/results/tables/suyama_fc_gate.json")))
print(f"Suyama canonical FC = purity-gated set: {len(gate56)}")
ml=ml&gate56
want=ml|pcl|astro
E=ex.reindex([c for c in ex.index if c in want])
grp=np.where(E.index.isin(ml),'FC',np.where(E.index.isin(pcl),'BG','OtherCBastro'))
sco,used=module_scores(E)   # canonical 56-cell set
xy=sa.set_index('cell').reindex(E.index); dpk=lam.set_index('cell').reindex(E.index)['d_pkj']
meta=pd.DataFrame({'cell':E.index,'group':grp,'section':'Adult','animal':'adult',
                   'x':xy['x'].values,'y':xy['y'].values,'depth_um':dpk.values})
meta=pd.concat([meta.reset_index(drop=True),sco.reset_index(drop=True)],axis=1)
meta.to_parquet(TAB+"/suyama_astro_cells.parquet"); E.to_parquet(TAB+"/suyama_astro_expr.parquet")
print("groups:",pd.Series(grp).value_counts().to_dict())
print("genes used:",{k:len(v) for k,v in used.items()})
print("depth available:",meta['depth_um'].notna().sum(),"/",len(meta))
print("saved suyama_astro_cells/expr")
