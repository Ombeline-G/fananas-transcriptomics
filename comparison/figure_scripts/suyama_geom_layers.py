#!/usr/bin/env python
"""Suyama adult cerebellar layers -> results/tables/suyama_geom_layers.csv (used by Figure 2 & 3).
Uses the shared signed-distance layering (see geom_layers.py): three strictly parallel bands
GCL | PCL-buffer | ML built from the granule-mass boundary, so the ML/GCL boundary is a clean line
that neither side crosses, and the ML/PCL contain 0 granule cells by construction.
"""
import json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from geom_layers import assign_layers
from project_paths import ROOT as PROJECT_ROOT
BASE=str(PROJECT_ROOT)
SU=BASE+"/datasets/suyama_2026_xenium"
cl2t=json.load(open(SU+"/results/tables/adult_cl2type.json"))
sa=pd.read_csv(SU+"/results/tables/adult_cells_scores.csv"); sa['cell']=sa['cell'].astype(str)
sa['ct']=sa['cl'].map(cl2t).fillna('Other')
lay,sdist=assign_layers(sa.x,sa.y,sa.ct=='Granule',sa.ct=='Purkinje',sa.ct=='Oligodendrocyte',umpp=1.0)
sa['layer_geom']=lay; sa['sdist_um']=sdist
print(pd.Series(lay).value_counts().to_dict())
print(f'granule cells inside ML: {int(((sa.ct=="Granule")&(lay=="ML")).sum())} | inside PCL: {int(((sa.ct=="Granule")&(lay=="PCL")).sum())} (target 0/0)')

# diagnostic (full + zoom on one folium)
COL={'ML':'#9ccbe0','PCL':'#e8894a','GCL':'#b89ad0','WM':'#9aa0a6','none':'#f0f0f0'}
fig,ax=plt.subplots(1,2,figsize=(15,8))
for L in ['none','WM','GCL','ML','PCL']:
    m=lay==L; ax[0].scatter(sa.x[m],sa.y[m],s=4,c=COL[L],linewidths=0,label=f'{L} ({m.sum()})')
ax[0].set_aspect('equal'); ax[0].invert_yaxis(); ax[0].legend(markerscale=3); ax[0].set_title('Suyama layers (granule-mass signed distance)')
x0,x1,y0,y1=8850,9550,900,1700
for L in ['none','WM','GCL','ML','PCL']:
    m=(lay==L)&(sa.x.between(x0,x1))&(sa.y.between(y0,y1)); ax[1].scatter(sa.x[m],sa.y[m],s=22,c=COL[L],linewidths=0)
mg=(sa.ct=='Granule')&(sa.x.between(x0,x1))&(sa.y.between(y0,y1)); ax[1].scatter(sa.x[mg],sa.y[mg],s=5,c='k',label='granule cell')
ax[1].set_aspect('equal'); ax[1].invert_yaxis(); ax[1].legend(); ax[1].set_title('zoom: black granule cells only in GCL; ML/GCL separated by PCL band')
fig.savefig(BASE+"/comparison/figures/SuppFig_suyama_geom_layers.png",dpi=110,bbox_inches='tight'); print('saved diagnostic')
sa[['cell','x','y','ct','layer_geom','sdist_um']].to_csv(SU+"/results/tables/suyama_geom_layers.csv",index=False)
print("saved suyama_geom_layers.csv")
