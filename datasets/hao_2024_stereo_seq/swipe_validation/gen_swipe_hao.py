#!/usr/bin/env python
"""Build Hao swipe data from the segmentation. Candidates = FC (soma-certain, >=45um,
declumped). Background cells colored by Hao's native layer annotation. Reads kept_sections.json
if present, else all sections with FC."""
import os, json, glob, numpy as np, pandas as pd
import anndata as ad
from pathlib import Path
HAO = str(Path(__file__).resolve().parents[1])
RAWDIR={'Mouse1':HAO+"/data_raw/stereo_seq_mouse1",'Mouse2':HAO+"/data_raw/stereo_seq_mouse2"}
LS=HAO+"/layer_segmentation"; OUT=HAO+"/swipe_validation"
RCODE={'molecular layer':0,'granular layer':1,'white matter':2,'purkinje layer':3}
A=pd.read_csv(LS+"/BG_assignments_all.csv")
keptf=LS+"/kept_sections.json"
KEEP=json.load(open(keptf)) if os.path.exists(keptf) else sorted(A[A['region']=='FC']['section'].unique())
KEEP=sorted(KEEP)   # slices in order (Mouse1_T167..T198, then Mouse2_T345..)
print("sections:",len(KEEP))
sections={}; cards=[]
for sec in KEEP:
    mouse='Mouse2' if sec.startswith('Mouse2') else 'Mouse1'
    f=f"{RAWDIR[mouse]}/{sec}.h5ad"
    if not os.path.exists(f): continue
    o=ad.read_h5ad(f, backed='r').obs
    ann=o['annotation'].astype(str).map(RCODE).fillna(1).astype(int).values
    x=o['coor_x'].values.astype(float); y=o['coor_y'].values.astype(float)
    seg=A[A['section']==sec]
    bg=set(seg[seg['region']=='BG']['cell_label'].astype(str))
    ids=o.index.astype(str).values
    typ=np.where(np.isin(ids,list(bg)),2,0)
    cand=seg[seg['region']=='FC']
    iscand=np.isin(ids,cand['cell_label'].astype(str).values)
    keep=(ann==3)|(typ>0)|iscand
    rest=np.where(~keep)[0]; cap=9000; nrest=max(0,cap-int(keep.sum()))
    rng=np.random.default_rng(0)
    if len(rest)>nrest: rest=rng.choice(rest,nrest,replace=False)
    sel=np.unique(np.concatenate([np.where(keep)[0],rest]))
    cells=[[round(float(x[k]),1),round(float(y[k]),1),int(ann[k]),int(typ[k])] for k in sel]
    sections[sec]=dict(cells=cells,xmin=round(float(x.min()),1),xmax=round(float(x.max()),1),
                       ymin=round(float(y.min()),1),ymax=round(float(y.max()),1))
    for _,r in cand.iterrows():
        cards.append(dict(id=str(r['cell_label']),dataset="Hao",section=sec,animal=mouse,
            cx=round(float(r['x']),1),cy=round(float(r['y']),1),
            conf=round(float(r['confidence']),2),dist_pcl=round(float(r['dist_to_PCL_um']),0),support=0))
    print(f"  {sec}: cells={len(cells)} FC={len(cand)}")
cards.sort(key=lambda c:(c['section'], -c['conf']))   # slices in order; within a slice by confidence
data=dict(meta=dict(dataset="Hao (segmentation)",zoom_um=250,n=len(cards),
                    source="SURE FC >=45um from PCL, declumped; masks from native annotation"),
          sections=sections,cards=cards)
open(OUT+"/swipe_data.js","w",encoding="utf-8").write("window.SWIPE_DATA="+json.dumps(data)+";")
open(OUT+"/swipe_data.json","w",encoding="utf-8").write(json.dumps(data))
print("candidates:",len(cards)," js MB:",round(os.path.getsize(OUT+'/swipe_data.js')/1e6,2))
