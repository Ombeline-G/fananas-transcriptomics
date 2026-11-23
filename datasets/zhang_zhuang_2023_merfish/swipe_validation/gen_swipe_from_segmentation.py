#!/usr/bin/env python
"""Build swipe data from the GCL-anchored segmentation (kept slices).
Candidates = confident FC (>=45um from PCL, declumped). Background = segmentation
masks (per-cell region looked up from the saved mask rasters) so the layers are clean."""
import os, sys, json, numpy as np, pandas as pd
from pathlib import Path
ZH = str(Path(__file__).resolve().parents[1])
LS=ZH+"/layer_segmentation"; OUT=ZH+"/swipe_validation"
allcb=pd.read_parquet(f"{ZH}/cb_cells_multivote.parquet")
# default: the 20 kept ABCA-2/3/4 sections (240 candidates). `abca1` argument: the fourth
# animal (review point 35) — candidates from FC_candidates_ABCA1_for_swipe.csv (406, already
# >=45um + declumped), BG points from BG_assignments_ABCA1.csv, same masks/section format.
# `band35` argument (review S. Alon, November 2026): the Bergmann-class ML somata lying 35-45 um
# above the PCL in the 20 kept sections + the 14 swiped ABCA-1 sections, i.e. the band a 35-um
# margin would add to the candidate pool. Same rule as the segmentation, margin lowered to 35 um:
# ML Bergmann somata >=35 um, declumped at that margin, of which only the 35-45 um band is new
# (the >=45 um candidates were already reviewed). Confidence = the segmentation's FC formula
# evaluated with the 35-um margin. Candidates written to FC_candidates_band35_45_for_swipe.csv.
MODE = sys.argv[1].lower() if len(sys.argv) > 1 else "kept"
SOURCE = "GCL-anchored segmentation, FC>=45um from PCL, declumped"
if MODE == "band35":
    from scipy.spatial import cKDTree
    MARGIN, DECLUMP_R, DECLUMP_NB = 35.0, 0.035, 2      # um, mm, n (segment_cerebellar_layers.py)
    kept=pd.read_csv(LS+"/BG_assignments_kept.csv")
    a1=pd.read_csv(LS+"/BG_assignments_ABCA1.csv")
    a1=a1[a1['section'].isin(set(pd.read_csv(LS+"/FC_candidates_ABCA1_for_swipe.csv")['section']))]
    A=pd.concat([kept,a1],ignore_index=True).drop_duplicates('cell_label').set_index('cell_label')
    amb=A['region']=='Ambiguous_BG'
    ml=A[(A['region']=='FC')|(amb&(A['dist_to_PCL_um']>=45))|(amb&(A['confidence']>0.30)&(A['dist_to_PCL_um']<45))]
    ml=ml[ml['dist_to_PCL_um']>=MARGIN]
    ok=[]
    for _,g in ml.groupby('section'):
        xy=g[['x','y']].values; tt=cKDTree(xy)
        nb=np.array([len(tt.query_ball_point(q,DECLUMP_R))-1 for q in xy])
        ok+=list(g.index[nb<=DECLUMP_NB])
    CAND=ml.loc[ok]; CAND=CAND[CAND['dist_to_PCL_um']<45].copy()
    depth=(CAND['dist_to_PCL_um']-MARGIN)/1000.0
    base=np.minimum(1.0,0.5+depth/0.05)
    geomf=np.where(CAND['purk_support']>=2,np.minimum(1.0,0.6+0.1*np.minimum(CAND['purk_support'],4)),0.6)
    CAND['confidence']=np.round(base*geomf,2)
    CAND.reset_index().to_csv(LS+"/FC_candidates_band35_45_for_swipe.csv",index=False)
    KEEP=sorted(set(CAND['section']))
    TAG, DSLABEL = "band35", "Zhuang ABCA-1 to -4, 35-45 um band (segmentation)"
    SOURCE = "GCL-anchored segmentation, Bergmann-class ML somata 35-45um from PCL, declumped at the 35-um margin"
elif MODE == "abca1":
    A=pd.read_csv(LS+"/BG_assignments_ABCA1.csv").set_index('cell_label')
    CAND=pd.read_csv(LS+"/FC_candidates_ABCA1_for_swipe.csv").set_index('cell_label')
    KEEP=sorted(set(CAND['section']))
    TAG, DSLABEL = "abca1", "Zhuang ABCA-1 (segmentation)"
else:
    A=pd.read_csv(LS+"/BG_assignments_kept.csv").set_index('cell_label')
    CAND=None
    KEEP=json.load(open(LS+"/kept_sections.json"))
    TAG, DSLABEL = "zhuang", "Zhuang (segmentation)"
RCODE_PCL,RCODE_ML,RCODE_GCL,RCODE_WM=3,0,1,2

def region_of(masks, ix, iy):
    if masks['PCL'][iy,ix]: return 3
    if masks['ML'][iy,ix]:  return 0
    if masks['WM'][iy,ix]:  return 2
    if masks['GCL'][iy,ix]: return 1
    return -1

sections={}; cards=[]
for sec in KEEP:
    sd=LS+f"/sections/{sec}"
    meta=json.load(open(sd+"/grid_meta.json"))
    masks={k:np.load(sd+f"/{k}_mask.npy") for k in ['ML','PCL','GCL','WM']}
    g=allcb[allcb['brain_section_label']==sec].copy()
    sc=g['subclass'].astype(str); g['is_purk']=sc.str.contains('Purkinje')
    x=g['x'].values; y=g['y'].values
    ix=np.clip(((x-meta['x0'])/meta['res']).astype(int),0,meta['nx']-1)
    iy=np.clip(((y-meta['y0'])/meta['res']).astype(int),0,meta['ny']-1)
    reg=np.array([region_of(masks,ix[k],iy[k]) for k in range(len(g))])
    # BG / candidate sets from assignments
    seg=A[A['section']==sec]
    bg_ids=set(seg.index[seg['region']=='BG'])
    cand=seg[seg['region']=='FC'] if CAND is None else CAND[CAND['section']==sec]
    typ=np.zeros(len(g),int)
    idset=set(g.index)
    isbg=g.index.isin(bg_ids)
    typ[(reg==3)&g['is_purk'].values]=1     # Purkinje on PCL
    typ[isbg]=2                           # BG
    iscand=g.index.isin(set(cand.index))
    # downsample for display: keep PCL-region, Purkinje, BG, candidates + sample of rest
    keep=(reg==3)|(typ>0)|iscand
    rest=np.where(~keep)[0]; cap=9000; nrest=max(0,cap-int(keep.sum()))
    rng=np.random.default_rng(0)
    if len(rest)>nrest: rest=rng.choice(rest,nrest,replace=False)
    sel=np.unique(np.concatenate([np.where(keep)[0],rest]))
    cells=[[round(float(x[k]),3),round(float(y[k]),3),int(reg[k]) if reg[k]>=0 else 1,int(typ[k])] for k in sel]
    sections[sec]=dict(cells=cells,xmin=round(float(x.min()),3),xmax=round(float(x.max()),3),
                       ymin=round(float(y.min()),3),ymax=round(float(y.max()),3))
    an=sec.split('.')[0]
    for cid,r in cand.iterrows():
        cards.append(dict(id=str(cid),dataset="Zhuang",section=sec,animal=an,
            cx=round(float(r['x']),3),cy=round(float(r['y']),3),
            conf=round(float(r['confidence']),2), dist_pcl=round(float(r['dist_to_PCL_um']),0),
            support=int(r['purk_support'])))
cards.sort(key=lambda c:-c['conf'])   # highest-confidence first
data=dict(meta=dict(dataset=DSLABEL,zoom_um=250,n=len(cards),
                    source=SOURCE),
          sections=sections, cards=cards)
open(OUT+f"/swipe_{TAG}.js","w",encoding="utf-8").write("window.SWIPE_DATA="+json.dumps(data)+";")
open(OUT+f"/swipe_{TAG}.json","w",encoding="utf-8").write(json.dumps(data))
print("sections:",len(sections)," candidates(FC):",len(cards))
print("js MB:",round(os.path.getsize(OUT+f'/swipe_{TAG}.js')/1e6,2))
