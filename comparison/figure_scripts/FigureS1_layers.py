#!/usr/bin/env python
"""Figure S1 - Definition and marker validation of cerebellar layers per dataset.
Same sections as Figure 1: Zhuang ABCA-3.005, Hao Mouse1_T176, Suyama Adult.
Row 1: per-cell layer map (ML/PCL/GCL/WM).
Row 2: layer x marker-panel enrichment matrix (mean z-scored panel score per assigned layer;
       hot diagonal = layers correctly defined). Columns: PCL/ML/GCL/WM marker panels + Astroglia."""
import os, sys, json, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec; import matplotlib.patches as mp
from matplotlib.lines import Line2D
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from geom_layers import assign_layers, layer_regions, draw_pcl_line
from figstyle import use_style, LAYER, panel_letter, save
use_style()

from project_paths import ROOT
BASE=str(ROOT)
ZH=BASE+"/datasets/zhang_zhuang_2023_merfish"; HAO=BASE+"/datasets/hao_2024_stereo_seq"; SU=BASE+"/datasets/suyama_2026_xenium"

LAYERS=['ML','PCL','GCL','WM']
LCOL={k:LAYER[k] for k in ['ML','PCL','GCL','WM']}
PANELS=['PCL','ML','GCL','WM','Astroglia']
# Layer-specific marker panels. ML uses molecular-layer-interneuron markers
# (Kit/Lypd6/Sorcs3/Nxph1/...) instead of broad pan-GABAergic signal.
MK={'PCL':['Pcp2','Calb1','Car8','Itpr1','Grid2','Aldoc','Prkcg','Homer3'],
    'ML':['Kit','Lypd6','Sorcs3','Nxph1','Cxcl14','Grm8','Prkcd','Nxph2'],
    'GCL':['Gabra6','Neurod1','Barhl1','Zic1','Zic2','Rbfox3','Cadps2','Etv1'],
    'WM':['Mbp','Plp1','Mog','Mag','Cnp'],
    'Astroglia':['Aldh1l1','Sox9','Slc1a3','Slc1a2','S100b','Gja1','Gjb6','Glul','Aqp4','Fabp7','Hopx','Gpr37l1','Ednrb','Notch2']}

def enrich(layer_lab, panel_score):
    # raw mean panel score per assigned layer, then z-score each panel (column) ACROSS the 4 layers
    # so every diagonal is on the same scale.
    lab=np.asarray(layer_lab)
    raw=np.full((4,5),np.nan)
    for i,L in enumerate(LAYERS):
        m=(lab==L)
        if m.sum():
            for j,p in enumerate(PANELS): raw[i,j]=np.nanmean(np.asarray(panel_score[p],float)[m])
    M=np.full((4,5),np.nan)
    for j in range(5):
        col=raw[:,j]; mu=np.nanmean(col); sd=np.nanstd(col)
        M[:,j]=(col-mu)/sd if (sd and not np.isnan(sd) and sd>0) else 0.0
    return M

# ================= ZHUANG (granule-mass signed-distance layers) =================
zsec='Zhuang-ABCA-3.005'
zall=pd.read_parquet(ZH+"/cb_cells_multivote.parquet")
zz=zall[zall.brain_section_label==zsec].copy()
zlab,_=assign_layers(zz.x,zz.y,zz.is_gran,zz.is_purk,
                     zz.subclass.astype(str).str.contains('Oligo'),umpp=1000.0,
                     BIN=22.0,smooth_um=45.0,D_ML=280.0)   # MERFISH: sparser -> smoother mask, wider ML reach
za=ad.read_h5ad(ZH+"/abc_cache/expression_matrices/Zhuang-ABCA-3/20230830/Zhuang-ABCA-3-log2.h5ad",backed='r')
sym2ens=dict(zip(za.var['gene_symbol'].astype(str), za.var_names))
zgenes={p:[g for g in MK[p] if g in sym2ens] for p in PANELS}
allg=sorted({g for p in PANELS for g in zgenes[p]})
ens=[sym2ens[g] for g in allg]
obset=set(za.obs_names)
ids=[c for c in zz.index if c in obset]
sub=za[ids].to_memory()[:,ens]   # rows first (single fancy index on backed h5ad), then columns in memory
X=sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
expr=pd.DataFrame(X,index=list(sub.obs_names),columns=allg).reindex(zz.index)
zscore={p:(expr[zgenes[p]].mean(1).values if zgenes[p] else np.zeros(len(zz))) for p in PANELS}
zxy=(zz.x.values, zz.y.values)
zM=enrich(zlab, zscore)                                   # marker validation: PCL = Purkinje-defined (clean)
zreg=layer_regions(zz.x,zz.y,zz.is_gran,zz.subclass.astype(str).str.contains('Oligo'),umpp=1000.0)
# MAP display: the Purkinje cell layer is Purkinje + Bergmann somata -> show Bergmann as PCL (visible line,
# esp. Zhuang where Purkinje are sparse). Kept OUT of the enrichment above to avoid diluting the PCL diagonal.
zlab_map=np.array(zlab,dtype=object); zlab_map[zz.subclass.astype(str).str.contains('Bergmann').values]='PCL'
print("Zhuang layer n:",pd.Series(zlab).value_counts().to_dict())
print("Zhuang genes used:",{p:zgenes[p] for p in PANELS})

# ================= HAO =================
hsec='Mouse1_T176'
ha=ad.read_h5ad(f"{HAO}/data_raw/stereo_seq_mouse1/{hsec}.h5ad")
hmap={'molecular layer':'ML','granular layer':'GCL','purkinje layer':'PCL','white matter':'WM'}
hlab=pd.Series(ha.obs['annotation'].astype(str).values).map(hmap).values
hvar=set(map(str,ha.var_names))
hgenes={p:[g for g in MK[p] if g in hvar] for p in PANELS}
allh=sorted({g for p in PANELS for g in hgenes[p]})
Xh=ha[:,allh].X; Xh=Xh.toarray() if sp.issparse(Xh) else np.asarray(Xh)
tot=np.asarray(ha.X.sum(1)).ravel(); tot[tot==0]=1
Xh=np.log1p(Xh/tot[:,None]*np.median(tot))
he=pd.DataFrame(Xh,columns=allh)
hscore={p:(he[hgenes[p]].mean(1).values if hgenes[p] else np.zeros(Xh.shape[0])) for p in PANELS}
hxy=(ha.obs['coor_x'].values.astype(float), ha.obs['coor_y'].values.astype(float))
hM=enrich(hlab, hscore)
hreg=layer_regions(ha.obs['coor_x'].astype(float),ha.obs['coor_y'].astype(float),
                   pd.Series(hlab).eq('GCL').values, pd.Series(hlab).eq('WM').values, umpp=10.0)
print("Hao genes used:",{p:hgenes[p] for p in PANELS})

# ================= SUYAMA (granule-mass signed-distance layers; see geom_layers.py) =================
sgl=pd.read_csv(SU+"/results/tables/suyama_geom_layers.csv"); sgl['cell']=sgl['cell'].astype(str)
sx=pd.read_csv(SU+"/results/tables/adult_marker_expr.csv"); sx['cell']=sx['cell'].astype(str)  # REAL normalized Xenium genes
sexpr=sx.set_index('cell').reindex(sgl['cell'].values)
svar=set(sexpr.columns)
sgenes={p:[g for g in MK[p] if g in svar] for p in PANELS}
sscore={p:(sexpr[sgenes[p]].mean(1).values if sgenes[p] else np.zeros(len(sgl))) for p in PANELS}
slab=sgl['layer_geom'].values.astype(object)
keep=slab!='none'   # drop the non-cortical DCN blob / sparse periphery
slab=slab[keep]
for p in PANELS: sscore[p]=sscore[p][keep]
sxy=(sgl['x'].values[keep], sgl['y'].values[keep])
sM=enrich(slab, sscore)                                   # marker validation: PCL = Purkinje-defined (clean)
slab_map=slab.copy(); slab_map[sgl['ct'].values[keep]=='Bergmann']='PCL'   # display: Bergmann somata shown as PCL
sreg=layer_regions(sgl['x'],sgl['y'],sgl['ct']=='Granule',sgl['ct']=='Oligodendrocyte',umpp=1.0)
print("Suyama genes used:",{p:sgenes[p] for p in PANELS})
print("Suyama layer n:",pd.Series(slab).value_counts().to_dict())

# ================= FIGURE =================
DS=[('A. MERFISH (Zhang et al., 2023)',zsec,zxy,zlab_map,zreg,zM,1000.0),
    ('B. Stereo-seq (Hao et al., 2024)',hsec,hxy,hlab,hreg,hM,10.0),
    ('C. Xenium 5K (Suyama et al., 2026)','Adult cerebellum',sxy,slab_map,sreg,sM,1.0)]
LFILL={'ML':'#8fc6e6','PCL':'#e8894a','GCL':'#c39bd8','WM':'#9aa0a6'}

fig=plt.figure(figsize=(15,11.5))
gs=GridSpec(2,3,figure=fig,height_ratios=[1.15,1.0],hspace=0.30,wspace=0.12)

def scalebar(ax,umpp,length_um=500):
    x0,x1=ax.get_xlim(); y0,y1=sorted(ax.get_ylim()); L=length_um/umpp
    xs=x0+(x1-x0)*0.06; ys=y0+(y1-y0)*0.06
    ax.add_line(Line2D([xs,xs+L],[ys,ys],color='k',lw=3,solid_capstyle='butt'))
    ax.text(xs+L/2,ys+(y1-y0)*0.02,f'{length_um} um',ha='center',va='bottom',fontsize=8)

for c,(title,sub,xy,lab,reg,M,umpp) in enumerate(DS):
    ax=fig.add_subplot(gs[0,c]); x,y=np.asarray(xy[0]),np.asarray(xy[1]); lab=np.asarray(lab)
    for L in ['GCL','WM','ML']:                                    # cell-by-cell (points), not zones
        m=(lab==L)
        if m.sum(): ax.scatter(x[m],y[m],s=3,c=LFILL[L],linewidths=0,rasterized=True,zorder=1)
    mp_=(lab=='PCL')
    if mp_.sum(): ax.scatter(x[mp_],y[mp_],s=5,c=LFILL['PCL'],linewidths=0,rasterized=True,zorder=4)
    ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    for s_ in ax.spines.values(): s_.set_visible(False)
    ax.set_title(title,fontsize=12,loc='left')
    panel_letter(ax, "abc"[c], dx=-0.04, dy=1.08)
    ax.text(0.5,-0.04,sub,transform=ax.transAxes,ha='center',fontsize=9,color='#555')
    scalebar(ax,umpp)
    axm=fig.add_subplot(gs[1,c])
    im=axm.imshow(M,cmap='RdBu_r',vmin=-1.6,vmax=1.6,aspect='auto')
    axm.set_xticks(range(5)); axm.set_xticklabels(PANELS,fontsize=8.5,rotation=45,ha='right')
    axm.set_yticks(range(4)); axm.set_yticklabels(LAYERS,fontsize=9)
    axm.set_xlabel('Marker panel',fontsize=9)
    if c==0: axm.set_ylabel('Assigned layer',fontsize=9)
    for i in range(4):
        for j in range(5):
            if not np.isnan(M[i,j]):
                axm.text(j,i,f'{M[i,j]:.1f}',ha='center',va='center',fontsize=7.5,
                         color='white' if abs(M[i,j])>0.7 else '#333')
    for i,L in enumerate(LAYERS):
        j=PANELS.index(L); axm.add_patch(plt.Rectangle((j-0.5,i-0.5),1,1,fill=False,ec='k',lw=2))

handles=[mp.Patch(color=LFILL[L],label=('PCL (Purkinje layer)' if L=='PCL' else
                  {'ML':'ML (molecular)','GCL':'GCL (granular)','WM':'WM (white matter)'}[L])) for L in LAYERS]
fig.legend(handles=handles,loc='upper center',ncol=4,fontsize=10,frameon=False,bbox_to_anchor=(0.5,0.945))
cax=fig.add_axes([0.92,0.12,0.014,0.28]); cb=fig.colorbar(im,cax=cax)
cb.set_label('panel score (z across layers)',fontsize=9)
# NOTE (not drawn on the figure): "assigned layer" = the per-cell layer label (ML/PCL/GCL/WM) that
# is the study's spatial reference — it lets us tell BG (Bergmann at the Purkinje layer) from
# FC (Bergmann displaced into the molecular layer). The bottom matrix checks each layer expresses
# its own markers (hot boxed diagonal).
save(fig, "FigureS1_layers")
plt.close(fig)

# ---- export enrichment matrices (reproducibility / supplementary) ----
rows=[]
for name,M in [('Zhuang2023',zM),('Hao2024',hM),('Suyama2026',sM)]:
    for i,L in enumerate(LAYERS):
        rows.append({'dataset':name,'assigned_layer':L,
                     **{f'panel_{PANELS[j]}':round(float(M[i,j]),3) for j in range(5)}})
pd.DataFrame(rows).to_csv(BASE+"/comparison/tables/FigureS1_layer_enrichment.csv",index=False)
with open(BASE+"/comparison/tables/FigureS1_markers_used.json","w") as f:
    json.dump({'Zhuang2023':zgenes,'Hao2024':hgenes,'Suyama2026':sgenes},f,indent=2)
print("saved enrichment csv + markers json")
