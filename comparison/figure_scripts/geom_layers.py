#!/usr/bin/env python
"""Shared cerebellar-layer assignment used by Figure 2 & Figure 3 (Suyama AND Zhuang).

The layers are three strictly PARALLEL bands built from a signed-distance field, so the boundary is
a clean line, not a fuzzy per-cell mixture:

  1. Build the granule MASS as a solid 2-D region: bin the section, local granule fraction,
     Gaussian-smooth, threshold, morphological close/open/fill, drop small islands.
  2. Signed distance to that region's boundary (EDT):  s > 0 inside the mass,  s < 0 outside.
  3. Assign by s, with a PCL BUFFER of half-thickness W_PCL that neither side may cross:
         GCL : s >  W_PCL                    (deep to the Purkinje line)
         PCL : |s| <= W_PCL                  (the Purkinje line, a band of fixed thickness)
         ML  : -D_ML < s < -W_PCL            (superficial to the Purkinje line)
     WM (oligodendrocyte-rich) overrides outside the mass; cells far from the boundary = none.
  Because assignment is monotone in s, NO GCL cell can sit on the ML side and NO ML cell on the GCL
  side — they are separated by the >=2*W_PCL-thick PCL band.

Granule cells are EXCLUDED from the ML and PCL bands (dropped to 'none'): the molecular layer and the
Purkinje line contain 0 granule cells by construction. Actual Purkinje somata are forced to PCL.
"""
import numpy as np
from scipy.ndimage import (gaussian_filter, binary_closing, binary_opening,
                           binary_fill_holes, label as ndlabel, distance_transform_edt,
                           binary_dilation)

def assign_layers(x, y, is_gran, is_purk, is_oli, umpp=1.0,
                  BIN=15.0, D_ML=175.0, W_PCL=12.0, gran_thr=0.35, oli_thr=0.30,
                  min_region_um2=13000.0, smooth_um=30.0,
                  force_purkinje=True, exclude_granule=True, wm_grow_um=0.0):
    # force_purkinje / exclude_granule (default True) apply per-cell-type overrides on top of the
    # geometry. Set BOTH False for a PURELY positional assignment (a cell's layer = its signed-distance
    # band only, no use of its own cell type) — used by Figure 1E. Defaults preserve prior behaviour.
    # wm_grow_um (default 0) dilates the WM (oligo) mask by that many µm INTO non-granule tissue only,
    # so WM claims the deep GCL<->WM transition zone that signed-distance-alone mislabels as ML
    # (both the superficial ML side and the deep WM side have s<0). 0 = unchanged.
    X=np.asarray(x,float)*umpp; Y=np.asarray(y,float)*umpp
    ig=np.asarray(is_gran,bool); ip=np.asarray(is_purk,bool); io=np.asarray(is_oli,bool)
    ix=((X-X.min())/BIN).astype(int); iy=((Y-Y.min())/BIN).astype(int)
    nx,ny=ix.max()+1,iy.max()+1
    def grid(v):
        G=np.zeros((nx,ny)); np.add.at(G,(ix,iy),v); return G
    sig=smooth_um/BIN
    tot=gaussian_filter(grid(np.ones(len(X))),sig)
    fracG=gaussian_filter(grid(ig.astype(float)),sig)/(tot+1e-6)
    fracO=gaussian_filter(grid(io.astype(float)),sig)/(tot+1e-6)
    gmask=fracG>gran_thr
    gmask=binary_closing(gmask,iterations=2); gmask=binary_opening(gmask,iterations=2)
    gmask=binary_fill_holes(gmask)
    lab,_=ndlabel(gmask); lsz=np.bincount(lab.ravel()); lsz[0]=0
    minbins=max(1,int(min_region_um2/(BIN*BIN)))
    gmask=np.isin(lab,np.where(lsz>=minbins)[0])
    # signed distance to the mass boundary, in µm (+inside / -outside)
    sdist=((distance_transform_edt(gmask)-distance_transform_edt(~gmask))*BIN)[ix,iy]
    wmask=fracO>oli_thr
    if wm_grow_um>0:                                 # grow WM into adjacent NON-granule tissue (deep transition zone)
        wmask=binary_dilation(wmask,iterations=max(1,int(round(wm_grow_um/BIN))))&(~gmask)
    wm_bin=wmask[ix,iy]
    lay=np.full(len(X),'none',dtype=object)
    lay[(sdist>W_PCL)&(~wm_bin)]='GCL'
    lay[(sdist<-W_PCL)&(sdist>-D_ML)&(~wm_bin)]='ML'
    lay[np.abs(sdist)<=W_PCL]='PCL'                 # fixed-thickness Purkinje band (separator)
    lay[wm_bin&(sdist<=W_PCL)]='WM'                 # white matter outside the mass
    if exclude_granule:
        lay[ig&np.isin(lay,['ML','PCL'])]='none'    # 0 granule cells in ML / PCL by construction
    if force_purkinje:
        lay[ip]='PCL'                                # Purkinje somata define the line
    return lay, sdist


def layer_regions(x, y, is_gran, is_oli, umpp=1.0, BIN=10.0, D_ML=175.0, W_PCL=12.0,
                  gran_thr=0.35, oli_thr=0.30, min_region_um2=13000.0, smooth_um=28.0,
                  tissue_close_um=45.0):
    """Grid fields for PAINTING the three layers as solid regions (not scattered cells).
    Returns bin-centre coords (in ORIGINAL units) and, per grid bin, the signed distance to the
    granule-mass boundary, a tissue mask, a white-matter mask, and distance-to-pia — everything
    needed to fill GCL | PCL band | ML as parallel regions and draw the Purkinje line at s=0."""
    X=np.asarray(x,float)*umpp; Y=np.asarray(y,float)*umpp
    ix=((X-X.min())/BIN).astype(int); iy=((Y-Y.min())/BIN).astype(int)
    nx,ny=ix.max()+1,iy.max()+1
    def grid(v):
        G=np.zeros((nx,ny)); np.add.at(G,(ix,iy),v); return G
    sig=smooth_um/BIN
    tot=grid(np.ones(len(X))); totS=gaussian_filter(tot,sig)
    fracG=gaussian_filter(grid(np.asarray(is_gran,float)),sig)/(totS+1e-6)
    fracO=gaussian_filter(grid(np.asarray(is_oli,float)),sig)/(totS+1e-6)
    gmask=fracG>gran_thr
    gmask=binary_closing(gmask,iterations=2); gmask=binary_opening(gmask,iterations=2); gmask=binary_fill_holes(gmask)
    lab,_=ndlabel(gmask); lsz=np.bincount(lab.ravel()); lsz[0]=0
    minbins=max(1,int(min_region_um2/(BIN*BIN))); gmask=np.isin(lab,np.where(lsz>=minbins)[0])
    sdist=(distance_transform_edt(gmask)-distance_transform_edt(~gmask))*BIN
    # tissue footprint (largest connected component), gaps closed
    cl=int(round(tissue_close_um/BIN))
    tissue=binary_closing(tot>0,iterations=cl); tissue=binary_fill_holes(tissue)
    labT,_=ndlabel(tissue); tsz=np.bincount(labT.ravel()); tsz[0]=0
    tissue=labT==(tsz.argmax() if tsz.size>1 else 0)
    d_ext=distance_transform_edt(tissue)*BIN             # distance to pia / tissue exterior
    wmask=fracO>oli_thr
    Xc=(X.min()+(np.arange(nx)+0.5)*BIN)/umpp; Yc=(Y.min()+(np.arange(ny)+0.5)*BIN)/umpp
    XX,YY=np.meshgrid(Xc,Yc,indexing='ij')
    return dict(XX=XX,YY=YY,sdist=sdist,tissue=tissue,wmask=wmask,d_ext=d_ext,gmask=gmask,D_ML=D_ML,W_PCL=W_PCL,BIN=BIN)


def draw_pcl_line(ax, reg, color='#b8560f', lw=1.4, zorder=6):
    """Draw ONLY the Purkinje line (granule-mass boundary, s=0) as a continuous curve — no fills."""
    import numpy as _np
    ax.contour(reg['XX'],reg['YY'],_np.where(reg['tissue'],reg['sdist'],_np.nan),
               levels=[0.0],colors=[color],linewidths=lw,zorder=zorder)


def paint_layers(ax, reg, colors=None, pcl_line='#b8560f', pcl_lw=1.8, alpha=1.0,
                 wm_min_um2=40000.0):
    """Fill GCL | PCL band | ML | WM as solid regions on ax and draw the Purkinje line at s=0.
    Every non-WM tissue bin is partitioned by the signed distance s to the granule-mass boundary:
    GCL (s>=W, inside) | PCL band (|s|<W) | ML (s<=-W, the whole superficial side incl. the
    inter-folia molecular layer). WM keeps only sizeable components adjacent to the cortex."""
    import numpy as _np
    from scipy.ndimage import binary_dilation as _dil
    C=dict(ML='#8fc6e6',PCL='#e8894a',GCL='#c39bd8',WM='#9aa0a6')
    if colors: C.update(colors)
    sd=reg['sdist']; T=reg['tissue']; WM=reg['wmask']; XX=reg['XX']; YY=reg['YY']
    W=reg['W_PCL']; BIN=reg['BIN']
    gnear=_dil(reg['gmask'],iterations=max(1,int(80.0/BIN)))   # within 80 µm of the granule mass
    wm=T&WM; labw,_=ndlabel(wm); minb=max(1,int(wm_min_um2/(BIN*BIN))); keep=[]
    for c in range(1,int(labw.max())+1):
        comp=labw==c
        if comp.sum()>=minb and (comp&gnear).any(): keep.append(c)   # sizeable AND touching the cortex
    wm=_np.isin(labw,keep)
    nz=T&(~wm)
    band=_np.where(nz,_np.where(sd>=W,2.0,_np.where(sd<=-W,-2.0,0.0)),_np.nan)   # GCL / ML / PCL band
    ax.contourf(XX,YY,band,levels=[-3,-1,1,3],colors=[C['ML'],C['PCL'],C['GCL']],alpha=alpha,antialiased=True)
    ax.contourf(XX,YY,_np.where(wm,1.0,_np.nan),levels=[.5,1.5],colors=[C['WM']],alpha=alpha,antialiased=True)
    ax.contour(XX,YY,_np.where(T,sd,_np.nan),levels=[0.0],colors=[pcl_line],linewidths=pcl_lw)   # Purkinje line

