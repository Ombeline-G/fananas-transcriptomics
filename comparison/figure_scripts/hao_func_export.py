#!/usr/bin/env python
"""Export Mouse1-only Hao functional expression and whole-transcriptome group summaries."""
import os,re,sys,json,numpy as np,pandas as pd,anndata as ad,scipy.sparse as sp
sys.path.insert(0,os.path.dirname(__file__)); from func_panels import ALL_FUNC
from project_paths import HAO as HAO_PATH, TABLES
HAO=str(HAO_PATH); TAB=str(TABLES)
rng=np.random.default_rng(0)
def strip(s): return re.sub(r'-\d+$','',str(s))
dec=pd.read_csv(HAO+"/swipe_validation/fc_decisions_hao.csv")
dec=dec[dec.animal.eq('Mouse1')].copy()
ml=dec[dec.decision=='yes'].copy(); ml['cl']=ml.cell_label.astype(str).map(strip)
seg=pd.read_csv(HAO+"/layer_segmentation/BG_assignments_all.csv")
seg=seg[seg.section.astype(str).str.startswith('Mouse1_')].copy()
seg['cl']=seg.cell_label.astype(str).map(strip)
pclseg=seg[seg.region=='BG']
def secpath(sec):
    m=sec.split('_')[0]; return f"{HAO}/data_raw/{'stereo_seq_mouse1' if m=='Mouse1' else 'stereo_seq_mouse2'}/{sec}.h5ad"
secs=sorted(set(ml.section)|set(pclseg.section)); frames=[]; metas=[]
stat_sum={g:pd.Series(dtype=float) for g in ('FC','BG')}
stat_nz={g:pd.Series(dtype=float) for g in ('FC','BG')}
stat_n={g:0 for g in ('FC','BG')}
for k,sec in enumerate(secs):
    p=secpath(sec)
    if not os.path.exists(p): continue
    a=ad.read_h5ad(p); a.var_names_make_unique(); idx=a.obs.index.astype(str); cl=np.array([strip(x) for x in idx])
    tot=np.asarray(a.X.sum(1)).ravel().astype(float); tot[tot==0]=1
    fg=[g for g in ALL_FUNC if g in set(map(str,a.var_names))]
    Xp=a[:,fg].X; Xp=Xp.toarray() if sp.issparse(Xp) else np.asarray(Xp); norm=np.log1p(Xp/tot[:,None]*1e4)/np.log(2)  # log2(1+CP10k)
    mlset=set(ml[ml.section==sec]['cl']); pclset=set(pclseg[pclseg.section==sec]['cl'])
    is_ml=np.isin(cl,list(mlset)); is_pcl=np.isin(cl,list(pclset))&~is_ml
    pcl_pos=np.where(is_pcl)[0]
    if len(pcl_pos)>50: pcl_pos=rng.choice(pcl_pos,50,replace=False)
    sel=np.concatenate([np.where(is_ml)[0],pcl_pos]).astype(int)
    if not len(sel): continue
    grp=np.where(np.isin(sel,np.where(is_ml)[0]),'FC','BG')
    cid=[f"{sec}|{idx[pp]}" for pp in sel]
    frames.append(pd.DataFrame(norm[sel],columns=fg,index=cid)); metas.append(pd.DataFrame({'cell':cid,'group':grp,'animal':sec.split('_')[0]}))
    # Whole-transcriptome means/fractions are accumulated without retaining a
    # large dense cell-by-gene matrix. Missing genes in a section count as zero.
    Xall=a.X[sel]
    if sp.issparse(Xall):
        Xall=Xall.multiply((1e4/tot[sel])[:,None]).tocsr()
        Xall.data=np.log1p(Xall.data)/np.log(2)  # log2(1+CP10k)
    else:
        Xall=np.log1p(np.asarray(Xall)/tot[sel,None]*1e4)/np.log(2)  # log2(1+CP10k)
    genes=pd.Index(a.var_names.astype(str))
    for group in ('FC','BG'):
        keep=np.where(grp==group)[0]
        if not len(keep): continue
        Xg=Xall[keep]
        sm=np.asarray(Xg.sum(axis=0)).ravel()
        nz=np.asarray((Xg>0).sum(axis=0)).ravel()
        stat_sum[group]=stat_sum[group].add(pd.Series(sm,index=genes),fill_value=0)
        stat_nz[group]=stat_nz[group].add(pd.Series(nz,index=genes),fill_value=0)
        stat_n[group]+=len(keep)
E=pd.concat(frames).fillna(0.0); meta=pd.concat(metas,ignore_index=True)
E=E.loc[~E.index.duplicated()]; meta=meta.drop_duplicates('cell').set_index('cell').reindex(E.index).reset_index()
out=pd.concat([meta.reset_index(drop=True),E.reset_index(drop=True)],axis=1)
out.to_parquet(TAB+"/hao_func_expr.parquet")
stats=[]
all_genes=stat_sum['FC'].index.union(stat_sum['BG'].index)
for group in ('FC','BG'):
    sums=stat_sum[group].reindex(all_genes,fill_value=0)
    nz=stat_nz[group].reindex(all_genes,fill_value=0)
    stats.append(pd.DataFrame({'dataset':'Stereo-seq','group':group,'gene':all_genes,
                               'mean':sums.values/stat_n[group],'frac':nz.values/stat_n[group]}))
pd.concat(stats,ignore_index=True).to_csv(TAB+"/hao_bg_gene_stats.csv",index=False)
print("Hao func expr:",E.shape,"groups",meta.group.value_counts().to_dict(),"genes",len([g for g in ALL_FUNC if g in E.columns]),"/",len(ALL_FUNC))
print("Hao whole-transcriptome summaries:",len(all_genes),"genes; n",stat_n)
