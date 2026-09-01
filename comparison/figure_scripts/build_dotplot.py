#!/usr/bin/env python
"""Combine FC/BG gene stats (Zhuang 1122 / Hao whole-tx / Suyama 5006) -> top most-expressed genes
common to all 3, scaled per dataset -> fig3_dotplot_data.csv (dataset,group,gene,mean_scaled,frac)."""
import numpy as np,pandas as pd
from project_paths import TABLES, SUYAMA
TAB=str(TABLES)
SU=str(SUYAMA)
def stats_from_expr(cells_pq,expr_pq):
    m=pd.read_parquet(cells_pq)
    if 'cell' not in m.columns: m=m.rename(columns={'index':'cell'})
    E=pd.read_parquet(expr_pq).reindex(m['cell'].values); rows=[]
    for grp in ['FC','BG']:
        sub=E[(m.group==grp).values]
        for g in E.columns: v=sub[g].values; rows.append(dict(group=grp,gene=g,mean=float(v.mean()),frac=float((v>0).mean())))
    return pd.DataFrame(rows)
# Zhuang (1122)
z=stats_from_expr(f"{TAB}/zhuang_astro_cells.parquet",f"{TAB}/zhuang_astro_expr.parquet"); z['dataset']='Zhuang'
# Suyama (5006 full, FC+BG)
e=pd.read_csv(SU+"/results/tables/bergmann_full_expr.csv"); rows=[]
for grp in ['FC','BG']:
    sub=e[e.group==grp].drop(columns=['cell','group'])
    for g in sub.columns: v=sub[g].values; rows.append(dict(dataset='Suyama',group=grp,gene=g,mean=float(v.mean()),frac=float((v>0).mean())))
s=pd.DataFrame(rows)
# Hao (whole-tx)
h=pd.read_csv(f"{TAB}/hao_bg_gene_stats.csv")
D=pd.concat([z,h,s],ignore_index=True)
# glob_scaled = mean / dataset's whole-panel max -> used ONLY to pick "most expressed" genes
D['glob_scaled']=D.groupby('dataset')['mean'].transform(lambda x:x/(x.max()+1e-9))
import re  # drop technical + ubiquitous housekeeping (mito/ribo/cytoskeleton/metallothionein)
TECH=re.compile(r'^(mt-|Rpl|Rps|Gm\d|Rn\d|Bc1$|Malat1$|Neat1$|Cox\d|Nduf|Uqcr|Atp5|Tmsb|Tub[ab]|Act[bg]|Fth1$|Ftl1$|Mt[12]$|Fau$|Eef|Eif|Hsp|Ubb$|Ubc$|Cox6|Cox7|Cox8|Ppia$|Gapdh$)',re.I)
# Exclude ONLY Pvalb (MLI marker) — it was the single most-expressed gene in the Hao Bergmann bins and
# dominated/flattened the colour scale. Other genes (incl. contamination markers) are kept.
CONTAM={'Pvalb'}
D=D[~D.gene.astype(str).str.contains(TECH)]
D=D[~D.gene.isin(CONTAM)]
per=D.groupby(['dataset','gene'])['glob_scaled'].mean().reset_index()
def rescale_within_displayed(sub):
    # colour scale: normalise each dataset by the max among the DISPLAYED genes (fairer; Hao not washed-out)
    sub=sub.copy(); sub['mean_scaled']=sub.groupby('dataset')['mean'].transform(lambda x:x/(x.max()+1e-9)); return sub
# (1) COMMON: genes present in all 3, top-20 most-expressed (dense, directly comparable)
common=set(z.gene)&set(h.gene)&set(s.gene)
keepc=D[D.gene.isin(common)].groupby('gene')['glob_scaled'].mean().sort_values(ascending=False).head(20).index.tolist()
rescale_within_displayed(D[D.gene.isin(keepc)]).to_csv(f"{TAB}/fig3_dotplot_common.csv",index=False)
print("COMMON genes (all 3):",len(common),"-> top",len(keepc))
# (2) UNION of each dataset's top-expressed genes (blanks where a gene isn't in a panel)
topset=set()
for ds in D.dataset.unique(): topset|=set(per[per.dataset==ds].sort_values('glob_scaled',ascending=False).head(12).gene)
keepu=D[D.gene.isin(topset)].groupby('gene')['glob_scaled'].mean().sort_values(ascending=False).head(26).index.tolist()
rescale_within_displayed(D[D.gene.isin(keepu)]).to_csv(f"{TAB}/fig3_dotplot_union.csv",index=False)
print("UNION genes:",len(keepu),"| colour re-scaled within displayed genes per dataset")
