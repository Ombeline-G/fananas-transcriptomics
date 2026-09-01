#!/usr/bin/env python
"""Build the Kozareva 2021 mouse atlas validation table and QC figure.

The normalized h5ad is the Kozareva atlas redistributed with Hao 2024; it is
not an independent Hao snRNA-seq dataset.
"""
import os,numpy as np,pandas as pd,anndata as ad,scanpy as sc,scipy.sparse as sp
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from project_paths import HAO, FIGURES, TABLES
F=str(HAO / "data_raw" / "snRNAseq" / "Mouse.sn.norm.h5ad")
FIG=str(FIGURES); TAB=str(TABLES)
SIG={'Bergmann_shared':['Hopx','Gpr37l1','A2m','Gdf10','Slc1a3','Timp4','Fabp7','Tnc','Ptprz1'],
     'Parenchymal_astro':['Aqp4','Gja1','Slc1a2','Glul','S100b','Aldh1l1'],
     'FC_astro_lean':['Gfap','Agt','Ntsr2','Aqp4']}
a=ad.read_h5ad(F,backed='r')
mask=a.obs['annotation'].isin(['Bergmann','Astrocyte']).values
A=a[mask].to_memory(); print("subset:",A.shape, A.obs['annotation'].value_counts().to_dict())
ct=A.obs['annotation'].astype(str).values
def score(M,genes):
    g=[x for x in genes if x in set(map(str,M.var_names))]
    X=M[:,g].X; X=X.toarray() if sp.issparse(X) else np.asarray(X)
    Z=(X-X.mean(0))/(X.std(0)+1e-9); return Z.mean(1),g
for k,gs in SIG.items(): A.obs['S_'+k],_=score(A,gs)
# UMAP on the subset
sc.pp.highly_variable_genes(A,n_top_genes=2000); Ah=A[:,A.var.highly_variable].copy()
sc.pp.scale(Ah,max_value=10); sc.tl.pca(Ah,n_comps=30); sc.pp.neighbors(Ah,n_neighbors=15,n_pcs=30); sc.tl.umap(Ah,random_state=0)
U=Ah.obsm['X_umap']
# recluster Bergmann only
bgmask=(ct=='Bergmann'); Ab=A[bgmask].copy()
sc.pp.highly_variable_genes(Ab,n_top_genes=2000); Abh=Ab[:,Ab.var.highly_variable].copy()
sc.pp.scale(Abh,max_value=10); sc.tl.pca(Abh,n_comps=30); sc.pp.neighbors(Abh,n_neighbors=15,n_pcs=30)
sc.tl.leiden(Abh,resolution=0.5,random_state=0,flavor="igraph",n_iterations=2,directed=False); sc.tl.umap(Abh,random_state=0)
Ab.obs['bsub']=Abh.obs['leiden'].values; bgU=Abh.obsm['X_umap']
for k in SIG: Ab.obs['S_'+k]=A.obs['S_'+k].values[bgmask]
sub_prof=Ab.obs.groupby('bsub')[['S_'+k for k in SIG]].mean(); sub_prof.to_csv(TAB+"/atlas_kozareva_bergmann_subclusters.csv")
print("Bergmann subclusters:",Ab.obs['bsub'].value_counts().to_dict()); print(sub_prof.round(3).to_string())

fig=plt.figure(figsize=(17,11)); gs=GridSpec(2,3,figure=fig,hspace=0.30,wspace=0.25,top=0.88,bottom=0.08)
ax=fig.add_subplot(gs[0,0])
for c,co in [('Bergmann','#d62728'),('Astrocyte','#2ca02c')]:
    m=ct==c; ax.scatter(U[m,0],U[m,1],s=2,c=co,linewidths=0,label=f'{c} ({m.sum()})',rasterized=True)
ax.set_title('A. Kozareva mouse cerebellum snRNA (611k nuclei)\nBergmann vs astrocyte',fontsize=11,loc='left'); ax.legend(fontsize=8,markerscale=3);ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP1',fontsize=8);ax.set_ylabel('UMAP2',fontsize=8)
axb=fig.add_subplot(gs[0,1]); scv=A.obs['S_Bergmann_shared'].values
p=axb.scatter(U[:,0],U[:,1],s=2,c=scv,cmap='viridis',vmin=np.percentile(scv,2),vmax=np.percentile(scv,98),linewidths=0,rasterized=True)
axb.set_title('B. Bergmann-shared signature',fontsize=11,loc='left');axb.set_xticks([]);axb.set_yticks([]); fig.colorbar(p,ax=axb,fraction=0.045)
axc=fig.add_subplot(gs[0,2]); scv2=A.obs['S_FC_astro_lean'].values
p2=axc.scatter(U[:,0],U[:,1],s=2,c=scv2,cmap='magma',vmin=np.percentile(scv2,2),vmax=np.percentile(scv2,98),linewidths=0,rasterized=True)
axc.set_title('C. Mouse FC astro-lean signature (Gfap/Agt/Ntsr2/Aqp4)',fontsize=10,loc='left');axc.set_xticks([]);axc.set_yticks([]); fig.colorbar(p2,ax=axc,fraction=0.045)
axd=fig.add_subplot(gs[1,0]); cols=plt.cm.tab10(np.linspace(0,1,10))
for i,cl in enumerate(sorted(Ab.obs['bsub'].unique(),key=int)):
    m=Ab.obs['bsub'].values==cl; axd.scatter(bgU[m,0],bgU[m,1],s=3,c=[cols[i%10]],linewidths=0,label=f'B{cl} ({m.sum()})')
axd.set_title(f'D. Bergmann reclustered ({Ab.obs.bsub.nunique()} subclusters)',fontsize=11,loc='left'); axd.legend(fontsize=7,markerscale=1.5,ncol=2);axd.set_xticks([]);axd.set_yticks([])
axe=fig.add_subplot(gs[1,1]); M=sub_prof.values
im=axe.imshow(M,cmap='RdBu_r',vmin=-0.5,vmax=0.5,aspect='auto')
axe.set_xticks(range(len(SIG))); axe.set_xticklabels(list(SIG.keys()),rotation=25,ha='right',fontsize=8)
axe.set_yticks(range(len(sub_prof))); axe.set_yticklabels([f'B{i}' for i in sub_prof.index],fontsize=8)
for i in range(M.shape[0]):
    for j in range(M.shape[1]): axe.text(j,i,f'{M[i,j]:.2f}',ha='center',va='center',fontsize=7)
axe.set_title('E. Signature per Bergmann subcluster (z)',fontsize=10,loc='left'); fig.colorbar(im,ax=axe,fraction=0.045)
axf=fig.add_subplot(gs[1,2])
data=[Ab.obs['S_FC_astro_lean'].values[Ab.obs['bsub'].values==cl] for cl in sorted(Ab.obs['bsub'].unique(),key=int)]
labs=[f'B{cl}' for cl in sorted(Ab.obs['bsub'].unique(),key=int)]
data.append(A.obs['S_FC_astro_lean'].values[ct=='Astrocyte']); labs.append('Astro')
axf.boxplot(data,labels=labs,showfliers=False)
axf.set_ylabel('FC astro-lean score (z)',fontsize=8); axf.set_title('F. Mouse FC signature: Bergmann subclusters vs astrocytes',fontsize=9.5,loc='left'); axf.tick_params(labelsize=8)
fig.suptitle('QC - Testing for an FC-like cluster in Kozareva 2021 mouse snRNA-seq',fontsize=14,y=0.945)
fig.text(0.5,0.915,'Same test as Fig. 8 in mouse: Bergmann and astrocytes are distinct; reclustering Bergmann tests whether a FC-like (astrocyte-lean) subpopulation exists among dissociated nuclei.',ha='center',fontsize=9,style='italic',color='#555')
out=FIG+"/QC_atlas_kozareva.png"; fig.savefig(out,dpi=200,bbox_inches='tight'); fig.savefig(out.replace('.png','.pdf'),bbox_inches='tight')
print("saved",out)
