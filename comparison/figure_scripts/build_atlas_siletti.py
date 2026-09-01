#!/usr/bin/env python
"""Build the Siletti 2023 non-spatial atlas validation table and QC figure.
Siletti 2023 human cerebellum snRNA: Bergmann glia (7684) + astrocytes (6578). Project the spatial
signatures; recluster Bergmann to test for a FC-like (astrocyte-lean) subpopulation."""
import os,numpy as np,pandas as pd,anndata as ad,scanpy as sc,scipy.sparse as sp
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from project_paths import SILETTI, FIGURES, TABLES
SIL=str(SILETTI)
FIG=str(FIGURES); TAB=str(TABLES)
# signatures (human orthologs, upper-case)
SIG={'Bergmann_shared':['HOPX','GPR37L1','A2M','GDF10','SLC1A3','TIMP4','FABP7','TNC','PTPRZ1'],
     'Parenchymal_astro':['AQP4','GJA1','SLC1A2','GLUL','S100B','ALDH1L1'],
     'FC_astro_lean':['GFAP','AGT','NTSR2','AQP4']}   # mouse Zhuang FC up-signature
a=ad.read_h5ad(SIL+"/results/combined.h5ad")
a.X=a.layers['lognorm'] if 'lognorm' in a.layers else a.X
ct=a.obs['cell_type'].astype(str).values; reg=a.obs['region'].astype(str).values
vn=set(map(str,a.var_names))
def score(A,genes):
    g=[x for x in genes if x in set(map(str,A.var_names))]
    X=A[:,g].X; X=X.toarray() if sp.issparse(X) else np.asarray(X)
    Z=(X-X.mean(0))/(X.std(0)+1e-9); return Z.mean(1)
for k,gs in SIG.items(): a.obs['S_'+k]=score(a,gs)
U=a.obsm['X_umap']
# ---- recluster Bergmann only ----
bg=a[a.obs['cell_type']=='Bergmann glial cell'].copy()
sc.pp.highly_variable_genes(bg,n_top_genes=2000); bgh=bg[:,bg.var.highly_variable].copy()
sc.pp.scale(bgh,max_value=10); sc.tl.pca(bgh,n_comps=30); sc.pp.neighbors(bgh,n_neighbors=15,n_pcs=30)
sc.tl.leiden(bgh,resolution=0.5,random_state=0,flavor="igraph",n_iterations=2,directed=False)
sc.tl.umap(bgh,random_state=0)
bg.obs['bsub']=bgh.obs['leiden'].values; bgU=bgh.obsm['X_umap']
for k in SIG: bg.obs['S_'+k]=bg.obs['S_'+k].values
nsub=bg.obs['bsub'].nunique()
# subcluster signature means
sub_prof=bg.obs.groupby('bsub')[['S_'+k for k in SIG]].mean()
sub_prof.to_csv(TAB+"/atlas_siletti_bergmann_subclusters.csv")
print("Bergmann leiden subclusters:",bg.obs['bsub'].value_counts().to_dict())
print(sub_prof.round(3).to_string())

# ================= FIGURE =================
fig=plt.figure(figsize=(17,11)); gs=GridSpec(2,3,figure=fig,hspace=0.30,wspace=0.25,top=0.88,bottom=0.08)
# A: UMAP cell_type
ax=fig.add_subplot(gs[0,0])
for c,co in [('Bergmann glial cell','#d62728'),('astrocyte','#2ca02c')]:
    m=ct==c; ax.scatter(U[m,0],U[m,1],s=3,c=co,linewidths=0,label=f'{c} ({m.sum()})',rasterized=True)
ax.set_title('A. Siletti human cerebellum (snRNA)\nBergmann vs astrocyte',fontsize=11,loc='left'); ax.legend(fontsize=8,markerscale=2); ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP1',fontsize=8);ax.set_ylabel('UMAP2',fontsize=8)
# B: UMAP region
axb=fig.add_subplot(gs[0,1])
for c,co in [('cortex','#1f77b4'),('deep_nuclei','#ff7f0e')]:
    m=reg==c; axb.scatter(U[m,0],U[m,1],s=3,c=co,linewidths=0,label=f'{c} ({m.sum()})',rasterized=True)
axb.set_title('B. Region (cortex vs deep nuclei)',fontsize=11,loc='left'); axb.legend(fontsize=8,markerscale=2); axb.set_xticks([]);axb.set_yticks([])
# C: UMAP Bergmann-signature score
axc=fig.add_subplot(gs[0,2]); scv=a.obs['S_Bergmann_shared'].values
p=axc.scatter(U[:,0],U[:,1],s=3,c=scv,cmap='viridis',vmin=np.percentile(scv,2),vmax=np.percentile(scv,98),linewidths=0,rasterized=True)
axc.set_title('C. Bergmann-shared signature (HOPX/GPR37L1/A2M/GDF10…)',fontsize=10,loc='left'); axc.set_xticks([]);axc.set_yticks([]); fig.colorbar(p,ax=axc,fraction=0.045)
# D: Bergmann-only UMAP by leiden subcluster
axd=fig.add_subplot(gs[1,0]); cols=plt.cm.tab10(np.linspace(0,1,10))
for i,cl in enumerate(sorted(bg.obs['bsub'].unique(),key=int)):
    m=bg.obs['bsub'].values==cl; axd.scatter(bgU[m,0],bgU[m,1],s=4,c=[cols[i%10]],linewidths=0,label=f'B{cl} ({m.sum()})')
axd.set_title(f'D. Bergmann reclustered ({nsub} subclusters)',fontsize=11,loc='left'); axd.legend(fontsize=7,markerscale=1.5,ncol=2); axd.set_xticks([]);axd.set_yticks([])
# E: Bergmann subcluster signature heatmap
axe=fig.add_subplot(gs[1,1]); M=sub_prof.values
im=axe.imshow(M,cmap='RdBu_r',vmin=-0.5,vmax=0.5,aspect='auto')
axe.set_xticks(range(len(SIG))); axe.set_xticklabels(list(SIG.keys()),rotation=25,ha='right',fontsize=8)
axe.set_yticks(range(len(sub_prof))); axe.set_yticklabels([f'B{i}' for i in sub_prof.index],fontsize=8)
for i in range(M.shape[0]):
    for j in range(M.shape[1]): axe.text(j,i,f'{M[i,j]:.2f}',ha='center',va='center',fontsize=7)
axe.set_title('E. Signature score per Bergmann subcluster (z)',fontsize=10,loc='left'); fig.colorbar(im,ax=axe,fraction=0.045)
# F: FC astro-lean score: Bergmann subclusters vs astrocytes
axf=fig.add_subplot(gs[1,2])
data=[bg.obs['S_FC_astro_lean'].values[bg.obs['bsub'].values==cl] for cl in sorted(bg.obs['bsub'].unique(),key=int)]
labs=[f'B{cl}' for cl in sorted(bg.obs['bsub'].unique(),key=int)]
data.append(a.obs['S_FC_astro_lean'].values[ct=='astrocyte']); labs.append('Astro')
axf.boxplot(data,labels=labs,showfliers=False)
axf.set_ylabel('FC astro-lean score (z)',fontsize=8); axf.set_title('F. Mouse FC signature in human Bergmann subclusters vs astrocytes',fontsize=9.5,loc='left'); axf.tick_params(labelsize=8)
fig.suptitle('Figure 8 - Testing for FC / BG as clusters in the Siletti (2023) human cerebellum snRNA atlas',fontsize=14,y=0.945)
fig.text(0.5,0.915,'Published atlas: Bergmann glia = ONE subcluster (3058), distinct from astrocytes. Reclustering Bergmann tests whether a FC-like (astrocyte-lean) subpopulation exists.',ha='center',fontsize=9,style='italic',color='#555')
out=FIG+"/QC_atlas_siletti.png"; fig.savefig(out,dpi=200,bbox_inches='tight'); fig.savefig(out.replace('.png','.pdf'),bbox_inches='tight')
print("saved",out)
