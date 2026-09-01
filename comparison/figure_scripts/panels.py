"""Shared module-score panels + helpers for Fig 4/5/6 (unified across datasets)."""
import numpy as np, pandas as pd
PANELS={
 'astrocyte':['Aqp4','Agt','Gja1','Slc1a2','Mt1','Mt2','Slc1a3','Glul','S100b','Aldh1l1','Gjb6'],
 'bergmann':['Gdf10','Hopx','A2m','Timp4','Slc6a11','Gpr37l1','Ntsr2','Tnc','Fabp7'],
 'radial_BGlike':['Nes','Vim','Sox2','Ptprz1','Fabp7','Slc1a3','Hopx','Gdf10'],
 'reactive':['Gfap','Serpina3n','Lcn2','C3','Cxcl10','Steap4','Hspb1','Vim'],
 'contam_Purkinje':['Pcp2','Calb1','Car8','Itpr1','Pcp4','Prkcg'],
 'contam_MLI':['Gad1','Gad2','Pvalb','Sst','Slc32a1','Kit'],
 'contam_granule':['Gabra6','Neurod1','Rbfox3','Cadps2','Gng13'],
 'contam_oligoWM':['Mbp','Plp1','Mog','Mag','Cnp','Mobp'],
 'contam_vascular':['Cldn5','Flt1','Pecam1','Pdgfrb','Dcn','Vtn','Rgs5'],
}
def module_scores(E):
    """E: cells x genes (log-space). Returns DataFrame of z-scored panel means + dict of genes used."""
    Z=(E-E.mean())/(E.std()+1e-9)
    out={}; used={}
    for p,gs in PANELS.items():
        g=[x for x in gs if x in Z.columns]; used[p]=g
        out[p]=Z[g].mean(1).values if g else np.zeros(len(E))
    return pd.DataFrame(out,index=E.index), used
