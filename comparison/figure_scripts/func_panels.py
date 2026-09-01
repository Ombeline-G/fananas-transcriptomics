"""Bergmann/astrocyte FUNCTIONAL modules (Fig 7) — established glial functions."""
FUNC={
 'Glutamate clearance':['Slc1a3','Slc1a2','Glul','Slc1a4'],
 'K+/water homeostasis':['Kcnj10','Kcnj16','Aqp4','Kcnn2'],
 'Gap-junction coupling':['Gja1','Gjb6'],
 'Ca2+/glutamate-receptor':['Gria1','Gria4','Itpr2','Grid2','Grm3','Grm5'],
 'Energy/lactate metabolism':['Ldha','Ldhb','Slc16a1','Aldoc','Pkm','Slc2a1'],
 'Lipid/cholesterol':['Apoe','Fabp7','Scd2','Apod','Clu'],
 'Radial fibre/cytoskeleton':['Gfap','Vim','Nes','Ptprz1'],
 'ECM/adhesion':['Tnc','Hepacam','Mlc1','Sparcl1','Ncan'],
 'Bergmann trophic/identity':['Gdf10','Gpr37l1','Ntsr2','Timp4','A2m','Hopx','Sox9'],
 'Reactive/immune':['Gfap','Serpina3n','C3','Cxcl10','Lcn2','Vim'],
}
ALL_FUNC=sorted({g for gs in FUNC.values() for g in gs})
