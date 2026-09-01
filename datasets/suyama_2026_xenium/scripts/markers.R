## Shared marker panels for Suyama2026 cerebellum Xenium analysis
## (mouse gene symbols, as in the Xenium 5K panel)

markers <- list(
  PCL_Purkinje   = c("Pcp2","Calb1","Car8","Itpr1","Grid2","Aldoc"),
  ML_interneuron = c("Gad1","Gad2","Pvalb","Sst","Reln","Slc32a1","Syn1","Snap25"),
  GCL_granule    = c("Gabra6","Gng13","Neurod1","Barhl1","Zic1","Zic2"),
  WM_oligo       = c("Mbp","Plp1","Mog","Mag","Cnp"),
  Astro_Bergmann = c("Aldh1l1","Sox9","Slc1a3","Slc1a2","S100b","Gja1","Gjb6",
                     "Glul","Aqp4","Fabp7","Hopx","Gpr37l1","Ednrb","Notch2"),
  Microglia      = c("Cx3cr1","Aif1","C1qa","C1qb","Csf1r","Tmem119","P2ry12"),
  Vascular       = c("Cldn5","Pecam1","Flt1","Pdgfrb","Rgs5","Cldn11"),
  Proliferation  = c("Mki67","Pcna","Top2a","Ccnb1")
)

## Module-score gene sets (signatures)
sig <- list(
  ## Pan-astrocyte identity
  astrocyte = c("Aldh1l1","Sox9","Slc1a3","Slc1a2","Gja1","Glul","Aqp4","Fabp7","Gjb6"),
  ## Bergmann glia (canonical, radial)
  bergmann  = c("Gdf10","Hopx","Fabp7","Slc1a3","Aqp4","Gria1","Sept4","A2m",
                "Timp4","Vim","Tnc","Gpr37l1"),
  ## Radial / BG-like progenitor-leaning
  bg_like   = c("Gdf10","Hopx","Fabp7","Sox2","Sox9","Nes","Tnc","Vim","Slc1a3"),
  ## Reactive astrocyte
  reactive  = c("Gfap","Vim","Serpina3n","Lcn2","Cp","Timp1","C3","S100b","Cd44"),
  ## IGL / parenchymal astrocyte (to distinguish from Bergmann)
  igl_astro = c("Gfap","Aqp4","Slc1a2","Mfge8","Agt","Slc6a11")
)

## Contamination signatures (other layers / cell types)
contam <- list(
  Purkinje    = c("Pcp2","Calb1","Car8","Itpr1","Pcp4"),
  Granule     = c("Gabra6","Neurod1","Rbfox3","Cntn2","Gng13"),
  Interneuron = c("Gad1","Gad2","Pvalb","Slc32a1","Reln"),
  Oligo       = c("Mbp","Plp1","Mog","Mag","Cnp"),
  Microglia   = c("Cx3cr1","C1qa","Csf1r","P2ry12"),
  Vascular    = c("Cldn5","Pecam1","Flt1","Rgs5","Pdgfrb")
)
