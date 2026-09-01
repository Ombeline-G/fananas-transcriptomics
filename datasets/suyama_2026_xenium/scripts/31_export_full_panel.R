## 31 - Export adult normalized expression for the FULL comparison marker panel (Fig 4/5/6)
suppressWarnings(suppressMessages({
library(Seurat); library(SeuratObject)
}))
project_root <- Sys.getenv("BGROOT")
if (!nzchar(project_root)) {
  script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  if (!length(script_arg)) stop("Set BGROOT or run this file with Rscript")
  project_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), "..", "..", "..") , mustWork=TRUE)
}
BASE <- file.path(project_root, "datasets/suyama_2026_xenium")

TAB  <- file.path(BASE,"results","tables")
cb <- readRDS(file.path(BASE,"data","adult_cerebellum.rds")); DefaultAssay(cb) <- "Xenium"
want <- c(
 "Pcp2","Calb1","Car8","Itpr1","Pcp4","Prkcg","Grid2","Aldoc","Homer3",
 "Gad1","Gad2","Pvalb","Sst","Reln","Slc32a1","Kit","Sorcs3","Nxph1","Kctd12","Lypd6",
 "Gabra6","Gng13","Neurod1","Barhl1","Zic1","Zic2","Rbfox3","Cadps2","Etv1",
 "Mbp","Plp1","Mog","Mag","Cnp","Mobp",
 "Aqp4","Agt","Gja1","Slc1a2","Mt1","Mt2","Slc1a3","Glul","S100b","Aldh1l1","Gjb6",
 "Gdf10","Hopx","A2m","Timp4","Slc6a11","Gpr37l1","Ntsr2","Tnc","Fabp7",
 "Nes","Vim","Sox2","Ptprz1","Sox9","Notch2","Ednrb",
 "Gfap","Serpina3n","Lcn2","C3","Cxcl10","Steap4","Hspb1",
 "Cldn5","Flt1","Pecam1","Pdgfrb","Dcn","Vtn","Rgs5","Col1a1","Igf2",
 "Cx3cr1","C1qa","C1qb","Csf1r","P2ry12")
g <- intersect(want, rownames(cb))
cat("panel genes present:", length(g), "/", length(want), "\n")
cat("MISSING:", paste(setdiff(want,rownames(cb)),collapse=", "), "\n")
M <- t(as.matrix(GetAssayData(cb, assay="Xenium", layer="data")[g, , drop=FALSE]))
df <- data.frame(cell=colnames(cb), M, check.names=FALSE)
write.csv(df, file.path(TAB,"adult_panel_expr.csv"), row.names=FALSE)
cat("wrote adult_panel_expr.csv (", nrow(df), "x", length(g), ")\nDONE_31\n")
