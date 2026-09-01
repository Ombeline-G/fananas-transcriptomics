## 30 - Export adult UMAP coords + normalized marker expression (for comparison Figures 1 & 2)
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
cb <- readRDS(file.path(BASE,"data","adult_cerebellum.rds"))
DefaultAssay(cb) <- "Xenium"
cat("adult object:", ncol(cb), "cells x", nrow(cb), "genes\n")
cat("reductions:", paste(Reductions(cb),collapse=", "), "\n")

## UMAP embeddings
um <- Embeddings(cb, "umap")
md <- data.frame(cell=colnames(cb), cl=as.character(cb$cl),
                 UMAP_1=um[,1], UMAP_2=um[,2], stringsAsFactors=FALSE)
write.csv(md, file.path(TAB,"adult_umap_export.csv"), row.names=FALSE)
cat("wrote adult_umap_export.csv\n")

## normalized expression for marker panels (only genes present in the 5K panel)
want <- c(
 # PCL
 "Pcp2","Calb1","Car8","Itpr1","Pcp4","Prkcg","Grid2","Aldoc","Homer3",
 # ML interneurons
 "Gad1","Gad2","Pvalb","Sst","Reln","Slc32a1","Kit","Sorcs3","Nxph1","Kctd12","Lypd6",
 # GCL
 "Gabra6","Gng13","Neurod1","Barhl1","Zic1","Zic2","Rbfox3","Cadps2","Etv1",
 # WM / oligo
 "Mbp","Plp1","Mog","Mag","Cnp","Mobp",
 # astroglia / Bergmann
 "Aldh1l1","Sox9","Slc1a3","Slc1a2","S100b","Gja1","Gjb6","Glul","Aqp4","Fabp7",
 "Hopx","Gpr37l1","Ednrb","Notch2","Gdf10","Agt","A2m","Timp4","Slc6a11",
 # other cell types (for Fig1 typing if needed)
 "Cldn5","Flt1","Pecam1","Pdgfrb","Vtn","Dcn","Cx3cr1","C1qa","C1qb","Csf1r",
 "Pdgfra","Cspg4","Olig1","Olig2","Sox10")
g <- intersect(want, rownames(cb))
cat("genes exported:", length(g), "/", length(want), "\n")
M <- t(as.matrix(GetAssayData(cb, assay="Xenium", layer="data")[g, , drop=FALSE]))
df <- data.frame(cell=colnames(cb), M, check.names=FALSE)
write.csv(df, file.path(TAB,"adult_marker_expr.csv"), row.names=FALSE)
cat("wrote adult_marker_expr.csv (", nrow(df), "cells x", length(g), "genes )\n")
cat("DONE_30\n")
