## 08 - Crop adult cerebellum, QC, cluster, validate layers
suppressWarnings(suppressMessages({
  library(Seurat); library(SeuratObject); library(ggplot2); library(dplyr)
  library(patchwork); library(viridis)
}))
set.seed(777)
project_root <- Sys.getenv("BGROOT")
if (!nzchar(project_root)) {
  script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  if (!length(script_arg)) stop("Set BGROOT or run this file with Rscript")
  project_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), "..", "..", "..") , mustWork=TRUE)
}
BASE <- file.path(project_root, "datasets/suyama_2026_xenium")

TAB  <- file.path(BASE,"results","tables"); FIG <- file.path(BASE,"results","figures")
source(file.path(BASE,"scripts","markers.R"))

full <- readRDS(file.path(BASE,"data","full_slide.rds"))
DefaultAssay(full) <- "Xenium"
co <- GetTissueCoordinates(full[["fov"]])
ccn <- intersect(c("cell","cells"),colnames(co))[1]; co$cell <- co[[ccn]]

## crop box around adult cerebellum (low-y, high-x dense foliated region)
BOX <- list(xmin=7800, xmax=10300, ymin=100, ymax=3400)
inbox <- co$cell[co$x>=BOX$xmin & co$x<=BOX$xmax & co$y>=BOX$ymin & co$y<=BOX$ymax]
cat("Cells in crop box:", length(inbox), "\n")
cb <- subset(full, cells = inbox)

## QC distribution
cat("nCount quantiles:\n"); print(quantile(cb$nCount_Xenium, c(.1,.25,.5,.75,.9)))
cat("nFeature quantiles:\n"); print(quantile(cb$nFeature_Xenium, c(.1,.25,.5,.75,.9)))
cb <- subset(cb, subset = nCount_Xenium>=150 & nFeature_Xenium>=80)
cat("After QC:", ncol(cb), "cells\n")

## process
cb <- NormalizeData(cb, verbose=FALSE)
cb <- FindVariableFeatures(cb, nfeatures=2000, verbose=FALSE)
cb <- ScaleData(cb, features=rownames(cb), verbose=FALSE)
cb <- RunPCA(cb, npcs=30, verbose=FALSE)
cb <- RunUMAP(cb, dims=1:30, verbose=FALSE)
cb <- FindNeighbors(cb, dims=1:30, verbose=FALSE)
cb <- FindClusters(cb, resolution=1.0, random.seed=777, verbose=FALSE)
cb$cl <- paste0("c", as.character(Idents(cb)))
Idents(cb) <- "cl"
cat("\n=== clusters ===\n"); print(table(cb$cl))

## module scores (signatures + contamination + layer)
for (nm in names(sig)) { gs<-sig[[nm]][sig[[nm]] %in% rownames(cb)]
  cb<-AddModuleScore(cb,features=list(gs),name=paste0("s_",nm),seed=777)
  names(cb@meta.data)[names(cb@meta.data)==paste0("s_",nm,"1")]<-paste0("score_",nm) }
for (nm in names(contam)) { gs<-contam[[nm]][contam[[nm]] %in% rownames(cb)]
  cb<-AddModuleScore(cb,features=list(gs),name=paste0("c_",nm),seed=777)
  names(cb@meta.data)[names(cb@meta.data)==paste0("c_",nm,"1")]<-paste0("contam_",nm) }
## layer scores
layer_sig <- list(PCL=markers$PCL_Purkinje, ML=markers$ML_interneuron,
                  GCL=markers$GCL_granule, WM=markers$WM_oligo)
for (nm in names(layer_sig)) { gs<-layer_sig[[nm]][layer_sig[[nm]] %in% rownames(cb)]
  cb<-AddModuleScore(cb,features=list(gs),name=paste0("L_",nm),seed=777)
  names(cb@meta.data)[names(cb@meta.data)==paste0("L_",nm,"1")]<-paste0("layer_",nm) }

scoreCols <- c(paste0("score_",names(sig)), paste0("layer_",names(layer_sig)),
               paste0("contam_",names(contam)))
sm <- cb@meta.data %>% group_by(cl) %>%
  summarise(n=n(), across(all_of(scoreCols), ~round(mean(.x),2)), .groups="drop")
cat("\n=== per-cluster mean scores ===\n"); print(as.data.frame(sm), row.names=FALSE)
write.csv(sm, file.path(TAB,"adult_cluster_scores.csv"), row.names=FALSE)

## cluster markers
mk <- FindAllMarkers(cb, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.5, verbose=FALSE)
write.csv(mk, file.path(TAB,"adult_cluster_markers.csv"), row.names=FALSE)
top <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC,n=12)
cat("\n=== top markers per cluster ===\n")
for (cl in sort(unique(top$cluster))) cat(cl,":",paste(top$gene[top$cluster==cl],collapse=", "),"\n")

## spatial coords for crop
spco <- co[match(colnames(cb), co$cell),]
md <- cb@meta.data; md$cell<-rownames(md)
sp <- left_join(md, spco[,c("cell","x","y")], by="cell")
write.csv(sp[,c("cell","x","y","cl",scoreCols)], file.path(TAB,"adult_cells_scores.csv"), row.names=FALSE)

## FIG: spatial layer markers in crop
spm <- function(g){ e<-FetchData(cb,vars=g,layer="data")[,1]; d<-sp; d$v<-e[match(d$cell,colnames(cb))]
  d<-d[order(d$v),]; ggplot(d,aes(x,y,color=v))+geom_point(size=0.6,stroke=0)+coord_fixed()+
  scale_color_viridis(option="magma")+theme_void(base_size=9)+ggtitle(g) }
lg <- intersect(c("Calb1","Pcp2","Gad1","Gabra6","Mog","Gdf10","Hopx","Aqp4","Slc1a3"),rownames(cb))
ggsave(file.path(FIG,"15_adult_layer_markers.png"),
       wrap_plots(lapply(lg,spm),ncol=3), width=15, height=5*ceiling(length(lg)/3), dpi=140, bg="white")

## FIG: clusters spatial + umap
psp <- ggplot(sp,aes(x,y,color=cl))+geom_point(size=0.6,stroke=0)+coord_fixed()+
  theme_void(base_size=9)+ggtitle("Adult CB clusters (spatial)")+
  guides(color=guide_legend(override.aes=list(size=2)))
pum <- DimPlot(cb,group.by="cl",label=TRUE)+ggtitle("Adult CB UMAP")
ggsave(file.path(FIG,"16_adult_clusters.png"), (psp|pum), width=15, height=7, dpi=140, bg="white")

saveRDS(cb, file.path(BASE,"data","adult_cerebellum.rds"))
cat("\nDONE_08\n")
