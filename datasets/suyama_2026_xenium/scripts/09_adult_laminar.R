## 09 - Adult: annotate cell types, spatial neighborhood, classify BG_PCL/FC/Ambiguous
suppressWarnings(suppressMessages({
  library(Seurat); library(SeuratObject); library(ggplot2); library(dplyr)
  library(patchwork); library(viridis); library(RANN)
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

cb <- readRDS(file.path(BASE,"data","adult_cerebellum.rds"))

## ---- cluster -> cell type ----
map <- c(
  c15="Purkinje",
  c0="Granule", c1="Granule", c2="Granule",
  c7="UBC",
  c9="ML_interneuron", c12="ML_interneuron", c16="ML_interneuron",
  c14="Bergmann", c10="Bergmann",
  c11="Astro_parenchymal", c13="Astro_parenchymal", c3="Astro_parenchymal",
  c6="Oligo", c17="Oligo", c18="OPC",
  c19="Microglia", c8="Endothelial", c23="Mural",
  c20="Meninges", c21="Choroid", c22="Ependymal",
  c4="Other_neuron", c5="Other_neuron", c24="Other_neuron")
cb$ctype <- unname(map[cb$cl])
cb$bgsub <- ifelse(cb$cl=="c14","BG_c14", ifelse(cb$cl=="c10","BG_c10", NA))
cat("=== cell types ===\n"); print(sort(table(cb$ctype), decreasing=TRUE))

## ---- spatial coords ----
sp <- read.csv(file.path(TAB,"adult_cells_scores.csv"))
sp <- sp[match(colnames(cb), sp$cell),]
xy <- as.matrix(sp[,c("x","y")])
ct <- cb$ctype

## ---- KNN neighborhood (k=31, drop self) ----
k <- 31
nn <- nn2(xy, xy, k=k)
idx <- nn$nn.idx[,-1]; dst <- nn$nn.dists[,-1]
neighType <- matrix(ct[idx], ncol=k-1)
frac <- function(type) rowMeans(neighType==type, na.rm=TRUE)
sp$frac_gran  <- frac("Granule")
sp$frac_ml    <- frac("ML_interneuron")
sp$frac_pkj   <- frac("Purkinje")
sp$frac_oligo <- frac("Oligo")
sp$frac_berg  <- frac("Bergmann")
sp$frac_astro <- frac("Astro_parenchymal")
## distance to nearest Purkinje cell
pkj_xy <- xy[ct=="Purkinje",,drop=FALSE]
dp <- nn2(pkj_xy, xy, k=1)$nn.dists[,1]
sp$d_pkj <- dp
sp$ctype <- ct; sp$bgsub <- cb$bgsub

## ---- diagnostics on Bergmann cells ----
bg <- sp[sp$ctype=="Bergmann",]
cat("\nBergmann cells:", nrow(bg), "\n")
cat("d_pkj quantiles (um):\n"); print(round(quantile(bg$d_pkj, c(.1,.25,.5,.75,.9)),1))
cat("frac_ml quantiles:\n");  print(round(quantile(bg$frac_ml, c(.1,.25,.5,.75,.9)),2))
cat("frac_gran quantiles:\n");print(round(quantile(bg$frac_gran,c(.1,.25,.5,.75,.9)),2))
cat("frac_oligo quantiles:\n");print(round(quantile(bg$frac_oligo,c(.1,.25,.5,.75,.9)),2))
cat("median NN distance (um):", round(median(dst[,1]),1), "\n")

## ---- laminar classification of Bergmann ----
classify <- function(r){
  if (r["frac_oligo"]>=0.25) return("Ambiguous")          # near WM
  if (r["d_pkj"]<=25 || r["frac_pkj"]>0) return("BG_PCL")  # adjacent to Purkinje soma
  if (r["frac_gran"]>=0.20 && r["d_pkj"]<=45) return("BG_PCL") # at PCL/GCL border
  if (r["frac_ml"]>=0.30 && r["frac_gran"]<=0.05 && r["d_pkj"]>=45) return("FC") # displaced into ML
  return("Ambiguous")
}
bg$bgclass <- apply(bg[,c("frac_oligo","d_pkj","frac_pkj","frac_gran","frac_ml")],1,classify)
cat("\n=== BG laminar classes ===\n"); print(table(bg$bgclass))
cat("\n=== BG class x transcriptomic subcluster ===\n"); print(table(bg$bgclass, bg$bgsub))
write.csv(bg, file.path(TAB,"adult_bergmann_laminar_purkinje.csv"), row.names=FALSE)

## ---- MAP: cortex context + BG classes ----
cortex <- sp[sp$ctype %in% c("Purkinje","Granule","ML_interneuron","Bergmann",
                             "Astro_parenchymal","Oligo"),]
base_cols <- c(Purkinje="#000000", Granule="grey80", ML_interneuron="#9ecae1",
               Astro_parenchymal="#c7e9c0", Oligo="#fdd0a2", Bergmann="grey50")
pA <- ggplot(cortex, aes(x,y,color=ctype)) + geom_point(size=0.35,stroke=0) + coord_fixed()+
  scale_color_manual(values=base_cols) + theme_void(base_size=9)+
  ggtitle("Adult cerebellar cortex (context)")+guides(color=guide_legend(override.aes=list(size=2)))
pB <- ggplot() +
  geom_point(data=cortex, aes(x,y), color="grey88", size=0.3, stroke=0)+
  geom_point(data=subset(sp,ctype=="Purkinje"), aes(x,y), color="black", size=0.4, stroke=0)+
  geom_point(data=bg, aes(x,y,color=bgclass), size=0.9, stroke=0)+
  scale_color_manual(values=c(BG_PCL="#1B998B", FC="#D7263D", Ambiguous="#F2C14E"))+
  coord_fixed()+theme_void(base_size=9)+ggtitle("Bergmann laminar class (BG_PCL/FC/Ambiguous)")+
  guides(color=guide_legend(override.aes=list(size=3)))
ggsave(file.path(FIG,"17_adult_BG_PCL_vs_ML_map.png"), (pA|pB), width=17, height=8, dpi=160, bg="white")

## ---- zoom to a couple folia for clarity ----
xr <- quantile(bg$x,c(.45,.8)); yr <- quantile(bg$y,c(.1,.5))
pz <- ggplot()+
  geom_point(data=subset(cortex,x>=xr[1]&x<=xr[2]&y>=yr[1]&y<=yr[2]),aes(x,y,color=ctype),size=0.8,stroke=0)+
  scale_color_manual(values=base_cols)+
  geom_point(data=subset(sp,ctype=="Purkinje"&x>=xr[1]&x<=xr[2]&y>=yr[1]&y<=yr[2]),aes(x,y),color="black",size=1)+
  geom_point(data=subset(bg,x>=xr[1]&x<=xr[2]&y>=yr[1]&y<=yr[2]),aes(x,y,fill=bgclass),
             shape=21,size=2.2,stroke=0.3,color="white")+
  scale_fill_manual(values=c(BG_PCL="#1B998B",FC="#D7263D",Ambiguous="#F2C14E"))+
  coord_fixed()+theme_void(base_size=10)+ggtitle("Zoom: Bergmann laminar position")
ggsave(file.path(FIG,"18_adult_BG_zoom.png"), pz, width=12, height=9, dpi=160, bg="white")

saveRDS(list(meta=sp, bg=bg), file.path(BASE,"data","adult_laminar.rds"))
cat("\nDONE_09\n")
