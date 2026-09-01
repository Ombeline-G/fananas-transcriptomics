## 19 - Two perfectly-aligned PNGs for swipe validation (anatomy vs BG classification)
suppressWarnings(suppressMessages({
  library(ggplot2); library(dplyr)
}))
project_root <- Sys.getenv("BGROOT")
if (!nzchar(project_root)) {
  script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  if (!length(script_arg)) stop("Set BGROOT or run this file with Rscript")
  project_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), "..", "..", "..") , mustWork=TRUE)
}
BASE <- file.path(project_root, "datasets/suyama_2026_xenium")

TAB <- file.path(BASE,"results/tables"); OUT <- file.path(BASE,"results/swipe")
dir.create(OUT, showWarnings=FALSE)

sp <- readRDS(file.path(BASE,"data","adult_laminar.rds"))$meta
bg <- read.csv(file.path(TAB,"adult_bergmann_laminar_granule.csv"))

## common extent (cerebellar cortex around the Bergmann)
xr <- range(bg$x)+c(-80,80); yr <- range(bg$y)+c(-80,80)
ctx <- sp[sp$x>=xr[1]&sp$x<=xr[2]&sp$y>=yr[1]&sp$y<=yr[2],]
pk  <- subset(ctx, ctype=="Purkinje")
mli <- subset(ctx, ctype=="ML_interneuron")
oth <- subset(ctx, !ctype %in% c("Purkinje","Bergmann"))
W<-11; H<-13; D<-150
base_layers <- function(){
  list(
    geom_point(data=oth, aes(x,y), color="grey85", size=0.45, stroke=0),
    geom_point(data=mli, aes(x,y), color="#9ecae1", size=0.5, stroke=0),   # ML interneurons mark the ML
    geom_point(data=pk,  aes(x,y), color="black",  size=0.7, stroke=0),
    coord_fixed(xlim=xr, ylim=yr, expand=FALSE),
    theme_void(base_size=10),
    theme(plot.background=element_rect(fill="white",color=NA),
          plot.margin=margin(2,2,2,2))
  )
}

## ---- IMAGE A : anatomy (Bergmann shown neutral grey) ----
pA <- ggplot()+ base_layers()+
  geom_point(data=bg, aes(x,y), color="#555555", size=1.1, stroke=0)
ggsave(file.path(OUT,"A_anatomy.png"), pA, width=W, height=H, dpi=D, bg="white")

## ---- IMAGE B : BG classification overlay ----
cols <- c(BG_PCL="#1B998B", Intermediate="#9970c0", BG_GCL="#F18701", FC="#D7263D")
pB <- ggplot()+ base_layers()+
  geom_point(data=subset(bg,laminar_class=="BG_PCL"),       aes(x,y), color=cols["BG_PCL"], size=1.1, stroke=0)+
  geom_point(data=subset(bg,laminar_class=="Intermediate"), aes(x,y), color=cols["Intermediate"], size=1.1, stroke=0)+
  geom_point(data=subset(bg,laminar_class=="BG_GCL"),       aes(x,y), color=cols["BG_GCL"], size=2.0, stroke=0)+
  geom_point(data=subset(bg,laminar_class=="FC"),        aes(x,y), color=cols["FC"], size=2.2, stroke=0)
ggsave(file.path(OUT,"B_classification.png"), pB, width=W, height=H, dpi=D, bg="white")

cat("Image px (A):", paste(dim(png::readPNG(file.path(OUT,"A_anatomy.png")))[1:2],collapse="x"),"\n")
cat("FC:",sum(bg$laminar_class=="FC")," BG_PCL:",sum(bg$laminar_class=="BG_PCL"),
    " Intermediate:",sum(bg$laminar_class=="Intermediate")," BG_GCL:",sum(bg$laminar_class=="BG_GCL"),"\n")
cat("DONE_19\n")
