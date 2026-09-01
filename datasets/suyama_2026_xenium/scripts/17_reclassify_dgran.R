## 17 - Classify Bergmann cells by distance to the granule layer (d_gran)
suppressWarnings(suppressMessages({
  library(ggplot2); library(dplyr); library(RANN)
}))
project_root <- Sys.getenv("BGROOT")
if (!nzchar(project_root)) {
  script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  if (!length(script_arg)) stop("Set BGROOT or run this file with Rscript")
  project_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), "..", "..", "..") , mustWork=TRUE)
}
BASE <- file.path(project_root, "datasets/suyama_2026_xenium")

TAB <- file.path(BASE,"results/tables"); FIG <- file.path(BASE,"results/figures")
sp <- readRDS(file.path(BASE,"data","adult_laminar.rds"))$meta   # x,y,ctype,d_pkj,frac_*,cell
xy <- as.matrix(sp[,c("x","y")]); ct <- sp$ctype

## distances to nearest cell of each anatomical type
dnear <- function(type){
  ref <- xy[ct==type,,drop=FALSE]
  RANN::nn2(ref, xy, k=1)$nn.dists[,1]
}
sp$d_gran <- dnear("Granule")
sp$d_ml   <- dnear("ML_interneuron")
sp$d_pkj2 <- dnear("Purkinje")

bg <- sp[sp$ctype=="Bergmann",]
cat("Bergmann:", nrow(bg), "\n")
cat("d_gran (dist to nearest granule, um) quantiles:\n")
print(round(quantile(bg$d_gran, c(.1,.25,.5,.75,.9,.95,.98)),1))

## Purkinje-based class retained as an anatomical covariate.
purkinje <- read.csv(file.path(TAB,"adult_bergmann_purkinje_classification.csv"))
bg$purkinje_class <- purkinje$ectopic[match(bg$cell, purkinje$cell)]
cat("\nDistance-to-granule among Purkinje-classified PCL Bergmann:\n")
pcl <- bg[bg$purkinje_class=="Canonical_PCL",]
for (thr in c(20,25,30,40,50)) cat(sprintf("  d_gran>%d um: %d cells (%.1f%%)\n",
    thr, sum(pcl$d_gran>thr), 100*mean(pcl$d_gran>thr)))

## ---- Laminar classification by ML depth (d_gran) ----
## ML has NO granule cells -> far from granule = in ML (regardless of d_pkj)
bg$laminar_class <- with(bg, ifelse(frac_oligo>=0.25, "WM/Ambiguous",
  ifelse(d_gran>=30, "FC",                       # clearly above granule layer = ML
  ifelse(d_gran<=16, "BG_PCL",                       # at PCL/GCL interface
         "Intermediate"))))
## flag GCL-sunk (below PCL, inside granule mass): close to granule but far from Purkinje & ML
bg$laminar_class[bg$laminar_class=="BG_PCL" & bg$d_pkj2>40 & bg$d_ml>60] <- "BG_GCL"
cat("\n=== Laminar classification (by d_gran) ===\n"); print(table(bg$laminar_class))
cat("\nPurkinje class vs laminar class:\n"); print(table(bg$purkinje_class, bg$laminar_class))
write.csv(bg, file.path(TAB,"adult_bergmann_laminar_granule.csv"), row.names=FALSE)

## ---- overview map ----
ctx <- sp[sp$ctype %in% c("Purkinje","Granule","ML_interneuron","Bergmann","Astro_parenchymal"),]
cols <- c(BG_PCL="#1B998B", Intermediate="#9ecae1", BG_GCL="#F18701",
          FC="#D7263D", `WM/Ambiguous`="grey60")
ov <- ggplot()+
  geom_point(data=ctx, aes(x,y), color="grey88", size=0.3, stroke=0)+
  geom_point(data=subset(sp,ctype=="Purkinje"), aes(x,y), color="black", size=0.35, stroke=0)+
  geom_point(data=subset(bg,laminar_class=="BG_PCL"), aes(x,y), color=cols["BG_PCL"], size=0.5, stroke=0)+
  geom_point(data=subset(bg,laminar_class=="Intermediate"), aes(x,y), color=cols["Intermediate"], size=0.5, stroke=0)+
  geom_point(data=subset(bg,laminar_class=="BG_GCL"), aes(x,y), color=cols["BG_GCL"], size=1.4, stroke=0)+
  geom_point(data=subset(bg,laminar_class=="FC"), aes(x,y), color=cols["FC"], size=1.5, stroke=0)+
  coord_fixed()+theme_void(base_size=11)+
  ggtitle(sprintf("Classification laminaire (d_gran): FC=%d, Intermediate=%d, BG_PCL=%d, BG_GCL=%d",
          sum(bg$laminar_class=="FC"),sum(bg$laminar_class=="Intermediate"),
          sum(bg$laminar_class=="BG_PCL"),sum(bg$laminar_class=="BG_GCL")))+
  theme(plot.title=element_text(face="bold",size=11))
ggsave(file.path(FIG,"bergmann_classification_overview.png"), ov, width=14, height=12, dpi=170, bg="white")

## ---- re-tile 3x3 with new classes ----
OUT <- file.path(FIG,"verify_tiles_granule_distance"); dir.create(OUT, showWarnings=FALSE)
bx <- range(bg$x)+c(-60,60); by <- range(bg$y)+c(-60,60)
xs <- seq(bx[1],bx[2],length.out=4); ys <- seq(by[1],by[2],length.out=4)
tid <- 0
for (i in 1:3) for (j in 1:3){
  xr<-c(xs[i],xs[i+1]); yr<-c(ys[j],ys[j+1])
  inb <- bg$x>=xr[1]&bg$x<=xr[2]&bg$y>=yr[1]&bg$y<=yr[2]
  if (sum(inb)<10) next
  tid<-tid+1
  c2<-subset(ctx,x>=xr[1]&x<=xr[2]&y>=yr[1]&y<=yr[2])
  pk<-subset(sp,ctype=="Purkinje"&x>=xr[1]&x<=xr[2]&y>=yr[1]&y<=yr[2])
  bb<-bg[inb,]
  p<-ggplot()+
    geom_point(data=c2,aes(x,y),color="grey86",size=1.4,stroke=0)+
    geom_point(data=pk,aes(x,y),color="black",size=2,stroke=0)+
    geom_point(data=subset(bb,laminar_class=="BG_PCL"),aes(x,y),color=cols["BG_PCL"],size=2.6,stroke=0)+
    geom_point(data=subset(bb,laminar_class=="Intermediate"),aes(x,y),color=cols["Intermediate"],size=2.6,stroke=0)+
    geom_point(data=subset(bb,laminar_class=="BG_GCL"),aes(x,y),shape=21,fill=cols["BG_GCL"],color="white",size=4,stroke=.5)+
    geom_point(data=subset(bb,laminar_class=="FC"),aes(x,y),shape=21,fill=cols["FC"],color="white",size=4.5,stroke=.5)+
    coord_fixed(xlim=xr,ylim=yr)+theme_void(base_size=12)+
    theme(plot.title=element_text(face="bold",size=12),panel.border=element_rect(color="grey60",fill=NA),
          plot.background=element_rect(fill="white",color=NA))+
    ggtitle(sprintf("T%d | ROUGE=FC(%d) orange=BG_GCL teal=PCL noir=Purkinje",
            tid, sum(bb$laminar_class=="FC")))
  ggsave(file.path(OUT,sprintf("tile_T%02d.png",tid)),p,width=11,height=11,dpi=200,bg="white")
}
cat("\nDONE_17\n")
