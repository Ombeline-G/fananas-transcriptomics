## 11 - Focused analysis of ECTOPIC Bergmann (soma displaced from the PCL)
suppressWarnings(suppressMessages({
  library(Seurat); library(SeuratObject); library(ggplot2); library(dplyr)
  library(patchwork); library(ggrepel); library(viridis)
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

bg <- read.csv(file.path(TAB,"adult_bergmann_laminar_purkinje.csv"))   # per Bergmann cell, has d_pkj, frac_*
cb <- readRDS(file.path(BASE,"data","adult_cerebellum.rds")); DefaultAssay(cb)<-"Xenium"

## ---- define ectopic groups (continuum on distance to Purkinje band) ----
bg$effdir <- with(bg, ifelse(frac_ml>0.25 & frac_gran<0.20, "ML",
                       ifelse(frac_gran>0.40, "GCL", "mixed")))
bg$ectopic <- with(bg, ifelse(d_pkj<20, "Canonical_PCL",
                       ifelse(d_pkj>40 & effdir=="ML",  "Ectopic_FC",
                       ifelse(d_pkj>40 & effdir=="GCL", "Ectopic_GCL", "Intermediate"))))
cat("=== ectopic groups ===\n"); print(table(bg$ectopic))

## ---- module scores already present in bg (from adult_cells_scores.csv) ----
cat("\n=== mean scores by ectopic group ===\n")
print(bg %>% group_by(ectopic) %>%
  summarise(n=n(), d_pkj=round(median(d_pkj),1),
            bergmann=round(mean(score_bergmann),2), astro=round(mean(score_astrocyte),2),
            reactive=round(mean(score_reactive),2), igl=round(mean(score_igl_astro),2),
            .groups="drop") %>% as.data.frame(), row.names=FALSE)
write.csv(bg, file.path(TAB,"adult_bergmann_purkinje_classification.csv"), row.names=FALSE)

## ---- DE: Ectopic_FC vs Canonical_PCL Bergmann ----
cb$ecto <- bg$ectopic[match(colnames(cb), bg$cell)]
Idents(cb) <- "ecto"
run_de <- function(g1,g2,file){
  d<-FindMarkers(cb, ident.1=g1, ident.2=g2, logfc.threshold=0, min.pct=0.1)
  d$gene<-rownames(d); d<-d%>%arrange(p_val_adj,desc(abs(avg_log2FC)))
  write.csv(d, file.path(TAB,file), row.names=FALSE); d }
deML <- run_de("Ectopic_FC","Canonical_PCL","DE_adult_ectopicFC_vs_canonicalPCL.csv")
nML <- sum(cb$ecto=="Ectopic_FC",na.rm=T); nPCL<-sum(cb$ecto=="Canonical_PCL",na.rm=T)
sig <- deML%>%filter(p_val_adj<0.05 & abs(avg_log2FC)>0.5)
cat(sprintf("\nEctopic_FC (n=%d) vs Canonical_PCL (n=%d): %d sig genes (padj<0.05,|lfc|>0.5)\n",
            nML,nPCL,nrow(sig)))
upML<-deML%>%filter(p_val_adj<0.05&avg_log2FC>0.5)%>%arrange(desc(avg_log2FC))
dnML<-deML%>%filter(p_val_adj<0.05&avg_log2FC< -0.5)%>%arrange(avg_log2FC)
cat("UP in Ectopic_FC :", paste(head(upML$gene,25),collapse=", "),"\n")
cat("UP in Canonical  :", paste(head(dnML$gene,25),collapse=", "),"\n")

## retain canonical Bergmann identity? check key markers in the DE table
bgmk <- c("Gdf10","Hopx","Fabp7","Slc1a3","Sept4","Gria1","Tnc","Gpr37l1","Aqp4","Timp4")
cat("\nCanonical Bergmann markers in Ectopic_FC vs PCL (log2FC; ~0 = retained):\n")
print(deML[deML$gene %in% bgmk, c("gene","avg_log2FC","p_val_adj")])

## ---- Volcano Ectopic_FC vs Canonical ----
deML$dir<-"ns"; deML$dir[deML$p_val_adj<0.05&deML$avg_log2FC>0.5]<-"Ectopic_FC"
deML$dir[deML$p_val_adj<0.05&deML$avg_log2FC< -0.5]<-"Canonical_PCL"
deML$mlp<- -log10(pmax(deML$p_val_adj,.Machine$double.xmin))
lg<-deML%>%filter(dir!="ns")%>%group_by(dir)%>%slice_max(abs(avg_log2FC),n=18)
v<-ggplot(deML,aes(avg_log2FC,mlp,color=dir))+geom_point(size=1,alpha=.7)+
  scale_color_manual(values=c(Ectopic_FC="#D7263D",Canonical_PCL="#1B998B",ns="grey80"))+
  geom_text_repel(data=lg,aes(label=gene),size=2.6,max.overlaps=40,show.legend=FALSE)+
  geom_vline(xintercept=c(-.5,.5),linetype="dashed",color="grey50")+
  geom_hline(yintercept=-log10(.05),linetype="dashed",color="grey50")+
  labs(title=sprintf("Ectopic-ML Bergmann (n=%d) vs canonical PCL Bergmann (n=%d)",nML,nPCL),
       x="avg log2FC (>0 ectopic-ML)",y="-log10 adj.P")+theme_classic(base_size=11)
ggsave(file.path(FIG,"22_volcano_ectopicFC_vs_canonical.png"),v,width=8.5,height=6.5,dpi=150,bg="white")

## ---- MAP: highlight ectopic Bergmann on the tissue ----
allsp <- readRDS(file.path(BASE,"data","adult_laminar.rds"))$meta
ctx <- allsp[allsp$ctype %in% c("Purkinje","Granule","ML_interneuron","Bergmann","Astro_parenchymal"),]
pmap <- ggplot()+
  geom_point(data=ctx,aes(x,y),color="grey88",size=0.3,stroke=0)+
  geom_point(data=subset(allsp,ctype=="Purkinje"),aes(x,y),color="black",size=0.35,stroke=0)+
  geom_point(data=subset(bg,ectopic=="Canonical_PCL"),aes(x,y),color="#1B998B",size=0.5,stroke=0)+
  geom_point(data=subset(bg,ectopic=="Intermediate"),aes(x,y),color="grey60",size=0.5,stroke=0)+
  geom_point(data=subset(bg,ectopic=="Ectopic_GCL"),aes(x,y),color="#F18701",size=1.4,stroke=0)+
  geom_point(data=subset(bg,ectopic=="Ectopic_FC"),aes(x,y),color="#D7263D",size=1.6,stroke=0)+
  coord_fixed()+theme_void(base_size=10)+
  ggtitle("Ectopic Bergmann: ML-displaced (red), GCL-displaced (orange), canonical PCL (teal)")
ggsave(file.path(FIG,"23_ectopic_bergmann_map.png"),pmap,width=13,height=10,dpi=160,bg="white")

cat("\nDONE_11\n")
