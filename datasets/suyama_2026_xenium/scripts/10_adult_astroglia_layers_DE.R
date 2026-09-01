## 10 - Adult: classify ALL astroglia by layer; ML-resident astroglia?; DE analyses
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

L   <- readRDS(file.path(BASE,"data","adult_laminar.rds")); sp <- L$meta
cb  <- readRDS(file.path(BASE,"data","adult_cerebellum.rds"))
## recreate cell-type annotation (not stored in rds)
map <- c(c15="Purkinje", c0="Granule", c1="Granule", c2="Granule", c7="UBC",
  c9="ML_interneuron", c12="ML_interneuron", c16="ML_interneuron",
  c14="Bergmann", c10="Bergmann",
  c11="Astro_parenchymal", c13="Astro_parenchymal", c3="Astro_parenchymal",
  c6="Oligo", c17="Oligo", c18="OPC", c19="Microglia", c8="Endothelial",
  c23="Mural", c20="Meninges", c21="Choroid", c22="Ependymal",
  c4="Other_neuron", c5="Other_neuron", c24="Other_neuron")
cb$ctype <- unname(map[cb$cl])
## bring transcriptomic scores into sp
mm <- cb@meta.data
for (c in c("score_astrocyte","score_bergmann","score_bg_like","score_igl_astro",
            "score_reactive","cl")) sp[[c]] <- mm[[c]][match(sp$cell, rownames(mm))]
## ctype from the clustered object
sp$ctype <- cb$ctype[match(sp$cell, colnames(cb))]

## ---- generic layer classifier for astroglia ----
layerOf <- function(d_pkj,frac_pkj,frac_gran,frac_ml,frac_oligo){
  out <- rep("Ambiguous", length(d_pkj))
  out[frac_oligo>=0.25] <- "WM"
  pcl <- (d_pkj<=25 | frac_pkj>0) & out=="Ambiguous"; out[pcl] <- "PCL"
  gcl <- frac_gran>=0.35 & out=="Ambiguous"; out[gcl] <- "GCL"
  ml  <- frac_ml>=0.25 & frac_gran<=0.10 & d_pkj>=40 & out=="Ambiguous"; out[ml] <- "ML"
  out
}
sp$layer <- layerOf(sp$d_pkj, sp$frac_pkj, sp$frac_gran, sp$frac_ml, sp$frac_oligo)

astro <- sp[sp$ctype %in% c("Bergmann","Astro_parenchymal"),]
cat("=== astroglia: ctype x layer ===\n"); print(table(astro$ctype, astro$layer))
cat("\n=== mean score_bergmann by ctype x layer ===\n")
print(astro %>% group_by(ctype,layer) %>%
      summarise(n=n(), bergmann=round(mean(score_bergmann,na.rm=T),2),
                astro=round(mean(score_astrocyte,na.rm=T),2),
                igl=round(mean(score_igl_astro,na.rm=T),2), .groups="drop") %>% as.data.frame(), row.names=FALSE)
write.csv(astro, file.path(TAB,"adult_astroglia_layers.csv"), row.names=FALSE)

## ML-resident astroglia of any kind
ml_astro <- astro[astro$layer=="ML",]
cat("\nML-resident astroglia total:", nrow(ml_astro),
    " (Bergmann:", sum(ml_astro$ctype=="Bergmann"),
    ", parenchymal:", sum(ml_astro$ctype=="Astro_parenchymal"),")\n")

## ================= DE analyses =================
DefaultAssay(cb) <- "Xenium"
cb$ctype <- cb$ctype; cb$layer <- sp$layer[match(colnames(cb), sp$cell)]
de_run <- function(obj, g1, g2, file, title){
  d <- FindMarkers(obj, ident.1=g1, ident.2=g2, logfc.threshold=0, min.pct=0.1)
  d$gene <- rownames(d); d <- d %>% arrange(p_val_adj, desc(abs(avg_log2FC)))
  write.csv(d, file.path(TAB,file), row.names=FALSE)
  d
}
volcano <- function(d, posname, negname, file, title){
  d$dir <- "ns"
  d$dir[d$p_val_adj<0.05 & d$avg_log2FC> 0.5] <- posname
  d$dir[d$p_val_adj<0.05 & d$avg_log2FC< -0.5] <- negname
  d$mlp <- -log10(pmax(d$p_val_adj,.Machine$double.xmin))
  labg <- d %>% filter(dir!="ns") %>% group_by(dir) %>% slice_max(abs(avg_log2FC),n=15)
  cols <- setNames(c("#D7263D","#1B998B","grey80"), c(posname,negname,"ns"))
  p <- ggplot(d, aes(avg_log2FC, mlp, color=dir)) + geom_point(size=1,alpha=.7)+
    scale_color_manual(values=cols)+
    geom_text_repel(data=labg,aes(label=gene),size=2.6,max.overlaps=30,show.legend=FALSE)+
    geom_vline(xintercept=c(-.5,.5),linetype="dashed",color="grey50")+
    geom_hline(yintercept=-log10(.05),linetype="dashed",color="grey50")+
    labs(title=title,x=sprintf("avg log2FC (>0 %s)",posname),y="-log10 adj.P")+
    theme_classic(base_size=11)
  ggsave(file.path(FIG,file), p, width=8,height=6.5,dpi=150,bg="white"); invisible(p)
}

## (A) within-PCL Bergmann subtypes: c10 vs c14
Idents(cb) <- "cl"
deA <- de_run(cb, "c10","c14","DE_adult_BGc10_vs_BGc14.csv","BG c10 vs c14")
cat(sprintf("\n(A) c10 vs c14: %d sig (|lfc|>0.5,padj<0.05)\n",
            sum(deA$p_val_adj<0.05 & abs(deA$avg_log2FC)>0.5)))
volcano(deA,"BG_c10","BG_c14","19_volcano_BGc10_vs_c14.png","Bergmann c10 vs c14 (both PCL)")

## (B) PCL Bergmann vs parenchymal astrocyte (reference contrast)
cb$grpB <- ifelse(cb$ctype=="Bergmann","Bergmann",
            ifelse(cb$ctype=="Astro_parenchymal","Astro_parench",NA))
Idents(cb) <- "grpB"
deB <- de_run(subset(cb,subset=!is.na(grpB)),"Bergmann","Astro_parench",
              "DE_adult_Bergmann_vs_parenchymal.csv","Bergmann vs parenchymal astro")
volcano(deB,"Bergmann","Astro_parench","20_volcano_Bergmann_vs_parenchymal.png",
        "Adult Bergmann vs parenchymal astrocyte")

## (C) ML-resident astroglia vs PCL Bergmann -- only if enough ML cells
n_ml <- nrow(ml_astro)
if (n_ml>=30){
  cb$grpC <- NA
  cb$grpC[colnames(cb) %in% ml_astro$cell] <- "ML_astroglia"
  cb$grpC[colnames(cb) %in% astro$cell[astro$ctype=="Bergmann"&astro$layer=="PCL"]] <- "PCL_Bergmann"
  Idents(cb) <- "grpC"
  deC <- de_run(subset(cb,subset=!is.na(grpC)),"ML_astroglia","PCL_Bergmann",
                "DE_adult_MLastroglia_vs_PCLBergmann.csv","ML-astroglia vs PCL-Bergmann")
  volcano(deC,"ML_astroglia","PCL_Bergmann","21_volcano_MLastroglia_vs_PCLBergmann.png",
          "ML-resident astroglia vs PCL Bergmann")
  cat(sprintf("(C) ML-astroglia vs PCL-Bergmann: %d sig\n",
              sum(deC$p_val_adj<0.05 & abs(deC$avg_log2FC)>0.5)))
} else cat(sprintf("\n(C) SKIPPED: only %d ML-resident astroglia (too few for DE)\n", n_ml))

cat("\nDONE_10\n")
