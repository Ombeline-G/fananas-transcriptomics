## 21 - Per-cell FC verification: one panel per cell with local tissue context
suppressWarnings(suppressMessages({
  library(ggplot2); library(dplyr); library(patchwork)
}))
project_root <- Sys.getenv("BGROOT")
if (!nzchar(project_root)) {
  script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  if (!length(script_arg)) stop("Set BGROOT or run this file with Rscript")
  project_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), "..", "..", "..") , mustWork=TRUE)
}
BASE <- file.path(project_root, "datasets/suyama_2026_xenium")

TAB <- file.path(BASE,"results/tables"); OUT <- file.path(BASE,"results/figures/FC_percell")
dir.create(OUT, showWarnings=FALSE)
sp <- readRDS(file.path(BASE,"data","adult_laminar.rds"))$meta
bg <- read.csv(file.path(TAB,"adult_bergmann_laminar_granule.csv"))
ml <- bg[bg$laminar_class=="FC",] %>% arrange(desc(d_gran))
ml$id <- seq_len(nrow(ml))
cat("FC cells:", nrow(ml), "\n")

xy <- sp; HALF <- 70   # 140 um window
typ_cols <- c(Purkinje="#000000", Granule="grey80", ML_interneuron="#3b8fd1",
              Bergmann="#1B998B", Astro_parenchymal="#74c476", Oligo="#F18701",
              UBC="grey60", Other_neuron="grey75")
panel <- function(r){
  box <- xy[xy$x>=r$x-HALF & xy$x<=r$x+HALF & xy$y>=r$y-HALF & xy$y<=r$y+HALF,]
  box$t <- ifelse(box$ctype %in% names(typ_cols), box$ctype, "other")
  ggplot()+
    geom_point(data=box, aes(x,y,color=t), size=1.5, stroke=0)+
    scale_color_manual(values=c(typ_cols, other="grey88"))+
    annotate("point", x=r$x, y=r$y, size=5.5, shape=21, fill="#D7263D", color="white", stroke=1.1)+
    coord_fixed(xlim=c(r$x-HALF,r$x+HALF), ylim=c(r$y-HALF,r$y+HALF), expand=FALSE)+
    theme_void(base_size=8)+theme(legend.position="none",
      plot.title=element_text(size=8,face="bold",hjust=.5),
      panel.border=element_rect(color="grey70",fill=NA),
      plot.background=element_rect(fill="white",color=NA))+
    ggtitle(sprintf("#%d  d_gran=%.0f", r$id, r$d_gran))
}
## paginate: 5 cols x 5 rows = 25 / page
PER<-25; NC<-5
np <- ceiling(nrow(ml)/PER)
for (pg in 1:np){
  idx <- ((pg-1)*PER+1):min(pg*PER, nrow(ml))
  pls <- lapply(idx, function(i) panel(ml[i,]))
  g <- wrap_plots(pls, ncol=NC) +
    plot_annotation(title=sprintf("FC verification %d/%d  —  cible=rouge | Purkinje=noir, granule=gris, interneurone ML=bleu, Bergmann=teal, astrocyte=vert",pg,np),
      theme=theme(plot.title=element_text(size=11,face="bold")))
  ggsave(file.path(OUT,sprintf("FC_page%02d.png",pg)), g, width=15, height=15.5, dpi=160, bg="white")
}
cat("pages:", np, "\n"); cat("DONE_21\n")
