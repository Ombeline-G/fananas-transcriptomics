## 22 - Prepare per-cell region/type data + FC cards for the interactive swipe reviewer
suppressWarnings(suppressMessages({
  library(dplyr); library(RANN)
}))
set.seed(0)
project_root <- Sys.getenv("BGROOT")
if (!nzchar(project_root)) {
  script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  if (!length(script_arg)) stop("Set BGROOT or run this file with Rscript")
  project_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), "..", "..", "..") , mustWork=TRUE)
}
BASE <- file.path(project_root, "datasets/suyama_2026_xenium")

TAB <- file.path(BASE,"results/tables"); OUT <- file.path(BASE,"results/swipe_validation")
dir.create(OUT, showWarnings=FALSE)
sp <- readRDS(file.path(BASE,"data","adult_laminar.rds"))$meta
bg <- read.csv(file.path(TAB,"adult_bergmann_laminar_granule.csv"))

## ---- per-cell layer REGION (0=ML,1=GCL,2=WM) by kNN vote among anatomical anchors ----
anchor_reg <- c(ML_interneuron=0, Granule=1, Oligo=2, OPC=2)
anc <- sp[sp$ctype %in% names(anchor_reg),]
anc$reg <- anchor_reg[anc$ctype]
kn <- RANN::nn2(as.matrix(anc[,c("x","y")]), as.matrix(sp[,c("x","y")]), k=9)
regmat <- matrix(anc$reg[kn$nn.idx], ncol=9)
sp$region <- apply(regmat, 1, function(v){ t<-tabulate(v+1,3); which.max(t)-1 })

## ---- type: 1=Purkinje (PCL line), 2=BG_PCL, 0=other ----
bgpcl <- bg$cell[bg$laminar_class=="BG_PCL"]
sp$type <- ifelse(sp$ctype=="Purkinje", 1L, ifelse(sp$cell %in% bgpcl, 2L, 0L))

## ---- downsample for display (keep Purkinje, BG_PCL, ML interneurons, candidates; sample rest) ----
cand <- bg$cell[bg$laminar_class=="FC"]
keepmask <- sp$type>0 | sp$ctype=="ML_interneuron" | sp$cell %in% cand
cap <- 11000
rest <- which(!keepmask); nrest <- max(0, cap - sum(keepmask))
if (length(rest) > nrest) rest <- sample(rest, nrest)
sel <- sort(unique(c(which(keepmask), rest)))
disp <- sp[sel, c("x","y","region","type")]
write.csv(disp, file.path(OUT,"swipe_cells.csv"), row.names=FALSE)
cat("display cells:", nrow(disp), " (purk:",sum(sp$type==1),"BG_PCL:",sum(sp$type==2),"ML-int:",sum(sp$ctype=="ML_interneuron"),")\n")

## ---- cards: 145 FC, sorted by d_gran desc; add d_pkj + Purkinje support ----
ml <- bg[bg$laminar_class=="FC",]
pk_xy <- as.matrix(sp[sp$ctype=="Purkinje", c("x","y")])
dpk <- RANN::nn2(pk_xy, as.matrix(ml[,c("x","y")]), k=1)$nn.dists[,1]
## Purkinje within 40 um (support)
within <- RANN::nn2(pk_xy, as.matrix(ml[,c("x","y")]), k=15)$nn.dists
ml$support <- rowSums(within <= 40)
ml$d_pkj1 <- round(dpk,0)
cards <- ml %>% arrange(desc(d_gran)) %>%
  transmute(id=cell, cx=round(x,1), cy=round(y,1),
            d_gran=round(d_gran,0), d_pkj=d_pkj1, support=support)
write.csv(cards, file.path(OUT,"swipe_cards.csv"), row.names=FALSE)
cat("FC cards:", nrow(cards), "\n")
cat("section bbox: x[",round(min(sp$x)),",",round(max(sp$x)),"] y[",round(min(sp$y)),",",round(max(sp$y)),"]\n")
cat("DONE_22\n")
