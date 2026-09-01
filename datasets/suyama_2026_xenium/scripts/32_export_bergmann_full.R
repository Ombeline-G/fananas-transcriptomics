suppressWarnings(suppressMessages({
  library(Seurat); library(SeuratObject)}))
project_root <- Sys.getenv("BGROOT")
if (!nzchar(project_root)) {
  script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  if (!length(script_arg)) stop("Set BGROOT or run this file with Rscript")
  project_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), "..", "..", "..") , mustWork=TRUE)
}
BASE <- file.path(project_root, "datasets/suyama_2026_xenium")
TAB <- file.path(BASE,"results","tables")
dec <- read.csv(file.path(BASE,"results","fc_decisions_suyama.csv")); ml <- as.character(dec$cell_label[dec$decision=="yes"])
lam <- read.csv(file.path(TAB,"adult_bergmann_laminar_purkinje.csv")); pcl <- setdiff(as.character(lam$cell[lam$bgclass=="BG_PCL"]), ml)
cb <- readRDS(file.path(BASE,"data","adult_cerebellum.rds")); DefaultAssay(cb)<-"Xenium"
cells <- intersect(c(ml,pcl), colnames(cb)); cat("cells:",length(cells),"\n")
M <- t(as.matrix(GetAssayData(cb,assay="Xenium",layer="data")[,cells]))
df <- data.frame(cell=cells, group=ifelse(cells %in% ml,"FC","BG_PCL"), M, check.names=FALSE)
write.csv(df, file.path(TAB,"bergmann_full_expr.csv"), row.names=FALSE)
cat("wrote bergmann_full_expr.csv",nrow(df),"x",ncol(M),"\nDONE\n")
