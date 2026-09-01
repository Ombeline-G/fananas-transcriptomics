## Environment check for Suyama2026 Xenium analysis
cat("R:", R.version.string, "\n")
cat("libPaths:\n"); print(.libPaths())
ip <- rownames(installed.packages())
cat("Total installed packages:", length(ip), "\n")
need <- c("Seurat","SeuratObject","ggplot2","dplyr","tidyr","patchwork",
          "presto","Matrix","data.table","ggrepel","viridis","cowplot",
          "scales","RColorBrewer")
for (p in need) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")
  cat(sprintf("  %-14s %s\n", p, v))
}
