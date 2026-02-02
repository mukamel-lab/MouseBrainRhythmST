library(Seurat)
library(tidyverse)
library(patchwork)
library(Matrix)
library(ggsci)
library(grid)

sr_outs  <- "~/PATH/TO/spaceranger/042-B/outs"  # download at NCBI accession GSE282203

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

# ---- load Space Ranger output (Visium) ----
seu <- Load10X_Spatial(
  data.dir = sr_outs,
  filename = "filtered_feature_bc_matrix.h5",
  slice = "042-B"
)
DefaultAssay(seu) <- "Spatial"

# ---- genes to plot ----
gns <- c("Drd1", "Trbc2", "Cux2", "C1ql2", "Agt", "Slc17a7", "Neurod6", "Lamp5")
cols_hi <- pals::glasbey(6)[4]

# ---- coordinates table (barcodes) ----
coords <- GetTissueCoordinates(seu) %>%
  rownames_to_column("barcode") %>%
  # keep names similar to your old code
  rename(col = imagecol, row = imagerow)

# ---- expression plots (raw counts from Space Ranger matrix) ----
counts <- GetAssayData(seu, slot = "counts")

plts <- lapply(gns, function(g) {
  if (!g %in% rownames(counts)) {
    # placeholder panel if gene not present
    return(ggplot() + theme_void() + ggtitle(paste0(g, " (not found)")))
  }
  
  df <- coords %>%
    mutate(gene = as.numeric(counts[g, barcode]))
  
  ggplot(df) +
    geom_point(aes(x = col, y = -row, col = gene), size = 0.01) +
    theme_void() +
    scale_color_gradient(
      name = "",
      low = "gray90",
      high = cols_hi,
      limits = c(0, as.numeric(quantile(df$gene, 0.99, na.rm = TRUE))),
      oob = scales::squish
    ) +
    ggtitle(g) +
    coord_fixed() +
    scale_y_reverse() +
    scale_x_reverse() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 12)
    )
})

p <- wrap_plots(plts, ncol = 3)