library(tidyverse)
library(ComplexHeatmap)
library(circlize)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_reg_int2/")

res_kegg = readRDS("enrich_res/enrich_res_kegg_NTG.rds")

res=res_kegg %>% dplyr::filter(p.adjust<0.1) %>% group_by(Description) %>% 
  dplyr::filter(n()>1) %>% mutate(Description=gsub(" - Mus musculus \\(house mouse\\)", "", Description))


str_to_doub=Vectorize(function(x) { eval(parse(text = x))})

res_wide = res %>%
  mutate(ltp=log2(str_to_doub(GeneRatio)/str_to_doub(BgRatio)))%>% 
  dplyr::select(cluster, Description, ltp) %>% 
  group_by(cluster, Description) %>%
  pivot_wider(
    names_from  = Description,
    values_from = ltp,
    values_fill = NA_real_
  )

mat = res_wide %>%
  column_to_rownames("cluster") %>%
  as.matrix()

## color scale
min_score=floor(min(mat, na.rm = TRUE))
max_score = ceiling(max(mat, na.rm = TRUE))

col_fun2 = circlize::colorRamp2(
seq(min_score, max_score, length.out = 100),
viridis::viridis(100)
)
gp6= gpar(fontsize = 6, col = "black",fontfamily="ArialMT")
## draw heatmap
ht = Heatmap(
  mat,
  name = "Fold\nEnrichment",
  col = col_fun2,
  na_col = "grey80",
  cluster_rows = F,
  cluster_columns = F,
  row_names_side = "left",
  column_names_rot = 90,
  column_names_gp = gp6,
  row_names_gp =  gp6,
  heatmap_legend_param = list(
    at = c(min_score, max_score),
    labels = as.character(c(min_score, max_score)),
    labels_gp = gp6,
    title_gp  = gp6
  )
)

draw(ht)

