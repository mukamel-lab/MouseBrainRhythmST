library(tidyverse)
library(glue)
library(ComplexHeatmap)
library(magrittr)
library(scico)
library(viridis)

res_sig <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/WT_harmonic/par_im/res.rds") %>% 
  filter(padj<0.01)
res_sig=split(res_sig, res_sig$cluster)
res_sig=lapply(res_sig, function(y) {y$gene})
names(res_sig)=gsub("Olfactory Tubercle", "Amygdala",names(res_sig))

cluster_order = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_order.rds")

jci=matrix(0, nrow=23, ncol=23)
rownames(jci)=cluster_order
colnames(jci)=cluster_order


for(i in 1:23) {
  for(j in 1:23) {
    jci[i,j]=length(intersect(res_sig[[cluster_order[i]]], res_sig[[cluster_order[j]]]))/
      length(unique(union(res_sig[[cluster_order[i]]], res_sig[[cluster_order[j]]])))
  }
}

for(i in 1:23){
  jci[i,i]=NA
}
col_fun = circlize::colorRamp2(seq(0,max(jci, na.rm = T), length.out=1000),viridis(1000) )

Heatmap(jci[c(1:11),c(1:11)], cluster_rows = F, cluster_columns = F, na_col = "white",
        col=col_fun, column_names_side = "bottom",row_names_side = "left",
        rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(i > j) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          }    })

rownames(jci)=gsub(cluster_order[6], "CTXsp",
                   gsub("Dentate Gyrus", 
                   "DG",
                   gsub("Cortex Layer ", 
                        "L", 
                        rownames(jci))), fixed=T)
colnames(jci)=gsub(cluster_order[6], "CTXsp",
                   gsub("Dentate Gyrus", 
                   "DG",
                   gsub("Cortex Layer ",
                        "L", 
                        colnames(jci))), fixed = T)

Heatmap(jci[c(1:11),c(1:11)], cluster_rows = F, cluster_columns = F, na_col = "white",col=col_fun,column_names_rot = 55)

