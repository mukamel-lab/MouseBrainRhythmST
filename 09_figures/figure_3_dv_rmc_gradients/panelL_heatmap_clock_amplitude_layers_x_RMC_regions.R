library(tidyverse)
library(ComplexHeatmap)
library(Matrix)
library(customfuncs)

setwd("/home/agelber/desp1/precast/precast_final_with_ros_caud/")


res = readRDS("analysis/deseq_rhyth/WT/res.rds")  
coefs = readRDS("analysis/deseq_rhyth/WT/coefs.rds") %>%
  bind_rows(.id = "cluster")

res=left_join( coefs,res) %>%
  mutate(
    phi = atan2(t_s, t_c) %% (2 * pi),
    amp = sqrt(t_c ^ 2 + t_s ^ 2),
    rel_amp = amp / log2(baseMean + 1),
    phi_hr = (12 / pi) * phi
  ) %>% 
  dplyr::select(gene, padj,cluster, amp, rel_amp,phi_hr) %>% 
  dplyr::mutate(cluster=gsub("23", "2\\/3", cluster))%>%
  na.omit()

clock_g = c("Arntl",
            "Cry1" ,
            "Cry2",
            "Dbp" ,
            "Per1" ,
            "Per2",
            "Per3",
            "Nr1d1",
            "Nr1d2",
            "Bhlhe41")

clust_ord=readRDS("agg_c.rds") %>% names()
clust_ord=clust_ord[c(3,6,9,12,15,18,2,5,8,11,14,17,1,4,7,10,13,16)]

#############################
# rel_amp heatmap
################




coefs_mat=res %>% dplyr::filter(padj<0.1, cluster %in% clust_ord) %>%
  dplyr::select(gene, cluster,  amp) %>% 
  # group_by(gene) %>%
  # mutate(rel_amp=(rel_amp-mean(rel_amp))/sd(rel_amp)) %>%
  pivot_wider(names_from = "gene", 
              values_from = "amp",
              id_cols = "cluster") %>% 
  column_to_rownames("cluster") %>% t()


coefs_mat=coefs_mat[clock_g, clust_ord]


regions=c(rep("Rostral",6 ),rep("Medial",6 ), rep("Caudal",6 ))

col_ha = HeatmapAnnotation(Region = regions,
                           col=list(Region=
                                      c(
                                        ggsci::pal_d3()(3)[1:3]) %>%
                                          setNames(c("Rostral", 
                                                     "Medial", 
                                                     "Caudal"))
                                      ),
                           annotation_legend_param = 
                             list(Region=list(
                                 at=c(
                                 "Rostral",
                                 "Medial",
                                 "Caudal"))),
                           height = unit(1, "mm") )



col_fun = circlize::colorRamp2(c(0,1), c("white", "red"))

short_names = readRDS("short_names.rds")
short_namesdf=data.frame(full=short_names,short=names(short_names))%>%
separate(full, into = c("clust", "region"), sep=" - ", remove = F)%>%
separate(short, into = c("short", "R"), sep="_")%>%
mutate(short=gsub("23","2/3", short))

colnames(coefs_mat)=plyr::mapvalues(colnames(coefs_mat), short_namesdf$full, short_namesdf$short)
col_ha@anno_list = lapply(col_ha@anno_list, function(ann) {
  ann@show_legend = FALSE
  ann
})
col_ha@anno_list$Region@label=""

ht2 = Heatmap(
  coefs_mat,
  name = "ht2",
  show_row_names = TRUE,
  column_names_rot = 90,
  top_annotation = col_ha,
  row_names_side = "left",
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  col = col_fun,
  
  column_names_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
  row_names_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
  
  heatmap_legend_param = list(
    title = "Log2\nAmplitude",
    legend_direction = "vertical",
    title_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
    labels_gp = gpar(fontsize = 6, fontfamily = "ArialMT")
  )
)



ht2

decorate_heatmap_body("ht2", {
  grid.lines(c(1/3, 1/3), c(0, 1), gp = gpar(lty = 1, lwd = 2))
}, slice = 1)

decorate_heatmap_body("ht2", {
  grid.lines(c(2/3, 2/3), c(0, 1), gp = gpar(lty = 1, lwd = 2))
}, slice = 1)

saveRDS(ht2, "~/desp1/precast/precast_final_with_ros_caud/figures/figure2/heatmap_clock_rel_amp_WT.rds")
       

