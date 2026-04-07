library(tidyverse)
library(ComplexHeatmap)

setwd("~/desp1/precast/precast_final_with_ros_caud/analysis/deseq_region_interaction/main_effect_res/gsea")

res=list.files(recursive = T)

res=grep("enrichment_results", res, value = T)

names(res)=gsub("/", "_",sub("/P.*$", "", res))

res=lapply(res, vroom::vroom)

res=bind_rows(res, .id="grp") %>%
  separate(grp, into = c( "comparison","cluster"))
res=res %>% mutate(up_in=ifelse(normalizedEnrichmentScore>0,
                                substr(comparison, 1, 1),
                                substr(comparison,2,2))) 
res_sig=res  %>% group_by(description) %>% dplyr::filter(FDR<0.2, n()>=6, description!="Ribosome")

res_sig_wide=res %>% ungroup() %>% dplyr::filter(cluster!="L4", description %in% r_ord) %>%
  pivot_wider(id_cols = "description", names_from = c("comparison", "cluster"),
              values_from = "normalizedEnrichmentScore", values_fill = 0)


res_sig_wide=res_sig_wide %>% column_to_rownames("description")
r_ord=c("Neuroactive ligand-receptor interaction", "Cytokine-cytokine receptor interaction",
        "Polycomb repressive complex","RNA degradation", "Circadian entrainment"
        
)

res_sig_wide=res_sig_wide[r_ord,sort(colnames(res_sig_wide))]
col_ha=str_split_fixed(colnames(res_sig_wide), "_",2)
col_ha=as.matrix(col_ha)  %>% as.data.frame() %>% t()
col_ha[,1]=plyr::mapvalues(col_ha[,1], c("RM","MC","RC"), c("R-M","M-C","R-C"))

col_ha = HeatmapAnnotation(Comp = col_ha[,1],clust=col_ha[,2],
                           col=list(Comp=
                                      c("R-M"="#6E984D",
                                        "M-C"="#505CA9",
                                        "R-C"="#CC2F27"
                                      )
                           ),
                           annotation_legend_param = 
                             list(Comp=list(
                               at=c("R-M","M-C","R-C"))),
                           height = unit(2, "mm") )

Heatmap(res_sig_wide, 
        show_row_names = T,
        bottom_annotation = col_ha,
        row_names_side = "left",
        show_column_names = F,
        cluster_rows = F,cluster_columns = F)
