library(tidyverse)
library(ComplexHeatmap)
library(Matrix)
library(customfuncs)

#par_readim("~/desp1/precast/prec_c25q25g3000/deseq_rhth/WT_harmonic/par_im/")
coefs=readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth/WT_harmonic/par_im/coefs.rds")
res_sig=readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth/WT_harmonic/par_im/res_sig.rds") %>% filter(padj<0.01)
setwd("~/desp1/precast/prec_c25q25g3000/figures/figure_2_plots/")
clock_g=readRDS("~/desp/clock_genes.rds")
clock_g=c(clock_g[c(1:3,5)], "Per3", "Hspa5", "Cirbp","Hspa8") 
clust_ord=readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_order.rds")
#############################
# phase heatmap
################

coefs_all=bind_rows(coefs, .id="cluster")
shared_sig=res_sig %>% group_by(gene) %>% tally %>% arrange(desc(n)) 

genes_low_phi_hr_sd=coefs_all %>% group_by(gene) %>% filter(gene %in% shared_sig$gene[1:100]) %>%
  summarise(phi_var=var(phi_hr)) %>% top_n(50, -1*phi_var)




#coefs_all_sig=coefs_all %>% dplyr::filter(gene %in% unique(res_sig$gene))
 #coefs_all_sig=coefs_all %>% dplyr::filter(gene %in% c(shared_sig$gene[1:50], clock_g))
 coefs_all_sig=coefs_all %>% dplyr::filter(gene %in% c(genes_low_phi_hr_sd$gene, clock_g))
 
# coefs_all_sig=coefs_all %>%
#   dplyr::filter(gene %in% c(clock_g,shared_sig$gene[shared_sig$n>=2]))
# coefs_all_sig=coefs_all %>% dplyr::filter(gene %in% c(clock_g))
# coefs_all_sig=coefs_all %>% dplyr::filter(gene %in% shared_sig$gene)
#coefs_all_sig=coefs_all %>% dplyr::filter(gene %in% wt_g10)


coefs_mat=coefs_all_sig %>% dplyr::select(gene, cluster, phi_hr) %>% group_by(gene) %>%
  pivot_wider(names_from = "gene", values_from = "phi_hr", id_cols = "cluster") %>% 
  column_to_rownames("cluster") %>% t()
coefs_mat=coefs_mat %>%na.omit()
col_fun = circlize::colorRamp2(seq(0,24, length.out=4), c(rainbow(2),rev(rainbow(2))))
clock_g_ind=sapply(intersect(clock_g, rownames(coefs_mat)), function(x) {grep(x, rownames(coefs_mat))})

row_ha2= rowAnnotation(genes = anno_mark(at = unlist(clock_g_ind),
                                         side="left",
                                         labels = intersect(clock_g, 
                                                            rownames(coefs_mat))))
gene_ord=coefs_all_sig %>% 
  dplyr::select(gene, cluster, phi_hr) %>% 
  filter(gene %in% rownames(coefs_mat)) %>%
  group_by(gene) %>% na.omit() %>%
  mutate(phi_hr=replace_na(phi_hr, mean(phi_hr))) %>% 
  summarise(mean_phi_hr=ifelse((median(phi_hr)-max(phi_hr))>5,max(phi_hr),median(phi_hr)) ) %>%
  arrange(mean_phi_hr)
colnames(coefs_mat)=gsub("Olfactory Tubercle", "Amygdala",colnames(coefs_mat))
coefs_mat=coefs_mat[gene_ord$gene, clust_ord]

clock_g_ind=sapply(intersect(clock_g, rownames(coefs_mat)), function(x) {grep(x, rownames(coefs_mat))})

row_ha2= rowAnnotation(genes = anno_mark(at = unlist(clock_g_ind),
                                         side="left",
                                         labels = intersect(clock_g, 
                                                            rownames(coefs_mat))))
ht2=Heatmap(coefs_mat, 
            name="ht2",
            show_row_names=F,
            column_names_rot = 75,
            col = col_fun,
            cluster_columns = F, 
            cluster_rows = F,
            left_annotation = row_ha2,
            heatmap_legend_param = list(
              title="Phase",
              at = c(0,12,24))
            )

ht2            

