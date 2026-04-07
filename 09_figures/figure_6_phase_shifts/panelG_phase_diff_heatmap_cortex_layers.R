library(tidyverse)
library(ComplexHeatmap)
library(patchwork)
library(scico)

cluster_color <- readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
temp_col=cluster_color



x=vroom::vroom("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth_int/for_share/full/results.csv.gz") %>%
  dplyr::filter(padj<0.1,cluster %in% names(cluster_color)[1:5])
phase_relamp <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhyth2/old_young_sep/phase_relamp.rds")

y=phase_relamp %>% dplyr::filter(gene %in% x$gene, 
                                 cluster %in% names(cluster_color)[1:5])

z=y %>% dplyr::select(genotype, cluster,gene,rel_amp,age, phi, amp, phi_hr) %>%
  pivot_wider(id_cols =c( gene, age, cluster),
              names_from = genotype, values_from = c(rel_amp, amp, phi_hr, phi) )

z$amp_dif=2*(z$amp_APP-z$amp_WT)
z$rel_amp_dif=(z$rel_amp_APP-z$rel_amp_WT)
z$rel_amp_l2fc=log2(z$rel_amp_APP/z$rel_amp_WT)
z$phi_dif=(12/pi)*atan2(sin(z$phi_APP-z$phi_WT),cos(z$phi_APP-z$phi_WT))
###################################################

to_plt_ral2fc=z %>% 
  dplyr::select(gene, cluster, age,rel_amp_l2fc ) %>% 
  pivot_wider(id_cols = gene, names_from = c(cluster, age), values_from =rel_amp_l2fc )

clord=paste0(rep(names(cluster_color)[1:5], times=2), rep(c("_young", "_old"), each=5))

to_plt_ral2fc=to_plt_ral2fc %>% column_to_rownames("gene")

to_plt_ral2fc=as.matrix(to_plt_ral2fc)



to_plt_ral2fc=to_plt_ral2fc[,clord]
###############################################
anno=data.frame(ROI=str_split_fixed(colnames(to_plt_ral2fc), "_",2)[,1],
                age=str_split_fixed(colnames(to_plt_ral2fc), "_",2)[,2])
anno$age=plyr::mapvalues(anno$age, c("old", "young"), c("14 months", "7 months"))
column_ha = HeatmapAnnotation(ROI=anno$ROI, Age =anno$age,
                              col=list(ROI=cluster_color[1:5],
                                       Age=c(RColorBrewer::brewer.pal(12, "Paired")[11],
                                             RColorBrewer::brewer.pal(3, "BrBG")[1]) %>%
                                         setNames(paste0(c(7, 14), " months"))),
                              annotation_legend_param = list(ROI = list(
                                at = unique(anno$ROI)
                              )))
################################################
mm=quantile(abs(to_plt_ral2fc),0.999,na.rm=T)

Heatmap(to_plt_ral2fc,name="ht1",
        cluster_columns = F, 
        top_annotation = column_ha,
        col=circlize::colorRamp2(seq(-mm, mm, length = 200), 
                                 scico(200,palette='bam')),
        show_column_names = F,
        heatmap_legend_param = list(
          title = "Log2FC\nRel. Amp"))
decorate_heatmap_body("ht1", 
                      {grid.lines(c(.5,.5),
                                  c(0, 1),
                                  gp = gpar(lty = 1, lwd = 4))
                      }, slice=1)



############################
to_plt_ral2fc=z %>% 
  dplyr::select(gene, cluster, age,phi_dif ) %>% 
  pivot_wider(id_cols = gene, names_from = c(cluster, age), values_from =phi_dif )


to_plt_ral2fc=to_plt_ral2fc %>% column_to_rownames("gene")

to_plt_ral2fc=as.matrix(to_plt_ral2fc)



to_plt_ral2fc=to_plt_ral2fc[,clord]


mm <- 12
k  <- 0.80  # how quickly it fades to black (0.7–0.9 is nice)
brks <- c(-mm, -mm*k, 0,  mm*k,  mm)
cols <- c("#000000",
          grDevices::hcl(h = 230, c = 70, l = 45),  # deep blue
          "#FFFFFF",
          grDevices::hcl(h =  10, c = 70, l = 45),  # maroon
          "#000000")


genes <- c(
  "Itpr1","Adrb1","Nr4a1","Dusp6","Csrnp1","Trib1","Cd200","Medag","Gemin6","Prdm2",
  "Fzd4","Smim1","Eif5","Bace2","Ccnd2","Igf1r","Map3k5","Rfesd","Cpsf4","Mier2",
  "Cd4","Clic5","Ppargc1a","Pik3r4","Med1","Strip2","Flrt1","Klf10","Rbm24","Sp2",
  "Dcaf4","Midn"
)
Heatmap(to_plt_ral2fc[genes,], name="ht1",
        cluster_columns = F,
        cluster_rows = F,
        
        col=circlize::colorRamp2(seq(-mm, mm, length =length(cols)), 
                                 cols),
        #scico(361,palette='vik') ),

        top_annotation = column_ha,
        show_column_names = F,
        heatmap_legend_param = list(
          title = "Phase Difference",
         at = c(-12, -6, 0, 6, 12),
         labels = c(-12, -6, 0, 6, 12))
)


decorate_heatmap_body("ht1", 
                      {grid.lines(c(.5,.5),
                                  c(0, 1),
                                  gp = gpar(lty = 1, lwd = 4))
                      }, slice=1)

