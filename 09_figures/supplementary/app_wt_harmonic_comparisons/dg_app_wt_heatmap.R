library(tidyverse)
library(ComplexHeatmap)

###############################
# sample order
################################
data_dir="/cndd2/agelber/hal/qc_aligned"

metadata = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) %>%
  dplyr::select(sample,time, age, genotype) %>% distinct() %>%
  mutate(age=factor(age, levels=c("7 months", "14 months")),
         time=factor(time, levels = paste0("ZT", 6*0:3)))%>%
  arrange(genotype,age,time)


wt_ncounts=readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/WT_harmonic/par_im/ncounts.rds") %>% 
  .[["Dentate Gyrus-sg"]]
app_ncounts=readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/APP_harmonic/par_im/ncounts.rds") %>% 
  .[["Dentate Gyrus-sg"]]

#choose genes
setwd("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/")

gene_df=readRDS("~/desp1/precast/prec_c25q25g3000/figures/fig3_app_wt_harm/both_sig_dgsg.rds")

gene_df=split(gene_df, gene_df$sig)
gene_df=lapply(gene_df, function(x) { x %>% arrange(wt_phi_hr)})

gene_df$`APP23-TG`=gene_df$`APP23-TG` %>% arrange(app_phi_hr)



gene_to_plot=c(gene_df$Both$gene, gene_df$NTG$gene ,gene_df$`APP23-TG`$gene)





############################
# expression dfs
###########################

wt= wt_ncounts %>%
  dplyr::filter(gene %in% gene_to_plot) %>% 
  group_by(gene) %>%
  dplyr::select(gene,sample,age,l2expr)%>%
  mutate(geno="NTG")
 
app=app_ncounts %>%
  dplyr::filter(gene %in% gene_to_plot) %>% 
  group_by(gene) %>%
  dplyr::select(gene,sample,age, l2expr) %>%
  mutate(geno="APP23-TG")

wt_app_df=rbind(wt,app)
wt_app_df=wt_app_df %>% ungroup() %>%
  group_by(gene, age) %>%
  mutate(z=(l2expr-mean(l2expr))/sd(l2expr))%>%
  ungroup() %>%
  pivot_wider(names_from = "gene", 
              values_from = "z",
              id_cols = c("sample", "geno")) %>%
  as.data.frame() %>%
  mutate(sample=factor(sample, levels=metadata$sample))%>%  
  arrange(desc(geno),sample)

wt_app_mat=as.matrix(wt_app_df %>% mutate(geno=NULL, sample=NULL))
rownames(wt_app_mat)=wt_app_df$sample
wt_app_mat=t(wt_app_mat)
wt_app_mat=wt_app_mat[gene_to_plot,]
wt_app=wt_app_mat %>% na.omit()


anno=data.frame(sample=colnames(wt_app),
                geno=wt_app_df$geno)

anno=left_join(anno, metadata)
gene_anno=c(rep("Both", sum(gene_df$Both$gene %in% rownames(wt_app))),
            rep("NTG", sum(gene_df$NTG$gene %in% rownames(wt_app))),
            rep("APP23-TG", sum(gene_df$`APP23-TG`$gene %in% rownames(wt_app))))

column_ha = HeatmapAnnotation(Genotype=anno$geno, Age =anno$age, Time =anno$time,
                              col=list(Genotype=c("#0072B5FF","#BC3C29FF") %>%
                                         setNames(c("NTG",
                                                    "APP23-TG")),
                                       Time=RColorBrewer::brewer.pal(4, "Set2") %>%
                                         setNames(paste0("ZT", 6*0:3)),
                                       Age=c(RColorBrewer::brewer.pal(12, "Paired")[11],
                                             RColorBrewer::brewer.pal(3, "BrBG")[1]) %>%
                                         setNames(paste0(c(7, 14), " months"))),
                              annotation_legend_param = list(Genotype = list(
                                at = unique(anno$geno)
                              )))

gene_df2=bind_rows(gene_df, .id="sig")
rownames(gene_df2)=gene_df2$gene
gene_df2=gene_df2[rownames(wt_app),]

mean_expr=gene_df2[,c("app_baseMean", "wt_baseMean")]
mean_expr$app_baseMean=log2(mean_expr$app_baseMean+1)
mean_expr$wt_baseMean=log2(mean_expr$wt_baseMean+1)

phase_mat=gene_df2[,c("app_phi_hr", "wt_phi_hr")]

gene_df2$relamp_lfc=log2(gene_df2$app_rel_amp/gene_df2$wt_rel_amp)

row_ha1= rowAnnotation(Gene=gene_anno, 
                       col=list(Gene=c(RColorBrewer::brewer.pal(3, "Dark2")[1],"#0072B5FF","#BC3C29FF") %>%
                                  as.character %>%
                                  setNames(c("Both", "NTG", "APP23-TG"))),
                       annotation_legend_param = list(Gene = list(
                         at = c("Both", "NTG", "APP23-TG"))
                       ))

row_ha2 = rowAnnotation(
  Rel_Amp_L2FC =
    anno_points(gene_df2$relamp_lfc,
                gp = gpar(
                  col = ifelse(gene_df2$relamp_lfc > 0,
                               "orange", "gray60")
                )), width =unit(3, "cm")
  
  # Mean_Expression = anno_points(as.matrix(mean_expr),
  #                               gp = gpar(col =
  #                                           c("#BC3C29FF",
  #                                                      "#0072B5FF"))),
    # Phase =  anno_points(as.matrix(phase_mat),
    #             gp = gpar(col =
    #                         c("#BC3C29FF",
    #                                      "#0072B5FF"))),
    #                                    width =unit(6, "cm")
    # 
)

Heatmap(wt_app,  show_row_names=F,
        name="ht1",
        show_column_names = F,
        cluster_rows = F, 
        cluster_columns = F, 
        row_order = 1:nrow(wt_app),
        column_order = 1:ncol(wt_app),
        top_annotation = column_ha,
        left_annotation = row_ha1,
        right_annotation = row_ha2,
        heatmap_legend_param = list(
          title = "Z-scored\nLog2 Expr"))

decorate_heatmap_body("ht1", 
                      {grid.lines(c((length(unique(wt$sample)))/ncol(wt_app), (length(unique(wt$sample)))/ncol(wt_app)),
                                  c(0, 1),
                                  gp = gpar(lty = 2, lwd = 4))
                      }, slice=1)

gsize1=sum(gene_df$`APP23-TG`$gene %in% rownames(wt_app))
gsize2=gsize1+sum(gene_df$NTG$gene %in% rownames(wt_app))

decorate_heatmap_body("ht1", 
                      {grid.lines(c(0, 1),
                                  c(gsize1/nrow(wt_app), 
                                    gsize1/nrow(wt_app)),
                                  gp = gpar(lty = 2, lwd = 4))
                      }, slice=1)


decorate_heatmap_body("ht1", 
                      {grid.lines(c(0, 1),
                                  c(gsize2/nrow(wt_app), 
                                    gsize2/nrow(wt_app)),
                                  gp = gpar(lty = 2, lwd = 4))
                      }, slice=1)



time_lines=as.data.frame(metadata)
rownames(time_lines)=metadata$sample
time_lines=time_lines[colnames(wt_app),]
prev="ZT0"
for(i in 1:nrow(time_lines)){
  if(time_lines$time[i]!=prev) {
    decorate_heatmap_body("ht1", 
                          {grid.lines(c((i-1)/ncol(wt_app), 
                                        (i-1)/ncol(wt_app)),
                                      c(0, 1),
                                      gp = gpar(lty = 1, lwd = 0.8, col="gray45"))
                          }, slice=1)
    
    
  }
  
  prev=time_lines$time[i]
  
}

