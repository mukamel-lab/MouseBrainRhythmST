library(tidyverse)
library(ComplexHeatmap)

###############################
# sample order
################################
data_dir="/cndd2/agelber/hal/qc_aligned"

metadata = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir), genotype=="WT") %>%
  dplyr::select(sample,time, age) %>% distinct() %>%
  mutate(age=factor(age, levels=c("7 months", "14 months")),
         time=factor(time, levels = paste0("ZT", 6*0:3)))%>%
  arrange(age,time)


#choose genes
setwd("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/WT_harmonic/")
ncounts=readRDS("par_im/ncounts.rds")

res_sig = readRDS("par_im/res_sig.rds") 
coefs=readRDS("par_im/coefs.rds") %>% bind_rows(.id="cluster") %>%  dplyr::select(gene, cluster, amp, phi_hr)
res_sig=left_join(res_sig, coefs) %>%
  mutate(rel_amp=amp/log2(baseMean+1)) %>% na.omit()



#choose top shared

res_l23_dgsg=res_sig %>% 
  dplyr::filter(cluster %in% c("Dentate Gyrus-sg", "Cortex Layer 2/3"))

shared=res_l23_dgsg %>% group_by(gene) %>% 
  dplyr::filter(n()>=2) %>% mutate(rel_amp_dif=abs(rel_amp[1]-rel_amp[2])) %>%
  dplyr::filter(rel_amp_dif<quantile(.$rel_amp_dif,.8))%>%
  mutate(phi_hr_dif=(phi_hr[1]-phi_hr[2]) %% 24) %>%
  dplyr::filter(phi_hr_dif<quantile(.$phi_hr_dif,.7))
# mutate(mean_rel_amp=mean(rel_amp)) %>% 
# ungroup() %>% top_n(150, mean_rel_amp)

shared_gn=unique(shared$gene)

# choose top dg

sig_dgsg_only=setdiff(res_l23_dgsg$gene[res_l23_dgsg$cluster!="Cortex Layer 2/3"], res_l23_dgsg$gene[res_l23_dgsg$cluster=="Cortex Layer 2/3"])
res=readRDS("par_im/res.rds") %>%
  dplyr::filter(cluster %in% c("Dentate Gyrus-sg", "Cortex Layer 2/3"))%>%
  left_join(coefs) %>% mutate(rel_amp=amp/log2(baseMean+1))

sig_dgsg=res %>% dplyr::filter(gene %in% sig_dgsg_only) %>% group_by(gene) %>%
  mutate(rel_amp_dif=abs(rel_amp[1]-rel_amp[2]), max_rel_amp=max(rel_amp),
         prod_dif_max=rel_amp_dif*max_rel_amp,
         pval_val_dif=padj[2]/padj[2]) %>%
  ungroup() %>% 
  top_n(200, max_rel_amp)

# choose top dg



sig_ctx_only=setdiff(res_l23_dgsg$gene[res_l23_dgsg$cluster=="Cortex Layer 2/3"], res_l23_dgsg$gene[res_l23_dgsg$cluster!="Cortex Layer 2/3"])
res=readRDS("par_im/res.rds") %>%
  dplyr::filter(cluster %in% c("Dentate Gyrus-sg", "Cortex Layer 2/3"))%>%
  left_join(coefs) %>% mutate(rel_amp=amp/log2(baseMean+1))

sig_ctl23=res %>% dplyr::filter(gene %in% sig_ctx_only) %>% group_by(gene) %>%
  mutate(rel_amp_dif=abs(rel_amp[1]-rel_amp[2]), max_rel_amp=max(rel_amp),
         prod_dif_max=rel_amp_dif*max_rel_amp,
         pval_val_dif=padj[2]/padj[2]) %>%
  ungroup() %>% 
  top_n(200, max_rel_amp)
############################

gene_ord_shared=readRDS("par_im/coefs.rds") %>% 
  bind_rows(.id="cluster") %>% 
  dplyr::select(gene, cluster, phi_hr) %>%
  dplyr::filter(cluster %in% c("Dentate Gyrus-sg", "Cortex Layer 2/3"),
                gene %in% unique(shared$gene)) %>%
  group_by(gene)%>%
  summarise(mean_phi_hr=mean(phi_hr)) %>%
  arrange(desc(mean_phi_hr))

gene_ord_cort=readRDS("par_im/coefs.rds") %>% 
  bind_rows(.id="cluster") %>% 
  dplyr::select(gene, cluster, phi_hr) %>%
  dplyr::filter(cluster %in% c("Cortex Layer 2/3"),
                gene %in% unique(sig_ctl23$gene)) %>%
  group_by(gene)%>%
  summarise(mean_phi_hr=mean(phi_hr)) %>%
  arrange(desc(mean_phi_hr))


gene_ord_dg=readRDS("par_im/coefs.rds") %>% 
  bind_rows(.id="cluster") %>% 
  dplyr::select(gene, cluster, phi_hr) %>%
  dplyr::filter(cluster %in% c("Dentate Gyrus-sg"),
                gene %in% unique(sig_dgsg$gene)) %>%
  group_by(gene)%>%
  summarise(mean_phi_hr=mean(phi_hr)) %>%
  arrange(desc(mean_phi_hr))

gene_to_plot=c(gene_ord_shared$gene,  gene_ord_cort$gene, gene_ord_dg$gene)





############################
# expression dfs
###########################

l23=ncounts$`Dentate Gyrus-sg` %>% 
  dplyr::filter(gene %in% gene_to_plot) %>% 
  group_by(gene) %>%
  dplyr::select(gene,sample,age,l2expr)%>%
  mutate(cluster="Dentate Gyrus-sg")

dg=ncounts$`Cortex Layer 2/3` %>% 
  dplyr::filter(gene %in% gene_to_plot) %>% 
  group_by(gene) %>%
  dplyr::select(gene,sample,age, l2expr) %>%
  mutate(cluster="Cortex Layer 2/3")

dg_l23_df=rbind(l23,dg)
dg_l23_df=dg_l23_df %>% ungroup() %>%
  group_by(gene, age, cluster) %>%
  mutate(z=(l2expr-mean(l2expr))/sd(l2expr))%>%
  ungroup() %>%
  pivot_wider(names_from = "gene", 
            values_from = "z",
            id_cols = c("sample", "cluster")) %>%
  as.data.frame() %>%
  mutate(sample=factor(sample, levels=metadata$sample))%>%  
  arrange(desc(cluster),sample)

dg_l23_mat=as.matrix(dg_l23_df %>% mutate(cluster=NULL, sample=NULL))
rownames(dg_l23_mat)=dg_l23_df$sample
dg_l23_mat=t(dg_l23_mat)
dg_l23_mat=dg_l23_mat[gene_to_plot,]
dg_l23=dg_l23_mat %>% na.omit()


anno=data.frame(sample=colnames(dg_l23),
                ROI=dg_l23_df$cluster)

anno=left_join(anno, metadata)
gene_anno=c(rep("Shared", sum(gene_ord_shared$gene %in% rownames(dg_l23))),
            rep("CTX L2/3", sum(gene_ord_cort$gene %in% rownames(dg_l23))),
            rep("DG-sg", sum(gene_ord_dg$gene %in% rownames(dg_l23))))

column_ha = HeatmapAnnotation(ROI=anno$ROI, Age =anno$age, Time =anno$time,
                              col=list(ROI=c("#984EA3","#4DAF4A") %>%
                                         setNames(c("Cortex Layer 2/3",
                                                    "Dentate Gyrus-sg")),
                                       Time=RColorBrewer::brewer.pal(4, "Set2") %>%
                                         setNames(paste0("ZT", 6*0:3)),
                                       Age=c(RColorBrewer::brewer.pal(12, "Paired")[11],
                                             RColorBrewer::brewer.pal(3, "BrBG")[1]) %>%
                                         setNames(paste0(c(7, 14), " months"))),
                              annotation_legend_param = list(ROI = list(
                                at = unique(anno$ROI)
                              )))

res=readRDS("par_im/res.rds") %>% filter(cluster %in% c("Dentate Gyrus-sg", "Cortex Layer 2/3"),
                                         gene %in% rownames(dg_l23))

relamp=left_join(coefs, res) %>%
  mutate(rel_amp=amp/log2(baseMean+1)) %>% na.omit()
relamp= relamp %>% group_by(gene) %>%
  summarise(rel_amp_fc=log2(rel_amp[2]/rel_amp[1]))
rownames(relamp)=relamp$gene
relamp=relamp[rownames(dg_l23),]

mean_expr=left_join(coefs, res) %>% 
  dplyr::filter(gene %in% rownames(dg_l23),
                cluster %in% c("Dentate Gyrus-sg","Cortex Layer 2/3")) %>%
  group_by(cluster, gene) %>%
  summarise(mean_expr=log2(mean(baseMean)+1)) %>% ungroup() %>%
  pivot_wider(names_from = "cluster", 
              values_from = "mean_expr",
              id_cols = "gene") %>% as.data.frame()

rownames(mean_expr)=mean_expr$gene
mean_expr$gene=NULL
mean_expr=mean_expr[rownames(dg_l23),]

coefs2=readRDS("par_im/coefs.rds") %>% 
  bind_rows(.id="cluster") %>% dplyr::select(gene, cluster, phi_hr)

phase_mat=left_join(coefs2, res) %>% 
  dplyr::filter(gene %in% rownames(dg_l23),
                cluster %in% c("Dentate Gyrus-sg","Cortex Layer 2/3")) %>%
  dplyr::select(cluster, gene, phi_hr) %>% distinct() %>%
  pivot_wider(names_from = "cluster", 
              values_from = "phi_hr",
              id_cols = "gene") %>% as.data.frame()

rownames(phase_mat)=phase_mat$gene
phase_mat$gene=NULL
phase_mat=phase_mat[rownames(dg_l23),]

row_ha1= rowAnnotation(Gene=gene_anno, 
                       col=list(Gene=RColorBrewer::brewer.pal(3, "Dark2")[1:3] %>%
                                  as.character %>%
                                  setNames(unique(gene_anno))),
                       annotation_legend_param = list(Gene = list(
                         at = unique(gene_anno)
                       )))

row_ha2 = rowAnnotation(
  Rel_Amp_L2FC =
    anno_points(relamp$rel_amp_fc,
                gp = gpar(
                  col = ifelse(relamp$rel_amp_fc > 0,
                               "orange", "gray60")
                )),
  
  Mean_Expression = anno_points(as.matrix(mean_expr),
                                gp = gpar(col =
                                            c("#984EA3",
                                              "#4DAF4A"))),
  Phase =
    anno_points(as.matrix(phase_mat),
                gp = gpar(col =
                            c("#984EA3",
                              "#4DAF4A" ))),
  width =unit(6, "cm")
  
)

Heatmap(dg_l23,  show_row_names=F,
        name="ht1",
        show_column_names = F,
        cluster_rows = F, 
        cluster_columns = F, 
        row_order = 1:nrow(dg_l23),
        column_order = 1:ncol(dg_l23),
        top_annotation = column_ha,
        left_annotation = row_ha1,
        right_annotation = row_ha2,
        heatmap_legend_param = list(
          title = "Z-scored\nLog2 Expr"))

decorate_heatmap_body("ht1", 
                      {grid.lines(c((length(unique(l23$sample)))/ncol(dg_l23), (length(unique(l23$sample)))/ncol(dg_l23)),
                                  c(0, 1),
                                  gp = gpar(lty = 2, lwd = 4))
                      }, slice=1)

decorate_heatmap_body("ht1", 
                      {grid.lines(c(0, 1),
                                  c(sum(gene_ord_dg$gene %in% rownames(dg_l23))/nrow(dg_l23), sum(gene_ord_dg$gene %in% rownames(dg_l23))/nrow(dg_l23)),
                                  gp = gpar(lty = 2, lwd = 4))
                      }, slice=1)

decorate_heatmap_body("ht1", 
                      {grid.lines(c(0, 1),
                                  c((sum(gene_ord_cort$gene %in% rownames(dg_l23))+sum(gene_ord_dg$gene %in% rownames(dg_l23)))/nrow(dg_l23), 
                                    (sum(gene_ord_cort$gene %in% rownames(dg_l23))+sum(gene_ord_dg$gene %in% rownames(dg_l23)))/nrow(dg_l23)),
                                  gp = gpar(lty = 2, lwd = 4))
                      }, slice=1)



time_lines=as.data.frame(metadata)
rownames(time_lines)=metadata$sample
time_lines=time_lines[colnames(dg_l23),]
prev="ZT0"
for(i in 1:nrow(time_lines)){
  if(time_lines$time[i]!=prev) {
    decorate_heatmap_body("ht1", 
                          {grid.lines(c((i-1)/ncol(dg_l23), 
                                        (i-1)/ncol(dg_l23)),
                                      c(0, 1),
                                      gp = gpar(lty = 1, lwd = 0.8, col="gray45"))
                          }, slice=1)
    
    
  }
  
  prev=time_lines$time[i]
  
}

