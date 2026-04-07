library(tidyverse)
library(glue)
library(ComplexHeatmap)
library(magrittr)
library(scico)
library(DESeq2)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_geno_sex_int/")

cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")

deseqs_old=readRDS("sex_geno_old_wald_full/deseqs.rds")
deseqs_young=readRDS("sex_geno_young_wald_full/deseqs.rds")
all_dsqs=readRDS("../deseq_app/all_wald/deseqs.rds")
sex_genes=c("Xist", "Tsix","Eif2s3y", "Ddx3y",  "Uty","Kdm5d", "Ddx3x","Eif2s3x",
            "Eif2sx")
meta.data=readRDS("../deseq_app/all_wald/meta.data.rds")

meta.data=meta.data %>% mutate(grp=paste0(sex, "_", age, "_", genotype))

sig_genes=lapply(list(old=deseqs_old, young=deseqs_young), function(deseqs) {
  res=lapply(deseqs, function(x) {
    results(x, 
            #contrast =list( c("genotype_APP23_vs_WT","genotypeAPP23.sexF")),
            #name="genotypeAPP23.sexF",
            contrast =list( c("sex_F_vs_M","genotypeAPP23.sexF")),
            tidy = T)
  }) %>% 
    bind_rows(.id="cluster") %>% dplyr::filter(!(row %in% sex_genes)) %>%
    dplyr::rename(gene=row)
  res_sig= res %>% dplyr::filter(padj<0.01, grepl("Cortex L|CTX", cluster))%>%
    group_by(cluster) %>%
    top_n(7, abs(log2FoldChange))
  res_sig
}) %>% 
  bind_rows(.id="age")


cnts=lapply(all_dsqs, function(df) {
  counts(df, normalized=T) %>% 
    as.data.frame() %>% 
    rownames_to_column("gene") %>%
    filter(gene %in% sig_genes$gene)
}) %>% bind_rows(.id="cluster")

cnts_sum=cnts %>% 
  pivot_longer(cols = meta.data$sample, 
               names_to = "sample",
               values_to = "cpm") %>%
  mutate(grp=plyr::mapvalues(sample, 
                             meta.data$sample, 
                             meta.data$grp)) %>%
  group_by(grp, gene, cluster) %>%
  summarise(cpm=mean(cpm, na.rm=T)) %>%
  mutate(cpm=log2(cpm+1))%>%
  ungroup() %>%
  group_by(gene, cluster) %>%
  mutate(cpm=scale(cpm)) %>%
  ungroup() %>%
  pivot_wider(names_from = c(cluster,grp),
              values_from = cpm)

colmn_ord=data.frame(cln=tail(colnames(cnts_sum),-1),ord=2:ncol(cnts_sum)) %>% 
  separate(cln, c("cluster", "sex", "age", "genotype"), sep="_") %>%
  mutate(cluster=factor(cluster, levels=names(cluster_color))) %>%
  mutate(age=relevel(factor(age), "7 months"),
         sex=relevel(factor(sex), "M"),
         genotype=relevel(factor(genotype), "WT")) %>%
  arrange(cluster, age,sex ,genotype ) %>% filter(age=="14 months")


to_plt=cnts_sum %>% ungroup %>%
  dplyr::select(-gene) %>% as.matrix()

rownames(to_plt)=cnts_sum$gene
gene_ord=sig_genes %>% group_by(gene) %>% summarise(lfc=mean(log2FoldChange, na.rm=T)) %>%
  arrange(-lfc)

to_plt=to_plt[,colmn_ord$ord-1]
to_plt=to_plt[gene_ord$gene,1:(8*3)]
to_plt=na.omit(to_plt)

colmn_ord=colmn_ord[1:24,]
colmn_ord$geno=plyr::mapvalues(colmn_ord$genotype, c("WT", "APP23"), c("NTG", "APP23-TG"))
column_ha1 = HeatmapAnnotation(Sex =colmn_ord$sex,
                               Genotype=colmn_ord$geno,
                               col=list(Sex=as.character(wesanderson::wes_palette("GrandBudapest2",2)) %>%
                                          setNames(c("F", "M")),
                                        Genotype=c("#0072B5FF","#BC3C29FF") %>%
                                          setNames(c("NTG",
                                                     "APP23-TG"))),
                               annotation_legend_param = list(
                                 Sex=list(
                                   at=c("M", "F")
                                 ),
                                 Genotype=list(at=c("NTG",
                                                    "APP23-TG"))
                               ),
                               show_legend = c(T,T)
)
clst_mark_loc=colmn_ord[duplicated(colmn_ord$cluster),]
clst_mark_loc=clst_mark_loc[!duplicated(clst_mark_loc$cluster),]

clst_names=HeatmapAnnotation(Cluster=colmn_ord$cluster,
                             clst_lab=
                               anno_mark(
                                 at = as.numeric(rownames(clst_mark_loc)),
                                 side="bottom",
                                 labels = clst_mark_loc$cluster,
                                 which = "column",
                                 labels_rot = 70,
                                 link_height = unit(.5, "mm")),
                             col=list(Cluster=cluster_color),
                             show_legend = F,
                             show_annotation_name = F)




#g_to_mark=c("B2m","Cst7","Crhbp","Msmp1","Gfap","Oxt", "Avp","Igfbp2", "Cort")
#g_to_mark=c("Agap2", "Vim")
g_to_mark=c("Serpina3n", "Vim", "Scg5","Ctss","Cartpt", "Nfix", "Resp18")
g_ind=sapply(g_to_mark, function(x) {grep(paste0(x,"$"), rownames(to_plt))})


row_ha2= rowAnnotation(genes = anno_mark(at = as.numeric(g_ind),
                                         side="right",
                                         labels = g_to_mark))

mm=quantile(abs(to_plt),0.98,na.rm=T)

ht1=Heatmap(
  name="ht1",
  to_plt, 
  cluster_columns = F,
  cluster_rows = F, 
  show_column_names = F,
  show_row_names = F,
  top_annotation = column_ha1 ,
  bottom_annotation = clst_names,
  col=circlize::colorRamp2(seq(-1*mm, mm, length = 50), 
                                 scico(50,palette='bam')),
  right_annotation = row_ha2,
  heatmap_legend_param = list(
          title = "Z-scored Expression",
          title_position = "leftcenter-rot",
          ncol=1
          
        ))





draw(ht1,
     heatmap_legend_side = "left")


prev=colmn_ord$cluster[1]
for(i in 1:nrow(colmn_ord)){
  if(colmn_ord$cluster[i]!=prev) {
    decorate_heatmap_body("ht1", 
                          {grid.lines(c((i-1)/ncol(to_plt), 
                                        (i-1)/ncol(to_plt)),
                                      c(0, 1),
                                      gp = gpar(lty = 1, lwd = 1.5, col="black"))
                          }, slice=1)
    
    
  }
  
  prev=colmn_ord$cluster[i]
  
}

prev=colmn_ord$sex[1]
for(i in 1:nrow(colmn_ord)){
  if(colmn_ord$sex[i]!=prev) {
    decorate_heatmap_body("ht1", 
                          {grid.lines(c((i-1)/ncol(to_plt), 
                                        (i-1)/ncol(to_plt)),
                                      c(0, 1),
                                      gp = gpar(lty = 2, lwd = 1.5, col="black"))
                          }, slice=1)
    
    
  }
  
  prev=colmn_ord$sex[i]
  
}


