library(tidyverse)
library(glue)
library(ComplexHeatmap)
library(magrittr)
library(scico)
library(DESeq2)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_geno_age_int/")

cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")

deseqs=readRDS("age_geno_wald_full/deseqs.rds")

meta.data=readRDS("age_geno_wald_full/meta.data.rds")

meta.data=meta.data %>% mutate(grp=paste0(age, "_", genotype))

  res=lapply(deseqs, function(x) {
    results(x, 
            contrast =list( c("age_14.months_vs_7.months","genotypeAPP23.age14.months")),
            #name="genotypeAPP23.sexF",
            tidy = T)
  }) %>% 
    bind_rows(.id="cluster") %>%
    dplyr::rename(gene=row) %>%
    filter(cluster %in% names(cluster_color)[7:11])
  res_sig= res %>% dplyr::filter(padj<0.01)
    
  

sig_genes=res_sig

cnts=lapply(deseqs, function(df) {
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
  separate(cln, c("cluster",  "age", "genotype"), sep="_") %>%
  mutate(cluster=factor(cluster, levels=names(cluster_color))) %>%
  mutate(age=relevel(factor(age), "7 months"),
         genotype=relevel(factor(genotype), "WT")) %>%
  arrange(cluster, genotype,age )


to_plt=cnts_sum %>% ungroup %>%
  dplyr::select(-gene) %>% as.matrix()

rownames(to_plt)=cnts_sum$gene
gene_ord=sig_genes %>% group_by(gene) %>% summarise(lfc=mean(log2FoldChange, na.rm=T)) %>%
  arrange(-lfc)

to_plt=to_plt[,colmn_ord$ord-1]
to_plt=to_plt[gene_ord$gene,25:44]
to_plt=na.omit(to_plt)

colmn_ord=colmn_ord[25:44,]
colmn_ord$geno=plyr::mapvalues(colmn_ord$genotype, c("WT", "APP23"), c("NTG", "APP23-TG"))
column_ha1 = HeatmapAnnotation(Age =colmn_ord$age,
                               Genotype=colmn_ord$geno,
                               
                               col=list(Age=c(RColorBrewer::brewer.pal(12, "Paired")[11],
                                              RColorBrewer::brewer.pal(3, "BrBG")[1]) %>%
                                          setNames(paste0(c(7, 14), " months")),
                                        
                                        
                                        Genotype=c("#0072B5FF","#BC3C29FF") %>%
                                          setNames(c("NTG",
                                                     "APP23-TG"))),
                               annotation_legend_param = list(Age = list(
                                 at = paste0(c(7, 14), " months")),
                                 
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
                                 labels = make.names(as.character(clst_mark_loc$cluster)),
                                 which = "column",
                                 labels_rot = 70,
                                 link_height = unit(.5, "mm")),
                             col=list(Cluster=cluster_color),
                             show_legend = F,
                             show_annotation_name = F)




g_to_mark=c("Bc1","C4b" , "Gfap","Crhr1" ,"S100b")

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
                                      gp = gpar(lty = 1, lwd = 1, col="gray25"))
                          }, slice=1)
    
    
  }
  
  prev=colmn_ord$cluster[i]
  
}

prev=colmn_ord$genotype[1]
for(i in 1:nrow(colmn_ord)){
  if(colmn_ord$genotype[i]!=prev) {
    decorate_heatmap_body("ht1", 
                          {grid.lines(c((i-1)/ncol(to_plt), 
                                        (i-1)/ncol(to_plt)),
                                      c(0, 1),
                                      gp = gpar(lty = 2, lwd = 0.9, col="gray45"))
                          }, slice=1)
    
    
  }
  
  prev=colmn_ord$genotype[i]
  
}


