library(tidyverse)
library(glue)
library(magrittr)
library(scico)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_app/")

cluster_order = readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/objects/cluster_order.rds")
cluster_color = readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")

grps=do.call(paste0, expand.grid(c("old", "young"), c("_M", "_F")))


res_comb={
  res=lapply(grps %>% setNames(grps), function(gr) {
    readRDS(glue(gr,"/res.rds")) %>% 
      mutate(cluster=gsub("Olfactory Tubercle", "Amygdala", cluster))
  })
  res_comb=bind_rows(res, .id="grp")
  res_comb=res_comb %>% separate(grp, into = c("age", "sex"), sep = "_", remove = F)
  res_comb
}

all_sig=readRDS("res/res.rds") %>% filter(padj<0.1)

all_sig$cluster=gsub("Olfactory Tubercle", "Amygdala", all_sig$cluster)

all_sig=split(all_sig, all_sig$cluster)

res_comb_splt=split(res_comb, res_comb$cluster)

gene_list=lapply(names(all_sig) %>% setNames(names(all_sig)), 
                 function(clst) {
                   
  x=all_sig[[clst]]
  y=res_comb_splt[[clst]] %>% filter(gene %in% x$gene) %>%
  mutate(dir=ifelse(log2FoldChange>0, "up", "down")) %>%
    group_by(gene) %>% mutate(concrd=max(table(dir)/n())) %>%
    filter(concrd==1, mean(abs(log2FoldChange))>.3,
           max(abs(log2FoldChange))>.7) %>% 
    mutate(mlfc=mean(abs(log2FoldChange)), 
           mpval=exp(mean(log(padj))))
  
  if(nrow(y)==0){return(NULL)} else{
    y=y %>%
      filter(age=="old", sex=="F") 
 return(data.frame(cluster=clst, 
                   gene=y$gene, 
                   dir=y$dir,
                   mlfc=y$mlfc,
                   mpval=y$mpval))
  }
})

sig_genes=bind_rows(gene_list) %>%
  mutate(cluster=factor(cluster, levels = cluster_order)) %>% 
  arrange(dir,cluster) %>% group_by(gene) %>% 
  filter(mlfc==max(mlfc)) %>%
  filter(!duplicated(gene))


to_plt=res_comb %>% dplyr::filter(gene %in% sig_genes$gene) %>% 
  mutate(grp_clust=paste0(cluster,"--", grp))%>% 
  group_by(cluster, gene) %>%
  ungroup() %>%
  dplyr::select(grp_clust, gene, log2FoldChange ) %>% 
  pivot_wider(names_from = "grp_clust", 
              values_from = log2FoldChange, 
              id_cols = "gene", values_fill = 0) %>%
  column_to_rownames("gene") %>% as.matrix()

colmn_ord=res_comb %>% dplyr::select(grp, cluster) %>%
  distinct() %>%
  mutate(cluster=factor(cluster, levels = cluster_order),
         grp=factor(grp, levels=rev(sort(grps)))) %>%
  arrange(cluster,grp) %>%
  mutate(grp_clust=paste0(cluster,"--", grp)) %>%
  separate(grp, into = c("age", "sex"), sep = "_", remove = F) %>%
  mutate(sex=plyr::mapvalues(sex,c("M", "F"), c("Male", "Female")),
         age=plyr::mapvalues(age,c("old", "young"), paste0(c(14, 7), " months")))
  
  
res_comb_filt=res_comb %>% dplyr::filter(gene %in% sig_genes$gene)

to_plt=to_plt[, colmn_ord$grp_clust]


to_plt=to_plt[sig_genes$gene,]
colmn_ord=colmn_ord[1:44,]
to_plt=to_plt[,1:44]

library(ComplexHeatmap)

column_ha1 = HeatmapAnnotation(Age =colmn_ord$age,
                              Sex =colmn_ord$sex,
                              col=list(Age=c(RColorBrewer::brewer.pal(12, "Paired")[11],
                                             RColorBrewer::brewer.pal(3, "BrBG")[1]) %>%
                                         setNames(paste0(c(7, 14), " months")),
                                       
                                       Sex=as.character(wesanderson::wes_palette("GrandBudapest2",2)) %>%
                                         setNames(c("Female", "Male"))),
                              annotation_legend_param = list(Age = list(
                                at = paste0(c(7, 14), " months")),
                                Sex=list(
                                  at=c("Male", "Female")
                                )
                              ),
                              show_legend = c(F,F)
                              )

column_ha2 = HeatmapAnnotation(Cluster=colmn_ord$cluster,
                               col=list(Cluster=cluster_color),
                               show_legend = F,
                               show_annotation_name = F
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

g_to_mark=c("B2m","C1qa","Cst7","Crh","Fabp7","Gfap","Syp" , "Pvalb", "Shh","Npy","Msmo1", "Crhbp")

g_ind=sapply(g_to_mark, function(x) {grep(paste0(x,"$"), rownames(to_plt))})


row_ha2= rowAnnotation(genes = anno_mark(at = as.numeric(g_ind),
                                         side="right",
                                         labels = g_to_mark))

mm=quantile(abs(to_plt),0.98,na.rm=T)

ht1=Heatmap(-1*to_plt, 
        name="ht1",
        cluster_columns = F, 
        col=circlize::colorRamp2(seq(-mm, mm, length = 200), 
                             scico(200,palette='vik', direction = 1)),
        cluster_rows = F, 
        show_row_names = F,
        show_column_names = F,
        top_annotation = column_ha1,
        bottom_annotation = clst_names,
        right_annotation = row_ha2,
        heatmap_legend_param = list(
          title = "Log2FC\nAPP23-TG vs. NTG",
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

prev=colmn_ord$age[1]
for(i in 1:nrow(colmn_ord)){
  if(colmn_ord$age[i]!=prev) {
    decorate_heatmap_body("ht1", 
                          {grid.lines(c((i-1)/ncol(to_plt), 
                                        (i-1)/ncol(to_plt)),
                                      c(0, 1),
                                      gp = gpar(lty = 2, lwd = 0.9, col="gray45"))
                          }, slice=1)
    
    
  }
  
  prev=colmn_ord$age[i]
  
}


