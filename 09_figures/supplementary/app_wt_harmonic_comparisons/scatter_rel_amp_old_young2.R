library(tidyverse)
library(ggrepel)
library(patchwork)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth/old_young_sep/")
clock_g=c(readRDS("~/desp/clock_genes.rds"), "Per3")
cluster_color=readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
cluster_order = readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/objects/cluster_order.rds")
cluster_color=colorspace::darken(cluster_color, .2) %>% setNames(names(cluster_color))
rhyth_genes=rbind(readRDS("../APP_harmonic/par_im/res.rds"),
                  readRDS("../WT_harmonic/par_im/res.rds")) %>% 
  filter(gene %in% clock_g|(baseMean>60 |  padj<0.05))
rhyth_genes=split(rhyth_genes, rhyth_genes$cluster)
rhyth_genes= lapply(rhyth_genes, function(df) {df %>% .[["gene"]]})

young_wt=left_join(readRDS("WT_young_harmonic/par_im/res.rds"),
                   readRDS("WT_young_harmonic/par_im/coefs.rds") %>%
                     bind_rows(.id="cluster")) %>%
  mutate(Age="7 months",
         Genotype="NTG",
         l2expr=log2(baseMean+1),
         rel_amp=amp/l2expr)

young_wt=split(young_wt, young_wt$cluster)

young_wt=lapply(names(young_wt), function(clst) {
  
  young_wt[[clst]] %>% filter(gene %in% rhyth_genes[[clst]])
  
}) %>% setNames(names(young_wt)) %>%
  bind_rows(.id="cluster")
young_app=left_join(readRDS("APP_young_harmonic/par_im/res.rds"),
                    readRDS("APP_young_harmonic/par_im/coefs.rds") %>%
                      bind_rows(.id="cluster")) %>%
  mutate(Age="7 months",
         Genotype="APP23-TG",
         l2expr=log2(baseMean+1),
         rel_amp=amp/l2expr)

young_app=split(young_app, young_app$cluster)

young_app=lapply(names(young_app), function(clst) {
  
  young_app[[clst]] %>% filter(gene %in% rhyth_genes[[clst]])
  
}) %>% setNames(names(young_app)) %>%
  bind_rows(.id="cluster")

young=bind_rows(list(young_app, young_wt)) 

young_sum=young %>% group_by(gene, cluster) %>% arrange(Genotype) %>%
  summarise(comp=paste0(Genotype[1],"/", Genotype[2]),
            relamp_fc=log2(rel_amp[1]/rel_amp[2])) %>% na.omit()

lfc_tally=young_sum %>% ungroup() %>% group_by(cluster) %>% summarise(up_prop=sum(relamp_fc>0)/n())
plt_young=split(young, young$cluster)




# rhyth_genes=readRDS("../WT_harmonic/par_im/res_sig.rds")
# rhyth_genes=split(rhyth_genes, rhyth_genes$cluster)
# rhyth_genes= lapply(rhyth_genes, function(df) {df %>% .[["gene"]]})

old_wt=left_join(readRDS("WT_old_harmonic/par_im/res.rds"),
                 readRDS("WT_old_harmonic/par_im/coefs.rds") %>%
                   bind_rows(.id="cluster")) %>%
  mutate(Age="7 months",
         Genotype="NTG",
         l2expr=log2(baseMean+1),
         rel_amp=amp/l2expr)

old_wt=split(old_wt, old_wt$cluster)

old_wt=lapply(names(old_wt), function(clst) {
  
  old_wt[[clst]] %>% filter(gene %in% rhyth_genes[[clst]])
  
}) %>% setNames(names(old_wt)) %>%
  bind_rows(.id="cluster")
old_app=left_join(readRDS("APP_old_harmonic/par_im/res.rds"),
                  readRDS("APP_old_harmonic/par_im/coefs.rds") %>%
                    bind_rows(.id="cluster")) %>%
  mutate(Age="7 months",
         Genotype="APP23-TG",
         l2expr=log2(baseMean+1),
         rel_amp=amp/l2expr)

old_app=split(old_app, old_app$cluster)

old_app=lapply(names(old_app), function(clst) {
  
  old_app[[clst]] %>% filter(gene %in% rhyth_genes[[clst]])
  
}) %>% setNames(names(old_app)) %>%
  bind_rows(.id="cluster")

old=bind_rows(list(old_app, old_wt))

old_sum=old %>% group_by(gene, cluster) %>% arrange(Genotype) %>%
  summarise(comp=paste0(Genotype[1],"/", Genotype[2]),
            relamp_fc=log2(rel_amp[1]/rel_amp[2])) %>% na.omit()

lfc_tally=old_sum %>% ungroup() %>% group_by(cluster) %>% summarise(up_prop=sum(relamp_fc>0)/n())
plt_old=split(old, old$cluster)

##############
#plot
###############

plts_comb = lapply(cluster_order[c(1:6)], function(clst) {
  genes_to_lab = c("Per2",  "Per3", "Arc")
  
  rel_amp_thresh=0.03
  
  df_old = plt_old[[clst]] %>% select(gene, Genotype, rel_amp) %>%
    pivot_wider(names_from = "Genotype", values_from = "rel_amp") %>%
    filter((`APP23-TG` > rel_amp_thresh & NTG > rel_amp_thresh) | gene %in% genes_to_lab)
  df_young = plt_young[[clst]] %>% select(gene, Genotype, rel_amp) %>%
    pivot_wider(names_from = "Genotype", values_from = "rel_amp") %>%
    filter((`APP23-TG` > rel_amp_thresh & NTG > rel_amp_thresh) | gene %in% genes_to_lab)
  # genes_to_lab=unique(c(df_young %>% filter(`APP23-TG`>0.1 | NTG>0.1) %>% .[["gene"]],
  #                       clock_g))
  #genes_to_lab=clock_g
  
  
  df_plt = rbind(df_old %>% mutate(Age = "14 months"),
                 df_young %>% mutate(Age = "7 months")) %>%
    mutate(Age = relevel(factor(Age), "7 months"))
  #genes_to_lab = c("Per2",  "Per3", "Arc", "Opalin")
  
  lmts = max(quantile(c(df_plt$`APP23-TG`, df_plt$NTG), .99),

             max(df_plt %>% filter(gene %in% genes_to_lab) %>%
                   select(-gene,-Age)))+0.02

  genes_to_lab_extra=c(df_plt %>% filter(`APP23-TG`<=lmts) %>%
                         top_n(2,`APP23-TG`) %>% .[["gene"]],
                       df_plt %>% filter(NTG<=lmts) %>%
                         top_n(2,NTG) %>% .[["gene"]])
  
  #genes_to_lab=grep("Rik",unique(c(gene_to_lab, genes_to_lab_extra)), value = T, invert = T)
  #print(genes_to_lab)
  p3 = suppressMessages({df_plt %>%
    ggplot(aes(y = `APP23-TG`, x = NTG, group=gene)) + geom_point(color = cluster_color[6], size =
                                                        .5) +
    geom_label_repel(
      data = . %>% filter(gene %in% genes_to_lab),
      aes(y = `APP23-TG`, x = NTG, label = gene),
      color = cluster_color[6],
      min.segment.length = unit(0, 'mm'),
      nudge_y = .04,
      nudge_x = .04,
      direction = "both"
    ) +
    #geom_density2d()+
    geom_abline(slope = 1, intercept = 0) + ggtitle(clst) +
    coord_fixed() + theme_bw(base_size = 14) +
    xlim(c(NA, lmts)) + ylim(c(NA, lmts)) +
    facet_wrap( ~ Age, nrow = 1) + theme(plot.title = element_text(hjust =
                                                                     .5))+
    scale_x_log10()+scale_y_log10()})

  
  p3
})
names(plts_comb) = cluster_order[c(1:6)]
plts_comb$`Cortex Layer 2/3`

wrap_plots(plts_comb[cluster_order[c(1:6)]])

names(plts_comb)=gsub(" \\(L6b, CLA, EP\\)", "",
gsub("Dentate Gyrus",  "DG",
gsub("Cortex Layer ", "CTX L",
     names(plts_comb))))

names(plts_comb)=make.names(names(plts_comb))
for(p1 in names(plts_comb)) {
  pdf(paste0(p1, ".pdf"), paper = "a4r")
  plot(plts_comb[[p1]])
  dev.off()
  }

dg={
  clst=cluster_order[11]
  
genes_to_lab = res_sig %>% filter(cluster==clst) %>% .[["gene"]]

rel_amp_thresh=0.01

df_old = plt_old[[clst]] %>% select(gene, Genotype, rel_amp) %>%
  pivot_wider(names_from = "Genotype", values_from = "rel_amp") %>%
  filter((`APP23-TG` > rel_amp_thresh & NTG > rel_amp_thresh) | gene %in% genes_to_lab)
df_young = plt_young[[clst]] %>% select(gene, Genotype, rel_amp) %>%
  pivot_wider(names_from = "Genotype", values_from = "rel_amp") %>%
  filter((`APP23-TG` > rel_amp_thresh & NTG > rel_amp_thresh) | gene %in% genes_to_lab)
# genes_to_lab=unique(c(df_young %>% filter(`APP23-TG`>0.1 | NTG>0.1) %>% .[["gene"]],
#                       clock_g))
#genes_to_lab=clock_g


df_plt = rbind(df_old %>% mutate(Age = "14 months"),
               df_young %>% mutate(Age = "7 months")) %>%
  mutate(Age = relevel(factor(Age), "7 months"))
#genes_to_lab = c("Per2",  "Per3", "Arc", "Opalin")

lmts = max(quantile(c(df_plt$`APP23-TG`, df_plt$NTG), .99),
           
           max(df_plt %>% filter(gene %in% genes_to_lab) %>%
                 select(-gene,-Age)))+0.02

genes_to_lab_extra=c(df_plt %>% filter(`APP23-TG`<=lmts) %>%
                       top_n(2,`APP23-TG`) %>% .[["gene"]],
                     df_plt %>% filter(NTG<=lmts) %>%
                       top_n(2,NTG) %>% .[["gene"]])

#genes_to_lab=grep("Rik",unique(c(gene_to_lab, genes_to_lab_extra)), value = T, invert = T)
#print(genes_to_lab)
p3 = suppressMessages({df_plt %>%
    ggplot(aes(y = `APP23-TG`, x = NTG, group=gene)) + geom_point(color = cluster_color[6], size =
                                                                    .5) +
    geom_label_repel(
      data = . %>% filter(gene %in% genes_to_lab),
      aes(y = `APP23-TG`, x = NTG, label = gene),
      color = cluster_color[6],
      min.segment.length = unit(0, 'mm'),
      nudge_y = .04,
      nudge_x = .04,
      direction = "both"
    ) +
    #geom_density2d()+
    geom_abline(slope = 1, intercept = 0) + ggtitle(clst) +
    coord_fixed() + theme_bw(base_size = 14) +
    xlim(c(NA, lmts)) + ylim(c(NA, lmts)) +
    facet_wrap( ~ Age, nrow = 1) + theme(plot.title = element_text(hjust =
                                                                     .5))+
    scale_x_log10()+scale_y_log10()})


p3
}
dg




plts_comb = lapply(cluster_order, function(clst) {
  try({genes_to_lab = c("Per2",  "Per3", "Arc")
  
  rel_amp_thresh=0.03
  
  df_old = plt_old[[clst]] %>% select(gene, Genotype, rel_amp) %>%
    pivot_wider(names_from = "Genotype", values_from = "rel_amp") %>%
    filter((`APP23-TG` > rel_amp_thresh & NTG > rel_amp_thresh) | gene %in% genes_to_lab)
  df_young = plt_young[[clst]] %>% select(gene, Genotype, rel_amp) %>%
    pivot_wider(names_from = "Genotype", values_from = "rel_amp") %>%
    filter((`APP23-TG` > rel_amp_thresh & NTG > rel_amp_thresh) | gene %in% genes_to_lab)
  # genes_to_lab=unique(c(df_young %>% filter(`APP23-TG`>0.1 | NTG>0.1) %>% .[["gene"]],
  #                       clock_g))
  #genes_to_lab=clock_g
  
  
  df_plt = rbind(df_old %>% mutate(Age = "14 months"),
                 df_young %>% mutate(Age = "7 months")) %>%
    mutate(Age = relevel(factor(Age), "7 months"))
  #genes_to_lab = c("Per2",  "Per3", "Arc", "Opalin")
  
  lmts = max(quantile(c(df_plt$`APP23-TG`, df_plt$NTG), .99),
             
             max(df_plt %>% filter(gene %in% genes_to_lab) %>%
                   select(-gene,-Age)))+0.02
  
  genes_to_lab_extra=c(df_plt %>% filter(`APP23-TG`<=lmts) %>%
                         top_n(2,`APP23-TG`) %>% .[["gene"]],
                       df_plt %>% filter(NTG<=lmts) %>%
                         top_n(2,NTG) %>% .[["gene"]])
  
  #genes_to_lab=grep("Rik",unique(c(gene_to_lab, genes_to_lab_extra)), value = T, invert = T)
  #print(genes_to_lab)
  p3 = suppressMessages({df_plt %>%
      ggplot(aes(y = `APP23-TG`, x = NTG, group=gene)) + geom_point(color = cluster_color[6], size =
                                                                      .5) +
      
      #geom_density2d()+
      geom_abline(slope = 1, intercept = 0) + ggtitle(clst) +
      coord_fixed() + theme_bw(base_size = 14) +
      xlim(c(NA, lmts)) + ylim(c(NA, lmts)) +
      facet_wrap( ~ Age, nrow = 1) + theme(plot.title = element_text(hjust =
                                                                       .5))+
      scale_x_log10()+scale_y_log10()})
  
  
  p3})
})
names(plts_comb) = cluster_order

