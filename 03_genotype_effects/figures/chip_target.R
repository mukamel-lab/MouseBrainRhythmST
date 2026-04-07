library(tidyverse)



# setwd('~/desp1/precast/precast_final_with_ros_caud/analysis/chipatlas/target_genes')
# 
# chip_tfs=pbmclapply(list.files(),function(y) {vroom::vroom(y, col_select = 1, n_max=100 ) %>% mutate(TF=gsub("\\.5\\.tsv$", "", y))})
# 
# names(chip_tfs)=gsub("\\.5\\.tsv$", "",list.files())
chip_tfs=readRDS("~/desp1/precast/precast_final_with_ros_caud/analysis/chipatlas/chiptargets.rds")
setwd("~/desp1/precast/prec_c25q25g3000")
region=readRDS("objects/region_order.rds")
cluster_order=readRDS("objects/cluster_order.rds")

# coefs_o = readRDS("deseq_rhth_int2/geno_int_new/old/interaction_results_readjusted_o.rds")%>% dplyr::filter(fdr<0.2)
# coefs_y=readRDS("deseq_rhth_int2/geno_int_new/young/interaction_results_readjusted_y.rds") %>% mutate(dif=abs(wt_amp-app_amp)) %>%
#   group_by(cluster) %>% slice_min(dif, prop = 0.1)
# coefs_y=readRDS("deseq_rhth_int2/geno_int_new/young/interaction_results_readjusted_y.rds") %>%
#   group_by(cluster) %>% slice_min(padj, prop = 0.1)

coefs_y=readRDS("deseq_rhth_int2/geno_int_new/young/interaction_results_readjusted_y.rds") %>% dplyr::filter(fdr<0.1) %>%
  mutate(region=plyr::mapvalues(cluster, cluster_order, region))

# coefs_y=readRDS("deseq_rhth_int2/geno_int_new/old/interaction_results_readjusted_o.rds") %>% dplyr::filter(fdr<0.2) %>%
#   mutate(region=plyr::mapvalues(cluster, cluster_order, region))

tf_target_mrv_y=pbmclapply(chip_tfs, function(tf) {
  cf=coefs_y %>% dplyr::filter(gene %in% tf$Target_genes) %>% group_by(region) %>%
    dplyr::filter(length(unique(gene))>5)%>%
    mutate(wtts = t_s + genotypeWT.t_s,
           wttc = t_c + genotypeWT.t_c) %>%
    summarise(
      tc = mean(t_c),
      ts = mean(t_s),
      wtts = mean(wtts),
      wttc = mean(wttc),
      num_genes=length(unique(gene)),
      num_deg=n()
    ) %>%
    mutate(
      amp_wt = sqrt(wtts ^ 2 + wttc ^ 2),
      phihr_wt = (12 / pi) * atan2(wtts, wttc) %% (2 * pi),
      amp_app = sqrt(ts ^ 2 + tc ^ 2),
      phihr_app = (12 / pi) * atan2(ts, tc) %% (2 * pi)
    )
}) %>% bind_rows(.id="tf")


tf_target_mrv_y=tf_target_mrv_y %>% group_by(region) %>% slice_max(order_by = num_genes,n=5 )%>% slice_max( order_by=num_deg, n=5)

tf_target_mrv_y=tf_target_mrv_y%>%
  pivot_longer(
    cols      = c(starts_with("amp_"), starts_with("phihr_")),
    names_to  = c(".value", "genotype"),
    names_sep = "_"
  )

pctx=ggplot(data=tf_target_mrv_y %>% dplyr::filter(region=="Cortex"))+
  geom_segment( aes(x = phihr, xend = phihr, y = 0, yend = amp, col = genotype),
                arrow = arrow(length = unit(0.05, "npc")), linewidth = .2) +
  coord_polar(start = 0) +facet_grid2(region~tf, render_empty = F, switch = "y")+
  scale_x_continuous(breaks = 0:3 * 6, limits = c(0, 24)) +
  scale_color_manual(name="",
                     values = ggsci::pal_nejm()(2) %>% setNames(c("app","wt")),
                     breaks=c("app","wt"),labels=c("APP23","NTG")
  )+
  theme_minimal(base_size = 6, base_family = "ArialMT") +
  theme(strip.placement = "outside",
    
    panel.border = element_blank(),       
    panel.spacing = unit(0.2, "lines"),
    axis.title = element_blank(),
    axis.text.x = element_blank())+scale_y_continuous(breaks=pretty_breaks(n=3))

phip=ggplot(data=tf_target_mrv_y %>% dplyr::filter(region=="Hippocampus"))+
  geom_segment( aes(x = phihr, xend = phihr, y = 0, yend = amp, col = genotype),
                arrow = arrow(length = unit(0.05, "npc")), linewidth = .2) +
  coord_polar(start = 0) +facet_grid2(region~tf, render_empty = F, switch = "y")+
  scale_x_continuous(breaks = 0:3 * 6, limits = c(0, 24)) +
  scale_color_manual(name="",
                     values = ggsci::pal_nejm()(2) %>% setNames(c("app","wt")),
                     breaks=c("app","wt"),labels=c("APP23","NTG")
  )+
  theme_minimal(base_size = 6, base_family = "ArialMT") +
  theme(strip.placement = "outside",
        
        panel.border = element_blank(),       
        panel.spacing = unit(0.2, "lines"),
        axis.title = element_blank(),
        axis.text.x = element_blank())+
  scale_y_continuous(breaks=pretty_breaks(n=3))
phip

tfeb_targets=readRDS("deseq_rhth_int2/geno_int_new/young/interaction_results_readjusted_y.rds") %>% 
  mutate(region=plyr::mapvalues(cluster, cluster_order, region))%>% 
  dplyr::filter(region=="Cortex", gene %in%  chip_tfs$Tfeb$Target_genes)


tfeb_long=tfeb_targets %>%
  dplyr::select(-ends_with("phi"))%>%
  pivot_longer(
    cols      = c(ends_with("_amp"), ends_with("phi_hr")),
    names_to  = c( "genotype",".value"),
    names_sep = "_"
  )
res_y=readRDS("deseq_rhth_int2/geno_int_new/young/res_y_aw.rds") %>% dplyr::filter(grepl("^Cortex|CTX", cluster), gene %in% tfeb_long$gene) %>%
  mutate(genotype=plyr::mapvalues(genotype, c("APP23","WT"),c("app","wt")))

tfeb_long=left_join(tfeb_long, res_y %>% dplyr::select(gene,cluster, padj, genotype) %>% 
                      dplyr::rename(rhth=padj), by = c("gene","cluster", "genotype"))
gtl=tfeb_long %>% dplyr::filter(rhth < .1) %>% slice_max(order_by = amp, n=10)

pctx = ggplot() +
  geom_point(data = tfeb_long %>% dplyr::filter(rhth >= .1),
             aes(x = phi, y = amp),
             col = "gray",size=0.15) +
  geom_point(data = tfeb_long %>% dplyr::filter(rhth < .1),
             aes(x = phi, y = amp, col = genotype),size=0.2) +
  geom_text_repel(
    data =tfeb_long %>% dplyr::filter(rhth < .1,fdr<0.1, gene!="Bhlhe40"),
    aes(x = phi, y = amp, label = gene),
    size = 6 * 5 / 14,
    family = "ArialMT",
    hjust = -0.2,
    vjust = 0.3,
    color = "black",
    force = 2,
    nudge_x = 0.03,
    nudge_y = 0.03,
    max.overlaps = 10,
    segment.size = 0.07,
    min.segment.length = .1
  ) +
  
  coord_polar(start = 0) +
  scale_x_continuous(breaks = 0:3 * 6, limits = c(0, 24)) +
  scale_color_manual(
    name = "",
    values = ggsci::pal_nejm()(2) %>% setNames(c("app", "wt")),
    breaks = c("app", "wt"),
    labels = c("APP23", "NTG")
  ) +
  theme_minimal(base_size = 6, base_family = "ArialMT") +
  theme(
    strip.placement = "outside",
    
    panel.border = element_blank(),
    panel.spacing = unit(0.2, "lines"),
    axis.title = element_blank(),
    axis.text.x = element_blank()
  ) + scale_y_continuous(breaks = pretty_breaks(n = 3))
pctx
