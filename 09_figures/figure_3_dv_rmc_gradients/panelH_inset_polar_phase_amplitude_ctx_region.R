library(tidyverse)
library(scattermore)
library(ggh4x)
library(ggrastr)
library(ggrepel)

setwd("~/desp1/precast/precast_final_with_ros_caud")

short_names = readRDS("~/desp1/precast/precast_final_with_ros_caud/short_names.rds")

short_namesdf=data.frame(full=short_names,short=names(short_names))%>%
  separate(full, into = c("full", "region"), sep=" - ")%>% 
  separate(short, into = c("short", "R"), sep="_")%>%
  mutate(short=gsub("23","2/3", short))

res=left_join(readRDS("analysis/deseq_rhyth/WT/res.rds"), 
              readRDS("analysis/deseq_rhyth/WT/coefs.rds")  %>% 
  bind_rows(.id="cluster")) %>%
  mutate(
    phi = atan2(t_s, t_c) %% (2 * pi),
    amp = sqrt(t_c ^ 2 + t_s ^ 2),
    rel_amp = amp / log2(baseMean + 1),
    phi_hr = (12 / pi) * phi
  ) %>% dplyr::filter(grepl("^Cortex", cluster))%>%
  
  separate(cluster, into = c("cluster", "region"), sep=" - ")%>%
  
  mutate(cluster=plyr::mapvalues(x=cluster, 
                                 from=short_namesdf$full,
                                 short_namesdf$short),
         region=str_to_title(region), 
         region=factor(region,
                       levels=rev(str_to_title(
                         unique(short_namesdf$region)))))
res=res %>% pivot_wider(id_cols = c("cluster","gene"), names_from = "region", values_from = c(-cluster,-gene
                                                                                              ))

res$rel_amp_l2fc_mc=res$amp_Medial-res$amp_Caudal

res$phi_dif_mc=(12/pi)*atan2(sin(res$phi_Medial-res$phi_Caudal),cos(res$phi_Medial-res$phi_Caudal))


res$rel_amp_l2fc_rc=res$amp_Rostral-res$amp_Caudal

res$phi_dif_rc=(12/pi)*atan2(sin(res$phi_Rostral-res$phi_Caudal),cos(res$phi_Rostral-res$phi_Caudal))


res$rel_amp_l2fc_rm=res$amp_Rostral-res$amp_Medial

res$phi_dif_rm=(12/pi)*atan2(sin(res$phi_Rostral-res$phi_Medial),cos(res$phi_Rostral-res$phi_Medial))


res=res %>% pivot_longer(
  cols = c("rel_amp_l2fc_mc", "phi_dif_mc", 
           "rel_amp_l2fc_rc", "phi_dif_rc", 
           "rel_amp_l2fc_rm", "phi_dif_rm"),
  names_to = c("metric", "condition"),
  names_pattern = "(.+?)_(mc|rc|rm)",
  values_to = "value"
) %>%
  pivot_wider(
    names_from = "metric",
    values_from = "value"
  )

res$condition=str_to_upper(res$condition)

interaction_res=readRDS("analysis/deseq_region_interaction/all_deseqs_WT.rds")

interaction_res=interaction_res[grep("^[^_]*_[^_]*_[^_]*_[^_]*$", names(interaction_res))]

interaction_res=lapply(names(interaction_res), function(i) {
  x=str_split(i, "_")[[1]]
  
  if(x[1]!=x[3] | x[2]==x[4]) { return(NULL)}
  
  df=interaction_res[[i]]$res

  clust=x[1]
  comp=paste0(sort(x[c(2,4)])%>% rev, collapse = "")
  df$cluster=gsub("23", "2/3",clust)
  df$condition=comp
  df
  
  })

interaction_res=bind_rows(interaction_res) %>%
  dplyr::select(gene, cluster, condition, baseMean, padj)


res=left_join(interaction_res,res  )
res=res %>%
  dplyr::select(cluster,
                gene,
                condition,
                baseMean,
                padj, 
                phi_dif, 
                rel_amp_l2fc)

# res=res %>% 
#   mutate(condition=plyr::mapvalues(condition,
#    c("MC", "RM","RC"),
#    c("Medial vs. Caudal",
#      "Rostral vs. Medial",
#      "Rostral vs. Caudal")))
res=res %>% dplyr::filter(cluster=="L2/3")

p1=ggplot(res %>% dplyr::filter(padj<0.1),
          aes(y = rel_amp_l2fc, x =phi_dif%%24)) +
  geom_point_rast(data=res %>%na.omit() %>% dplyr::filter(padj>=0.1),
                  aes(y = rel_amp_l2fc, x =phi_dif%%24),
                  col="grey70",alpha=0.3,
                  inherit.aes = F,size=.1)+
  geom_point(aes(col=cluster), size=.3) + 
  theme_bw(base_size = 16) +
  
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.2, colour="black") + 
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.2, colour="black")+
  scale_x_continuous(breaks = 0:3*6, limits = c(0, 24)) +
 scale_y_continuous(limits=c(-0.5,.5), oob=scales::squish, breaks = pretty_breaks(n=3))+
  coord_polar(start =0)+
  facet_grid(.~condition, switch = "y")+
  ggsci::scale_color_aaas()+
  theme_minimal(base_size = 6, base_family = "ArialMT")+
  theme(
    plot.title = element_text(
      angle = 0,
      family = "ArialMT",
      size = 6,
      face = "bold",
      vjust = 1,
      hjust = .5
    ),
    
    plot.caption = element_text(size = 7, family = "ArialMT",
                                hjust =0,  margin = margin(t = -2, b = 0, unit = "pt")),
    
    legend.text = element_text(size = 6, family = "ArialMT"),
    strip.background = element_blank(),
    strip.text = element_text(size = 6, family = "ArialMT"),
    strip.placement = "outside",
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6, family = "ArialMT",color="black"),
    panel.grid.minor.y = element_blank(),
    legend.key.size = unit(0.5, "lines"),
    plot.margin = margin(t = 1, r = 1, b = 0.2, l = 1, unit = "lines")
    
  )

p1
pdf("~/desp1/precast/precast_final_with_ros_caud/figures/figure2/plots/polar_volcano_RMC.pdf",height = 2, width =6 )
plot(p1)
dev.off()
