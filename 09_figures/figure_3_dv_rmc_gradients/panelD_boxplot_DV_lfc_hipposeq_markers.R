library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(ggh4x)  

setwd("~/desp1/precast/precast_final_with_RMC_DV/")

res=readRDS("analysis/deseq_APP_main_effect/res.rds")

arial_theme = readRDS("~/desp1/precast/precast_final_with_RMC_DV/arial6_theme.rds")

base_theme = theme_bw(base_size = 6, base_family = "ArialMT") +
  arial_theme +
  theme(
    text=element_text(size = 6, family = "ArialMT", colour = "black"),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),  
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    legend.position = "bottom", 
    plot.title = element_blank()
  )
cluster_color=readRDS("cluster_color.rds")
short_names = readRDS("short_names.rds")
short_namesdf=data.frame(full=short_names,short=names(short_names))%>%
  separate(full, into = c("full", "region"), sep=" - ")%>%
  separate(short, into = c("short", "R"), sep="_")%>%
  mutate(short=gsub("23","2/3", short))

res=res %>%
  mutate(age=factor(age, levels=paste0(c(7, 14), " months")))%>%
  ungroup() %>%left_join(cluster_color) %>%
  mutate(cluster=str_split_fixed(cluster, " - ",2)[,1],
                           cluster=plyr::mapvalues(cluster,
                                                   short_namesdf$full, 
                                                   short_namesdf$short))
res$region=str_split_fixed(res$region, " - ", 2)[,2]

res=res %>% dplyr::filter(region!="")
ctx=res %>% dplyr::filter(region %in% c("Rostral", "Medial", "Caudal"))
hip=res %>% dplyr::filter(cluster %in% c("CA3so" ,"CA3sp", "DGmo" , "DGsg"))

ctx_wide=ctx %>% ungroup() %>% pivot_wider(id_cols = c("gene","cluster","age"), 
                                           names_from = "region",
                                           values_from = c("padj", 
                                                           "log2FoldChange"))
hip_wide=hip %>% ungroup() %>% pivot_wider(id_cols = c("gene","cluster", "age"), 
                                           names_from = "region",
                                           values_from = c("padj", "log2FoldChange"))


###################################################################
#BOX plots
##################################################################

#################################
#CTX Plots
################################


df_filtered = ctx_wide %>%
  filter(padj_Caudal < 0.1 | padj_Medial < 0.1 |padj_Rostral < 0.1 ) %>%  
   dplyr::filter(if_all(starts_with("log2"), ~ !is.na(.))) %>%
  mutate(
    sig_Caudal = padj_Caudal < 0.1,
    sig_Medial = padj_Medial < 0.1,
    sig_Rostral = padj_Rostral < 0.1,
    sig_label = case_when(
      sig_Caudal & sig_Medial & sig_Rostral ~ "Both",
      sig_Caudal & sig_Medial ~ "Both",
      sig_Caudal & sig_Rostral ~ "Both",
      sig_Medial & sig_Rostral ~ "Medial and/or\nRostral",
      sig_Medial ~ "Medial and/or\nRostral",
      sig_Rostral ~ "Medial and/or\nRostral",
      sig_Caudal ~ "Caudal"
      
    )
  ) %>%
  mutate(
    across(starts_with("padj_"), ~ replace_na(., 1)),

    min_region = c("Caudal", "Medial", "Rostral")[apply(select(., padj_Caudal, padj_Medial, padj_Rostral), 1, which.min)],
    
    direction = case_when(
      min_region == "Caudal"  & log2FoldChange_Caudal  > 0 ~ "Upregulated",
      min_region == "Caudal"  & log2FoldChange_Caudal  < 0 ~ "Downregulated",
      min_region == "Medial"  & log2FoldChange_Medial  > 0 ~ "Upregulated",
      min_region == "Medial"  & log2FoldChange_Medial  < 0 ~ "Downregulated",
      min_region == "Rostral" & log2FoldChange_Rostral > 0 ~ "Upregulated",
      min_region == "Rostral" & log2FoldChange_Rostral < 0 ~ "Downregulated")) 
  
    
 

res_all_genes=readRDS("analysis/deseq_APP_main_effect/region_interaction/cortex/res_all_genes.rds")

int_effect=lapply(res_all_genes, 
                  function(dsq_list) {lapply(dsq_list, function(dsq) {
                    coefficients(dsq, ) %>% as.data.frame %>%
                      rownames_to_column("gene")
}) %>% bind_rows(.id="cluster")}) %>% bind_rows(.id="age")

int_effect$cluster=gsub("23","2/3", int_effect$cluster)

df_filtered2=left_join(df_filtered, int_effect )

signif_data =df_filtered2 %>%
  filter(age == "14 months") %>%
  group_by(cluster, direction) %>%
  summarise(
    pval = tryCatch(
      t.test(regionR.genotypeAPP23, mu = 0)$p.value,
      error = function(e) NA
    ),
    max_y = max(regionR.genotypeAPP23, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_position = .8*ifelse(direction=="Upregulated", 1, -1),  # space above the data
    annotation = ifelse(pval < 0.001, "***",
                        ifelse(pval < 0.01, "**",
                               ifelse(pval < 0.05, "*", "")))
  )

# Rostral vs  Caudal
p2b =ggplot(df_filtered2 %>% dplyr::filter(age=="14 months"), 
            aes(x = cluster, y = regionR.genotypeAPP23)) +
  geom_point(aes( col=direction),size = 0.02) +
  geom_boxplot( fill = NA, linewidth = 0.3, outlier.shape = NA, 
                col="black") +
  geom_hline(yintercept = c(-.75, .75), linetype = "dotted", color = "gray30", linewidth = 0.3) + 
  labs(x = "", y = "Interaction Effect (Rosral APP23/Caudal APP23)", title = "") +
  base_theme +
  theme(
    strip.placement = "outside",
    strip.text = element_text(size = 6, family = "ArialMT", colour = "black"),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(-.75, .75), breaks = c(-.5, 0, .5), 
                     oob = scales::squish)+
scale_color_manual(name="", values=ggsci::pal_nejm()(2)[2:1])+
  geom_text(
    data = signif_data,
    aes(
      x = cluster,
      group=direction,
      label=annotation,
      y = y_position
    ), 
    position = position_dodge2(padding = .3,width = 0.4), 
  
    color = "black"
  )+facet_grid(~direction)



pdf("figures/figure3/plots/interaction_effect_size_ros_v_caud.pdf",
    height = 2,width = 3.5)


p2b
dev.off() 


###################
#HIP
################

df_hip_filtered = hip_wide %>%
  filter(padj_Ventral < 0.1 | padj_Dorsal < 0.1) %>%
  mutate(
    sig_Ventral = padj_Ventral < 0.1,
    sig_Dorsal = padj_Dorsal < 0.1,
    sig_label = case_when(
      sig_Ventral & sig_Dorsal ~ "Both",
      sig_Ventral ~ "Ventral",
      sig_Dorsal ~ "Dorsal"
    )
  ) %>%  dplyr::filter(if_all(starts_with("log2"), ~ !is.na(.)))

