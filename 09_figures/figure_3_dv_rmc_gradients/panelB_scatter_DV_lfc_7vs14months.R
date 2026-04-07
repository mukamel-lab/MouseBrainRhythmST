library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(ggh4x)  

setwd("~/desp1/precast/precast_final_with_RMC_DV/")

res=readRDS("analysis/deseq_APP_main_effect/res.rds")

arial_theme=readRDS("arial6_theme.rds")

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

ctx_wide=ctx %>% ungroup() %>% pivot_wider(id_cols = c("gene","cluster","age"), names_from = "region", values_from = c("padj", "log2FoldChange"))
hip_wide=hip %>% ungroup() %>% pivot_wider(id_cols = c("gene","cluster", "age"), 
                                           names_from = "region",
                                           values_from = c("padj", "log2FoldChange"))


#################################
#CTX Plots
################################


df_filtered = ctx_wide %>%
  filter(padj_Caudal < 0.1 | padj_Medial < 0.1 | padj_Rostral < 0.1) %>%
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
  ) %>%  dplyr::filter(if_all(starts_with("log2"), ~ !is.na(.)))

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

#  Medial
p1 = ggplot(df_filtered, aes(x = log2FoldChange_Caudal, y = log2FoldChange_Medial, color = sig_label)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.3) +
  facet_grid(age~cluster,switch = "y") +
  labs(x = "Caudal Log2FC (APP23/NTG)", y = "Medial Log2FC (APP23/NTG)", title = "Caudal vs Medial") +
  base_theme+scale_x_continuous(limits = c(-2,2), breaks = c(-2,0,2))+scale_y_continuous(limits = c(-2,2), breaks = c(-2,0,2))+
  theme(strip.placement = "outside",
        strip.text = element_text(size = 6, family = "ArialMT", colour = "black"))+coord_fixed()

# Caudal vs Rostral
p2 = ggplot(df_filtered, aes(x = log2FoldChange_Caudal, y = log2FoldChange_Rostral, color = sig_label)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.3) +
  facet_grid(age~cluster,switch = "y") +
  labs(x = "Caudal Log2FC (APP23/NTG)", y = "Rostral Log2FC (APP23/NTG)", title = "Caudal vs Rostral") +
  base_theme+scale_x_continuous(limits = c(-2,2), breaks = c(-2,0,2))+scale_y_continuous(limits = c(-2,2), breaks = c(-2,0,2))+
  theme(strip.placement = "outside",
        strip.text.x=element_blank(),
        strip.text.y = element_text(size = 6, family = "ArialMT", colour = "black"))+coord_fixed()

#Medial vs Rostral
p3 = ggplot(df_filtered, aes(x = log2FoldChange_Medial, y = log2FoldChange_Rostral, color = sig_label)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.3) +
  facet_grid(age~cluster,switch = "y") +
  labs(x = "Medial Log2FC (APP23/NTG)", y = "Rostral Log2FC (APP23/NTG)", title = "Medial vs Rostral") +
  base_theme+scale_x_continuous(limits = c(-2,2), breaks = c(-2,0,2))+scale_y_continuous(limits = c(-2,2), breaks = c(-2,0,2))+
  theme(strip.placement = "outside",
        strip.text.x=element_blank(),
        strip.text = element_text(size = 6, family = "ArialMT", colour = "black"))+coord_fixed()

final_plot = ((p1 / p2/ p3)+ 
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") & 
  theme(legend.position = "none", legend.title = element_blank()))&
  ggsci::scale_color_jama()

pdf("figures/figure3/plots/scatter_cortex_RMC.pdf", height = 4,width = 4)

print(final_plot)
dev.off()

#####################
# HIP plots
###########################
df_hip_filtered <- hip_wide %>%
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
p_dv <- ggplot(df_hip_filtered, aes(x = log2FoldChange_Dorsal, y = log2FoldChange_Ventral, color = sig_label)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.3) +
  facet_grid(age ~ cluster, switch = "y") +
  labs(
    x = "Dorsal Log2FC (APP23/NTG)",
    y = "Ventral Log2FC (APP23/NTG)",
    title = "Dorsal vs Ventral"
  ) +
  base_theme +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  scale_y_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  theme(
    strip.placement = "outside",
    strip.text = element_text(size = 6, family = "ArialMT", colour = "black")
  ) +
  coord_fixed()+ggsci::scale_color_igv()

##########
pdf("figures/figure3/plots/scatter_hippocampus_DV.pdf", height = 2, width = 4)
print(p_dv)
dev.off()
############

###################################################
# Check up/downregulate concordance for box plots
# for whichever regions a gene is significant in (FDR<0.1)
# the log2FC directions should concordant
# throws error if not
######################################################

#check ctx

stopifnot({ local({
df_ctx_concordance = ctx_wide %>%
  filter(padj_Caudal < 0.1 | padj_Medial < 0.1 | padj_Rostral < 0.1) %>%
  mutate(across(starts_with("padj_"), ~ replace_na(., 1))) %>%  # treat NA padj as nonsig
  mutate(across(starts_with("log2"), ~ replace_na(., 0))) %>%   # optional, to be safe
  rowwise() %>%
  mutate(
    
    min_region = c("Caudal", "Medial", "Rostral")[which.min(c(padj_Caudal, padj_Medial, padj_Rostral))],
    
    main_direction = case_when(
      min_region == "Caudal" ~ sign(log2FoldChange_Caudal),
      min_region == "Medial" ~ sign(log2FoldChange_Medial),
      min_region == "Rostral" ~ sign(log2FoldChange_Rostral)
    ),
    
    sig_Caudal = padj_Caudal < 0.1,
    sig_Medial = padj_Medial < 0.1,
    sig_Rostral = padj_Rostral < 0.1,
    
    # Direction in each sig region
    caudal_dir = ifelse(sig_Caudal, sign(log2FoldChange_Caudal), NA),
    medial_dir = ifelse(sig_Medial, sign(log2FoldChange_Medial), NA),
    rostral_dir = ifelse(sig_Rostral, sign(log2FoldChange_Rostral), NA)
  ) %>%
  mutate(
    sig_directions = list(na.omit(c(caudal_dir, medial_dir, rostral_dir))),
    unique_dirs = length(unique(sig_directions)),
    concordant = unique_dirs == 1
  ) %>%
  ungroup()
length(unique(df_ctx_concordance$concordant))==1})
})

# check hip

stopifnot({ local({
df_hip_concordance = hip_wide %>%
  filter(padj_Ventral < 0.1 | padj_Dorsal < 0.1) %>%
  mutate(
    sig_Ventral = padj_Ventral < 0.1,
    sig_Dorsal = padj_Dorsal < 0.1,
    sig_label = case_when(
      sig_Ventral & sig_Dorsal ~ "Both",
      sig_Ventral ~ "Ventral",
      sig_Dorsal ~ "Dorsal"
    ),
    
    # direction
    ventral_dir = ifelse(sig_Ventral, sign(log2FoldChange_Ventral), NA),
    dorsal_dir  = ifelse(sig_Dorsal,  sign(log2FoldChange_Dorsal),  NA),
    
    concordant = case_when(
      sig_Ventral & sig_Dorsal ~ (ventral_dir == dorsal_dir),
      TRUE ~ TRUE
    )
  ) %>%
  filter(if_all(starts_with("log2"), ~ !is.na(.)))
length(unique(df_hip_concordance$concordant))==1 })
})

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
      min_region == "Caudal"  & log2FoldChange_Caudal  > 0 ~ "Up",
      min_region == "Caudal"  & log2FoldChange_Caudal  < 0 ~ "Down",
      min_region == "Medial"  & log2FoldChange_Medial  > 0 ~ "Up",
      min_region == "Medial"  & log2FoldChange_Medial  < 0 ~ "Down",
      min_region == "Rostral" & log2FoldChange_Rostral > 0 ~ "Up",
      min_region == "Rostral" & log2FoldChange_Rostral < 0 ~ "Down"))
    
 


#  Medial
p1b = ggplot(df_filtered, aes(x = cluster, y = abs(log2FoldChange_Medial) - abs(log2FoldChange_Caudal))) +
  geom_point(aes(color = sig_label), size = 0.02, position = position_jitter(width = 0.2, height = 0)) +
  geom_boxplot(color = "black", fill = NA, linewidth = 0.3, outlier.shape = NA, width = 0.3) +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "gray30", linewidth = 0.3) +  
  facet_grid(~age, switch = "y") +
  labs(x = "", y = "Medial - Caudal", title = "") +
  base_theme +
  theme(
    strip.placement = "outside",
    strip.text = element_text(size = 6, family = "ArialMT", colour = "black"),
    legend.position = "bottom",
    axis.text.x = element_blank()
  ) +
  scale_y_continuous(limits = c(-.6, .6), breaks = c(-.5, 0, .5), oob = scales::squish)

# Caudal vs Rostral
p2b =ggplot(df_filtered, aes(x = cluster, y = abs(log2FoldChange_Rostral) - abs(log2FoldChange_Caudal))) +
  geom_point(aes(color = sig_label), size = 0.02, position = position_jitter(width = 0.2, height = 0)) +
  geom_boxplot(color = "black", fill = NA, linewidth = 0.3, outlier.shape = NA, width = 0.3) +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "gray30", linewidth = 0.3) +  
  facet_grid(~age, switch = "y") +
  labs(x = "", y = "Rostral - Caudal", title = "") +
  base_theme +
  theme(
    strip.placement = "outside",
    strip.text = element_text(size = 6, family = "ArialMT", colour = "black"),
    legend.position = "bottom",
    axis.text.x = element_blank()
  ) +
  scale_y_continuous(limits = c(-.6, .6), breaks = c(-.5, 0, .5), oob = scales::squish)

  
  


#Medial vs Rostral
p3b = ggplot(df_filtered, aes(x = cluster, y = abs(log2FoldChange_Rostral) - abs(log2FoldChange_Medial))) +
  geom_point(aes(color = sig_label), size = 0.02, position = position_jitter(width = 0.2, height = 0)) +
  geom_boxplot(color = "black", fill = NA, linewidth = 0.3, outlier.shape = NA, width = 0.3) +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "gray30", linewidth = 0.3) +  
  facet_grid(~age, switch = "y") +
  labs(x = "", y = "Rostral - Medial", title = "") +
  base_theme +
  theme(
    strip.placement = "outside",
    strip.text = element_text(size = 6, family = "ArialMT", colour = "black"),
    legend.position = "bottom"
  ) +
  scale_y_continuous(limits = c(-.6, .6), breaks = c(-.5, 0, .5), oob = scales::squish)

final_plotb = (((p1b / p2b/ p3b)+ 
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") & 
  theme(legend.position = "right", legend.title = element_blank()))&
    guides(color = guide_legend(override.aes = list(size = 1.5))))&
  ggsci::scale_color_jama()

pdf("figures/figure3/plots/box_cortex_RMC_log2FC_diff.pdf", height = 4,width =4.5)

print(final_plotb)
dev.off()

#####################
# HIP plots
###########################
df_hip_filtered <- hip_wide %>%
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


p_dvb=ggplot(df_hip_filtered, aes(x = cluster, y = abs(log2FoldChange_Dorsal)- abs(log2FoldChange_Ventral))) +
  geom_point(aes(color = sig_label), size = 0.02, position = position_jitter(width = 0.2, height = 0)) +
  geom_boxplot(color = "black", fill = NA, linewidth = 0.3, outlier.shape = NA, width = 0.3) +
  geom_hline(yintercept = c(-1.3, 1.3), linetype = "dotted", color = "gray30", linewidth = 0.3) +  
  facet_grid(~age, switch = "y") +
  labs(x = "", y = "Dorsal - Ventral", title = "") +
  base_theme +
  theme(
    strip.placement = "outside",
    strip.text = element_text(size = 6, family = "ArialMT", colour = "black"),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(-1.3, 1.3), breaks = c(-1, 0, 1), oob = scales::squish)+
  ggsci::scale_color_igv()



##########
pdf("figures/figure3/plots/box_hippocampus_DV_log2FC_diff.pdf", height = 1.3, width = 2.7)
print(p_dvb)
dev.off()
############


