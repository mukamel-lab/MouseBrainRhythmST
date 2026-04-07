library(tidyverse)
library(vroom)
library(grid)    # for unit()
library(ggpp)    # for arrow()
library(ggsci)   # for the color scale

setwd("~/desp1/precast/prec_c25q25g3000")

# 1) Load your interaction results and region mapping
coefs_y <- readRDS("deseq_rhth_int2/geno_int_new/young/coefs_y.rds") %>% dplyr::filter(padj<0.05)
coefs_o <- readRDS("deseq_rhth_int2/geno_int_new/old/coefs_o.rds") %>% dplyr::filter(padj<0.05)
cluster_order <- readRDS("objects/cluster_order.rds")


region_map    <- readRDS("objects/region_order.rds")

coefs_y <- coefs_y %>%
  mutate(region = plyr::mapvalues(cluster, cluster_order, region_map)) 
coefs_o <- coefs_o %>%
  mutate(region = plyr::mapvalues(cluster, cluster_order, region_map))

# 2) Read in TF target lists from your ChIP‐atlas files
chip_tfs=readRDS("~/desp1/precast/precast_final_with_ros_caud/analysis/chipatlas/chiptargets.rds")
tf_targets_list=lapply(chip_tfs,function(x) {x$Target_genes})

# 3) Plotting function for one TF in Cortex
plot_tf_segments <- function(tf_name) {
  targets <- tf_targets_list[[tf_name]]
  
  # helper to summarize amplitude & phase per age
  summarize_age <- function(coefs, age_label) {
    coefs %>%
      filter(region == "Cortex", gene %in% targets) %>%
      group_by(region) %>%
      summarise(
        ts     = mean(t_s),
        tc     = mean(t_c),
        ts_wt  = mean(t_s  + genotypeWT.t_s),
        tc_wt  = mean(t_c  + genotypeWT.t_c),
        .groups = "drop"
      ) %>%
      mutate(age = age_label)
  }
  
  df7  <- summarize_age(coefs_y, "7 months")
  df14 <- summarize_age(coefs_o, "14 months")
  
  df <- bind_rows(df7, df14) %>%
 
    mutate(
      amp_APP = sqrt(ts^2 + tc^2),
      phi_APP = (12/pi) * ( atan2(ts,  tc) %%2*pi ),
      amp_WT  = sqrt(ts_wt^2 + tc_wt^2),
      phi_WT  = (12/pi) * ( atan2(ts_wt,tc_wt) %% 2*pi)
    ) %>%
    # pivot into long form
    pivot_longer(
      cols = c(starts_with("amp_"), starts_with("phi_")),
      names_to  = c(".value", "geno"),
      names_sep = "_"
    ) %>%
    # make a combo factor for coloring
    mutate(
      combo = paste0(geno, "_", age),
      combo = factor(
        combo,
        levels = c(
          "WT_7 months", "WT_14 months",
          "APP_7 months", "APP_14 months"
        )
      )
    )
  full_age <- function(coefs, age_label) {
    coefs %>%
      filter(region == "Cortex", gene %in% targets) %>%
    
      mutate(age = age_label)
  }
  
  df_full= rbind(full_age(coefs_y, "7 months"),
                 full_age(coefs_o, "14 months")) %>%
    pivot_longer(
      cols = c(ends_with("_amp"), ends_with("_phi")),
      names_to  = c("geno",".value"),
      names_sep = "_"
    ) %>%
    # make a combo factor for coloring
    mutate(
      geno=str_to_upper(geno),
      combo = paste0(str_to_upper(geno), "_", age),
      combo = factor(
        combo,
        levels = c(
          "WT_7 months", "WT_14 months",
          "APP_7 months", "APP_14 months"
        )
      )
    )
  
  
  # 4) Build the polar‐segment plot
  ggplot(df, aes(x = phi, xend = phi, y = 0, yend = amp, color = combo)) +
    geom_segment(
      arrow = arrow(length = unit(0.1, "cm")),
      linewidth = 0.3
    ) +
    #geom_point(data=df_full, aes(x=12*phi/pi, y=amp, color=combo), inherit.aes = F,size=.6)+
    coord_polar(start = 0) +
    scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(0, 24))+
    #scale_y_continuous(trans="log1p")+
    scale_color_manual(
      name   = "Genotype × Age",
      values = c(
        "WT_7 months"    = "#a6cee3",
        "WT_14 months"   = "#1f78b4",
        "APP_7 months" = "#fb9a99",
        "APP_14 months"= "#e31a1c"
      ),
      labels = c(
        "WT, 7 mo", "WT, 14 mo",
        "APP23, 7 mo", "APP23, 14 mo"
      )
    ) +
    labs(title = tf_name) +
    theme_minimal(base_size = 6, base_family = "ArialMT") +
    theme(
      axis.title      = element_blank(),
     
      strip.text       = element_text(face = "italic"),
      legend.position = "bottom"
    )
}

# # 5) Generate two separate plots
# p_Tfeb   <- plot_tf_segments("Tfeb")
# p_Crebbp <- plot_tf_segments("Crebbp")
# 
# # Display them
# print(p_Tfeb)
# print(p_Crebbp)
# tfs <- c("Crebbp", "Crtc3", "Cry1", "Foxc1", "Tfeb")

# Generate a list of plots via your helper
plots <- lapply(tfs, plot_tf_segments)
plots[[5]]#+facet_grid(age~geno, as.table = F)
# Arrange them in a single row
combined_plot <- wrap_plots(plots, nrow = 1)


#(combined_plot+plot_layout(guide="collect"))&theme(legend.position = "bottom")
