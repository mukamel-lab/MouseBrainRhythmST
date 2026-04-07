library(tidyverse)
library(scales)

# Load once
coefs_b <- readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/joint/coefs_b.rds")
counts_b <- readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/joint/counts_b.rds")

# Helpers
mean_logcpm <- function(x) log2(mean(2^x))
std_dev <- function(x) {
  m <- mean_logcpm(x)
  s <- sd(x)
  list(y = m, ymin = m - s, ymax = m + s)
}

plot_gene_cluster <- function(gene, cluster) {
  # 1) Subset coefficients for that gene+cluster
  cf <- coefs_b %>%
    filter(gene == !!gene, cluster == !!cluster)
  if (nrow(cf) != 1) stop("Expected exactly one row of coefficients for that gene+cluster.")
  
  # 2) Subset observed counts
  obs <- counts_b %>%
    filter(gene == !!gene, cluster == !!cluster) %>%
    mutate(
      time = as.numeric(str_remove(time, "ZT")),
      combo = paste(genotype, age, sep = "_")
    )
  samples_df = obs %>% dplyr::rename(sample=sample_name)%>%
    distinct(sample, genotype,sex,age) 
  
  
  # 3) Build prediction grid
  pred <- expand.grid(
    sample = samples_df$sample,
    time   = seq(0, 24, length.out = 100)
  ) %>%
    left_join(samples_df, by = "sample")  %>%
    mutate(
      combo   = paste(genotype, age, sep = "_"),
      time_pi = time * 2*pi/24,
      pred_lin = pmap_dbl(list(genotype, age, time_pi, sex), function(gno, ag, tp,sx) {
        # intercept + main effects
        base <- cf$Intercept +
          if_else(gno == "WT", cf$genotype_WT_vs_APP23, 0) +
          if_else(ag == "7 months", cf$age_7.months_vs_14.months, 0) +
          if_else(sx == "M", cf$sex_M_vs_F, 0)
        # sine/cosine terms
        ts <- cf$t_s +
          if_else(ag == "7 months", cf$age7.months.t_s, 0) +
          if_else(gno == "WT", cf$genotypeWT.t_s, 0)
        tc <- cf$t_c +
          if_else(ag == "7 months", cf$age7.months.t_c, 0) +
          if_else(gno == "WT", cf$genotypeWT.t_c, 0)
        
        base + tc * cos(tp) + ts * sin(tp)
      }),
      pred_log2 = log2(2^pred_lin + 1)
    ) %>%
    group_by(genotype, age, combo, time) %>%
    summarise(pred_mean = mean_logcpm(pred_log2), .groups="drop") %>%
    group_by(age,genotype,combo) %>%
    mutate(pred_mean = pred_mean-mean(pred_mean), .groups="drop")
  
  # 4) Plot
  ggplot() +
    # observed means ± SD
    # geom_point(data = obs,
    #            aes(x = time, y = l2expr, color = combo),
    #            stat = "summary", fun = mean_logcpm, size = 1) +
    # stat_summary(data = obs,
    #              aes(x = time, y = l2expr, color = combo),
    #              fun.data = std_dev, geom = "errorbar", width = 0.5) +
    # 
    geom_line(data = pred,
              aes(x = time, y = pred_mean, color = combo),
              size = 0.8) +
    # styling
   scale_color_manual( name = "Genotype × Age",
  values = c(
    # WT in blues
    "WT_14 months"  = "#1f78b4",   # darker blue
    "WT_7 months"   = "#a6cee3",   # lighter blue
    # APP23 in reds
    "APP23_14 months" = "#e31a1c", # darker red
    "APP23_7 months"  = "#fb9a99"  # lighter red
  ),
  breaks=c("WT_7 months"  ,   # darker blue
           "WT_14 months"  ,   # lighter blue
           # APP23 in reds
           "APP23_7 months", # darker red
           "APP23_14 months"   ),
  labels = c(
    "WT, 7 mo",
    "WT, 14 mo",
    "APP23, 7 mo",
    "APP23, 14 mo"
  )
  ) +
    scale_x_continuous(breaks = seq(0,24,6), labels = seq(0,24,6)) +
    labs(
      x = "Zeitgeber time (h)",
      y = "Log₂ CPM",
      title = paste0(gene, " in ", cluster),
      subtitle = "Colored by genotype_age (NTG vs APP23 × 14 mo vs 7 mo)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position  = "bottom",
      panel.grid.minor = element_blank()
    )
}

# Example:
plot_gene_cluster("Per2", "Cortex Layer 2/3")

