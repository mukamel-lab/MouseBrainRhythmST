library(tidyverse)
library(circular)
library(smplot2)  
library(ggrastr)
setwd("~/desp1/precast/prec_c25q25g3000/")

region = readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_order2.rds")

short_names = readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")
fdr_cut=0.1
# 1) Load data and filter based on FDR
res1 = left_join(
  readRDS("deseq_rhth_int2/geno_int_new/young/coefs_y.rds"),
   readRDS("deseq_rhth_int2/geno_int_new/young/res_y_aw.rds") %>%
    dplyr::select(gene, cluster, padj, genotype) %>% pivot_wider(names_from = genotype,
                                                                 values_from = padj, 
                                                                 values_fill = 1, 
                                                                 names_prefix = "padj_")
) 

res2 = left_join(
  readRDS("deseq_rhth_int2/age_int_new/coefs_wt.rds"),
  readRDS("deseq_rhth_int2/age_int_new/res_wt_yo.rds") %>%
    dplyr::select(gene, cluster, padj, age) %>% pivot_wider(names_from = age,
                                                                 values_from = padj, 
                                                                 values_fill = 1, 
                                                                 names_prefix = "padj_")
) 

res3 = left_join(res1, res2, by = c("cluster", "gene"))
res3=res3 %>%filter(
    # compare each padj to fdr_cut, coerce TRUE→1, FALSE→0, then require ≥2 trues
    ( (padj_APP23    < fdr_cut) +
      (`padj_14 months` < fdr_cut) +
      (`padj_7 months`  < fdr_cut)
    ) >= 3
  ) 

# 2) Build dataframe
df <- tibble(
  gene = res3$gene,
  cluster = res3$cluster,
  sig = res3$sig,
  phi_WT7 = res3$y_phi,    # 7mo WT phase (radians)
  phi_WT14 = res3$o_phi,   # 14mo WT phase (radians)
  phi_APP7 = res3$app_phi  # 7mo APP23 phase (radians)
)

# 3) Convert phases to circular objects
df2 <- df %>%
  mutate(
    phi_WT7  = circular(phi_WT7, units = "radians", modulo = "2pi"),
    phi_WT14 = circular(phi_WT14, units = "radians", modulo = "2pi"),
    phi_APP7 = circular(phi_APP7, units = "radians", modulo = "2pi")
  )

# 4) Calculate circular differences (WT14−WT7 and APP7−WT7)
delta_age  <- df2$phi_WT14 - df2$phi_WT7
delta_geno <- df2$phi_APP7 - df2$phi_WT7

# 5) Represent all circular variables as (sin, cos) components
df3 <- df2 %>%
  mutate(
    sin_wt7   = sin(phi_WT7),
    cos_wt7   = cos(phi_WT7),
    sin_wt14 = sin(phi_WT14),
    cos_wt14 = cos(phi_WT14),
    sin_app7= sin(phi_APP7),
    cos_app7= cos(phi_APP7),
    delta_age  = phi_WT14 - phi_WT7,
    delta_geno = phi_APP7 - phi_WT7,
    sin_delta_age=sin(delta_age),
    cos_delta_age=cos(delta_age),
    sin_delta_geno=sin(delta_geno),
    cos_delta_geno=cos(delta_geno),
  )

# 6) Regress out baseline phase (sin and cos separately)
fit_sin_age  <- lm(sin_delta_age  ~0+ sin_wt7, data = df3)
fit_cos_age  <- lm(cos_delta_age  ~0+ cos_wt7, data = df3)

fit_sin_geno <- lm(sin_delta_geno ~0+ sin_wt7 , data = df3)
fit_cos_geno <- lm(cos_delta_geno ~0+ cos_wt7, data = df3)

df3 <- df3 %>%
  mutate(
    resid_sin_age  = residuals(fit_sin_age),
    resid_cos_age  = residuals(fit_cos_age),
    resid_sin_geno = residuals(fit_sin_geno),
    resid_cos_geno = residuals(fit_cos_geno),
    region = plyr::mapvalues(cluster, region, names(region))
  )

# 7) Reconstruct residual angles
resid_angle_age  <- atan2(df3$resid_sin_age,  df3$resid_cos_age)
resid_angle_geno <- atan2(df3$resid_sin_geno, df3$resid_cos_geno)



df3 <- df3 %>% mutate (
  

  resid_age  = resid_angle_age,
  resid_geno = resid_angle_geno

  
)
df3$region=factor(df3$region, levels=names(region)[!duplicated(names(region))])

df4=split(df3, df3$region)

cors=lapply(df4,function(dfspl){
  cor.circular(dfspl$resid_age, dfspl$resid_geno, test = T)
})

stat_df <- tibble(
  region = factor(names(cors),
                  levels=names(region)[!duplicated(names(region))]),
  cor = map_dbl(cors, ~.$cor),
  p   = map_dbl(cors, ~.$p.value)
) %>%
  mutate(
    label = sprintf("r = %.3f\np = %.2g", cor, p),
    # pick a nice corner inside the panel:
    x = -pi * 0.8,
    y =  pi* 0.8
  )
cor_all=cor.circular(df_plot$resid_age, df_plot$resid_geno, test = T)
stat_df2 <- tibble(
  cor = cor_all$cor,
  p   =cor_all$p.value
) %>%
  mutate(
    label = sprintf("r = %.3f\np = %.2g", cor, p),
    # pick a nice corner inside the panel:
    x = -pi * 0.8,
    y =  pi * 0.8
  )
pdf("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/figures/phase_shared.pdf", height = 2.8, width = 5)

ggplot(df3, aes(x = resid_age, y = resid_geno)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_density_2d(size = .2) +
  geom_point(size=0.1, col="gray30")+# geom_scattermore(col="gray30",
  #                  pointsize   = 11.2,
  #                  alpha = .7,
  #                  pixels      = c(1500,1500),
  #                  interpolate =F ) +
  facet_grid(~ region) +
  geom_label(
    data        = stat_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size        = 2,
    hjust       = 0
  ) +
  scale_x_continuous(
    breaks = seq(-pi, pi, by = pi/2),
    labels = c("-12", "-6", "0", "6", "12"),
    limits=c(-pi,pi),
    oob=scales::squish
  ) +
  scale_y_continuous(
    breaks = seq(-pi, pi, by = pi/2),
    labels = c("-12", "-6", "0", "6", "12"),
    limits=c(-pi,pi),
    oob=scales::squish
  ) +
  labs(
    x     = "Residualized age phase (rad)",
    y     = "Residualized genotype phase (rad)",
    title = "Age vs. Genotype Phase Reprogramming (2D density)",
    subtitle = "Contour lines show regions of high point density"
  ) +
  theme_bw(base_size = 6, base_family = "ArialMT")+
  theme(strip.background = element_blank())+coord_fixed()
dev.off()

