library(tidyverse)
library(ppcor)

setwd("~/desp1/precast/prec_c25q25g3000/")
region=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_order2.rds")
fdr_cut=0.1
short_names=readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")
res1=left_join(readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/young/coefs_y.rds"),
  readRDS("deseq_rhth_int2/geno_int_new/young/interaction_results_readjusted_y.rds") %>%
    dplyr::select(gene,cluster, fdr)) %>%
  mutate(fdr=replace_na(fdr,1))

res2=left_join( readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/age_int_new/coefs_wt.rds"),
  readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/age_int_new/interaction_results_readjusted.rds")%>% dplyr::select(gene,cluster, fdr)) %>%
  mutate(fdr=replace_na(fdr,1)) 

res3=left_join(res1, res2, by=c("cluster", "gene")) 

amp_ef_g=readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/comparison_with_age_effect/all_genes.rds")

res3=left_join(res3,amp_ef_g) %>% na.omit() 

res3=res3 %>% mutate(sig=
                       case_when(fdr.x<fdr_cut &fdr.y<fdr_cut ~"Both",
                                 fdr.x<fdr_cut & fdr.y>=fdr_cut ~"Genotype",
                                 fdr.x>=fdr_cut & fdr.y<fdr_cut ~ "Age",
                                 .default = NA_character_)) %>% na.omit

#saveRDS(res3 %>% dplyr::select(cluster, gene, sig), "~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/comparison_with_age_effect/young_app_genes.rds")

df <- tibble(
  gene       = res3$gene,
  cluster   = res3$cluster,
  sig       = res3$sig,
  amp_WT7    = res3$y_amp,
  amp_WT14   = res3$o_amp,
  amp_APP7   = res3$app_amp
)

# 2) Compute raw contrasts
df2 <- df %>%
  mutate(
    delta_age  = amp_WT14  - amp_WT7,  # aging effect in NTG
    delta_geno = amp_APP7  - amp_WT7   # APP23 vs WT effect at 7mo
  )

fit_age  <- lm(delta_age  ~ 0+amp_WT7,  data = df2)
fit_geno <- lm(delta_geno ~ 0+amp_WT7,  data = df2)

df3 <- df2 %>%
  mutate(
    resid_age  = residuals(fit_age),
    resid_geno = residuals(fit_geno),
    region=plyr::mapvalues(cluster, region, names(region))
  )

# 4) Correlate the residuals
cor_res <- cor.test(df3$resid_age, df3$resid_geno)
print(cor_res)
df3$sig=factor(df3$sig, levels= c("Age","Both" ,   "APP23 7 mo." ,     "APP23 14 mo."))

# 5) Scatterplot of those residuals
p=ggplot(df3, aes(x = resid_age, y = resid_geno, colour = sig, group = sig)) +
  geom_point(alpha = 0.6, size=0.2) +
  #geom_smooth(aes(group = sig),method = "lm", color = "firebrick", se = FALSE) +
  sm_statCorr( fullrange = T,
    linetype = "dashed",text_size = 6*5/14
  )+theme_classic(base_size = 6, base_family = "ArialMT")+
  labs(
    x = "Age effect residual (14 mo–7 mo)",
    y = "Genotype effect residual (APP23–WT at 7 mo)",
    title = "Residualized comparison of aging vs. APP23 amplitude changes",
    subtitle = paste0("partial r = ", round(cor_res$estimate, 2),
                      ", p = ", signif(cor_res$p.value, 2))
  )+guides(col=guide_legend(title ="Differentially rhythmic by:"))+
  geom_abline(intercept = 0, slope = 1)+
  coord_fixed()+ggsci::scale_color_bmj()+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+ guides(
    colour = guide_legend(
      title = "Differentially rhythmic by:",
      override.aes = list(
        linetype = 0,   # remove any lines
        shape    = 16   # solid dot
      )
    )
  )

saveRDS(p, "~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/comparison_with_age_effect/plot_amp_youngapp.rds")
