library(tidyverse)
library(ggrepel)
library(patchwork)



x=vroom::vroom("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth_int/for_share/young/results.csv.gz") %>%
  dplyr::filter(padj<0.05)

phase_relamp <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhyth2/old_young_sep/phase_relamp.rds")
clock_g=c("Dbp", "Per2", "Nr1d1", "Per3", "Nfil3", "Arntl","Cry1","Cry2", "Per1","Clock",
          "Nr1d2","Ciart", "Tef")
clock_g=c("Dbp",   "Arntl",
          "Nr1d2")

names(clock_g)=clock_g
plts=lapply(clock_g, function(gn){
#y=phase_relamp %>% dplyr::filter(paste0(gene, cluster) %in% paste0(x$gene, x$cluster))
 y=phase_relamp %>% dplyr::filter(gene%in% gn, grepl("Cortex L", cluster))

z=y %>% dplyr::select(genotype, cluster,gene,rel_amp,age, phi, amp, phi_hr) %>%
pivot_wider(id_cols =c( gene, age, cluster),
names_from = genotype, values_from = c(rel_amp, amp, phi_hr, phi) )

z$amp_dif=2*(z$amp_APP-z$amp_WT)
z$phi_dif=atan2(sin(z$phi_APP-z$phi_WT),cos(z$phi_APP-z$phi_WT))
cluster_color <- readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
temp_col=cluster_color

p2=ggplot(z %>% mutate(
  age = relevel(factor(plyr::mapvalues(
    age,
    c("young", "old"),
    paste0(1:2 * 7, " months")
  )), "7 months"),
  cluster = factor(cluster, levels = names(cluster_color))
),
aes(x = amp_dif, y = phi_dif, col = cluster)) +
  geom_point() + facet_wrap( ~ age) +
  geom_label_repel(aes(label = gene), inherit.aes = T) +
  scale_color_manual(
    values = colorspace::darken(temp_col,
                                amount = .3) %>%
      setNames(names(cluster_color)),
    name = ""
  ) + theme_bw(base_size = 16) +
  xlab("Peak-to-Trough Difference (log2(Normalized Expression))") +
  ylab("Difference in Phase (hours)") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

p2
})

dta=lapply(plts, function(df) {df$data})
dta=bind_rows(dta)
p2=ggplot(dta,
          aes(x = amp_dif, y = phi_dif, col = gene)) +
  geom_point() + facet_wrap( ~ age) +
  geom_label_repel(aes(label = gene), inherit.aes = T) +
  theme_bw(base_size = 16) +
  xlab("Peak-to-Trough Difference (log2(Normalized Expression))") +
  ylab("Difference in Phase (hours)") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

p2
 