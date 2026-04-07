library(tidyverse)
library(ggrepel)
library(patchwork)



x=vroom::vroom("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth_int/for_share/full/results.csv.gz") %>%
  dplyr::filter(padj<0.05)
#  x=vroom::vroom("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth_int/for_share/young/results.csv.gz") 
phase_relamp <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhyth2/old_young_sep/phase_relamp.rds")

y=phase_relamp %>% dplyr::filter(paste0(gene, cluster) %in% paste0(x$gene, x$cluster))
#y=phase_relamp %>% dplyr::filter(gene %in% "Dbp")

z=y %>% dplyr::select(genotype, cluster,gene,rel_amp,age, phi, amp, phi_hr) %>%
pivot_wider(id_cols =c( gene, age, cluster),
names_from = genotype, values_from = c(rel_amp, amp, phi_hr, phi) )

z$amp_dif=2*(z$amp_APP-z$amp_WT)
z$rel_amp_dif=2*(z$rel_amp_APP-z$rel_amp_WT)
z$rel_amp_l2fc=log2(z$rel_amp_APP/z$rel_amp_WT)

z$phi_dif=(12/pi)*atan2(sin(z$phi_APP-z$phi_WT),cos(z$phi_APP-z$phi_WT))
cluster_color <- readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
temp_col=cluster_color

p2=ggplot(z %>% mutate(
  age = relevel(factor(plyr::mapvalues(
    age,
    c("young", "old"),
    paste0(1:2 * 7, " months")
  )), "7 months"),
  cluster = factor(cluster, levels = names(cluster_color))
) %>% dplyr::filter(cluster %in% names(cluster_color)[1:6]),
aes(x = amp_dif, y = phi_dif, col = cluster)) +
  geom_point() + facet_wrap( ~ age) +
  geom_label_repel(aes(label = gene), inherit.aes = T, show.legend = F) +
  scale_color_manual(
    values = colorspace::darken(temp_col,
                                amount = .3) %>%
      setNames(names(cluster_color)),
    name = ""
  ) + theme_bw(base_size = 16) +
  xlab("Peak-to-Trough Difference (log2(Normalized Expression))") +
  ylab("Difference in Phase (hours)") + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

p2


################
p3=ggplot(z %>% mutate(
  age = relevel(factor(plyr::mapvalues(
    age,
    c("young", "old"),
    paste0(1:2 * 7, " months")
  )), "7 months"),
  cluster = factor(cluster, levels = names(cluster_color))
) %>% dplyr::filter(cluster %in% names(cluster_color)[1:6]),
aes(x = rel_amp_dif, y = phi_dif, col = cluster)) +
  geom_point() + facet_wrap( ~ age) +
  geom_label_repel(aes(label = gene), inherit.aes = T, show.legend = F) +
  scale_color_manual(
    values = colorspace::darken(temp_col,
                                amount = .3) %>%
      setNames(names(cluster_color)),
    name = ""
  ) + theme_bw(base_size = 16) +
  xlab("Relative Amplitude Difference") +
  ylab("Difference in Phase (hours)") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p3
#####################################################
p4=ggplot(z %>% mutate(
  age = relevel(factor(plyr::mapvalues(
    age,
    c("young", "old"),
    paste0(1:2 * 7, " months")
  )), "7 months"),
  cluster = factor(cluster, levels = names(cluster_color))
) %>% dplyr::filter(cluster %in% names(cluster_color)[1:6]),
aes(y = amp_dif, x = phi_dif, col = cluster)) +
  geom_point() + facet_wrap( ~ age) +
  geom_label_repel(aes(label = gene), inherit.aes = T, show.legend = F) +
  scale_color_manual(
    values = colorspace::darken(temp_col,
                                amount = .3) %>%
      setNames(names(cluster_color)),
    name = ""
  ) + theme_bw(base_size = 16) +
  ylab("Peak-to-Trough Difference (log2(Normalized Expression))") +
  xlab("Difference in Phase (hours)") + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+coord_polar()

p4
#############################

ctx_vals=z %>% mutate(
  age = relevel(factor(plyr::mapvalues(
    age,
    c("young", "old"),
    paste0(1:2 * 7, " months")
  )), "7 months"),
  cluster = factor(cluster, levels = names(cluster_color))
) %>% dplyr::filter(cluster %in% names(cluster_color)[1:6])


ctx_vals_phi=ctx_vals %>% pivot_wider(id_cols = c("cluster","gene"),
                                      names_from = "age", 
                                      values_from = "phi_dif")
ctx_vals_rel_amp=ctx_vals %>% pivot_wider(id_cols = c("cluster","gene"),
                                      names_from = "age", 
                                      values_from = "rel_amp_dif")

ctx_vals_amp=ctx_vals %>% pivot_wider(id_cols = c("cluster","gene"),
                                          names_from = "age", 
                                          values_from = "amp_dif")
#########################################
ggplot(ctx_vals_phi,aes(x=`7 months`,y=`14 months`,col=cluster))+geom_point()+geom_abline(intercept = 0,slope = 1)+
  coord_fixed()+
  theme_bw(base_size = 16)+
  geom_label_repel(data=ctx_vals_phi %>% dplyr::filter(gene %in% x$gene[x$padj<0.05 & x$cluster %in% unique(ctx_vals_phi$cluster)]),
                   aes(x=`7 months`,y=`14 months`,col=cluster,
                       label = gene),
                   inherit.aes = F, 
                   show.legend = F) +
  scale_color_manual(
    values = colorspace::darken(temp_col,
                                amount = .3) %>%
      setNames(names(cluster_color)),
    name = ""
  )+xlab("7 month Phase Difference (hours)")+
  ylab("14 month Phase Difference (hours)")+
  ggtitle("Phase Difference FDR<0.1")

############################
ggplot(ctx_vals_amp,aes(x=`7 months`,y=`14 months`,col=cluster))+geom_point()+geom_abline(intercept = 0,slope = 1)+
  coord_fixed()+
  theme_bw(base_size = 16)+
  geom_label_repel(data=ctx_vals_amp %>% dplyr::filter(gene %in% x$gene[x$padj<0.05 & x$cluster %in% unique(ctx_vals_phi$cluster)]),
                   aes(x=`7 months`,y=`14 months`,col=cluster,
                       label = gene),
                   inherit.aes = F, 
                   show.legend = F) +
  scale_color_manual(
    values = colorspace::darken(temp_col,
                                amount = .3) %>%
      setNames(names(cluster_color)),
    name = ""
  ) +xlab("7 month Peak-to-Trough Difference")+
  ylab("14 month Peak-to-Trough Difference")+
  ggtitle("Amplitude Difference FDR<0.1")

##############################



ggplot(ctx_vals_rel_amp,aes(x=`7 months`,y=`14 months`,col=cluster))+geom_point()+geom_abline(intercept = 0,slope = 1)+
  coord_fixed()+
  theme_bw(base_size = 16)+
  geom_label_repel(aes(label = gene),
                   inherit.aes = T, 
                   show.legend = F) +
  scale_color_manual(
    values = colorspace::darken(temp_col,
                                amount = .3) %>%
      setNames(names(cluster_color)),
    name = ""
  )  +xlab("7 month Peak-to-Trough Difference")+
  ylab("14 month Peak-to-Trough Difference")+
  ggtitle("Amplitude Difference FDR<0.1")



##############################################
ctx_vals_rel_amp_lfc=ctx_vals %>% pivot_wider(id_cols = c("cluster","gene"),
                                          names_from = "age", 
                                          values_from = "rel_amp_l2fc")

ggplot(ctx_vals_rel_amp_lfc,aes(x=`7 months`,y=`14 months`,col=cluster))+geom_point()+geom_abline(intercept = 0,slope = 1)+
  coord_fixed()+
  theme_bw(base_size = 16)+
  geom_label_repel(data=ctx_vals_rel_amp_lfc %>% dplyr::filter(gene %in% x$gene[x$padj<0.05 & x$cluster %in% unique(ctx_vals_phi$cluster)]),
                   aes(x=`7 months`,y=`14 months`,col=cluster,
                       label = gene),
                   inherit.aes = F, 
                   show.legend = F) +
  scale_color_manual(
    values = colorspace::darken(temp_col,
                                amount = .3) %>%
      setNames(names(cluster_color)),
    name = ""
  ) +xlab("7 month Log2FC Rel. Amp")+
  ylab("14 month Log2FC Rel. Amp")+
  ggtitle("Log2FC Rel. Amp FDR<0.1")