x=vroom::vroom("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth_int/for_share/app_age/results.csv.gz") %>%
dplyr::filter(padj<0.05)
# x=rbind(x,x=vroom::vroom("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth_int/for_share/young/results.csv.gz") %>%
#           dplyr::filter(padj<0.05))
phase_relamp <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhyth2/old_young_sep/phase_relamp.rds")

y=phase_relamp %>% dplyr::filter(paste0(gene, cluster) %in% paste0(x$gene, x$cluster))

z=y %>% dplyr::select(genotype, cluster,gene,rel_amp,age, phi, amp, phi_hr) %>%
pivot_wider(id_cols =c( gene, age, cluster),
names_from = genotype, values_from = c(rel_amp, amp, phi_hr, phi) )

z$amp_dif=2*(z$amp_APP-z$amp_WT)
z$phi_dif=atan2(sin(z$phi_APP-z$phi_WT),cos(z$phi_APP-z$phi_WT))


ggplot(z %>% mutate(
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


#############################
phase_relamp <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhyth2/phase_relamp.rds")

y=phase_relamp %>% dplyr::filter(paste0(gene, cluster) %in% paste0(x$gene, x$cluster))
z=y %>% dplyr::select(genotype, cluster,gene,rel_amp, phi, amp, phi_hr) %>%
  pivot_wider(id_cols =c( gene, cluster),
              names_from = genotype, values_from = c(rel_amp, amp, phi_hr, phi) )

z$amp_dif=2*(z$amp_APP-z$amp_WT)
z$phi_dif=atan2(sin(z$phi_APP-z$phi_WT),cos(z$phi_APP-z$phi_WT))

cluster_color <- readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
temp_col=cluster_color
ggplot(z %>% mutate(cluster=factor(cluster, levels = names(cluster_color))),
       aes(x=amp_dif, y=phi_dif, col=cluster))+
  geom_point()+
  geom_label_repel(aes(label=gene), inherit.aes = T)+
  scale_color_manual(values=colorspace::darken(temp_col,
                                               amount = .3) %>%
                       setNames(names(cluster_color)), name="")+theme_bw()+
  xlab("Peak-to-Trough Difference (log2(Normalized Expression))")+ 
  ylab("Difference in Phase (hours)")+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)


