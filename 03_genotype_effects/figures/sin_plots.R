library(tidyverse)
library(patchwork)
library(pbmcapply)
library(ggh4x)
library(scales)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/young/")

arial_theme=readRDS("~/desp1/precast/precast_final_with_RMC_DV/arial6_theme.rds")
short_names = readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")

counts=readRDS("counts_y.rds")
coefs=readRDS("coefs_y.rds")
counts$cluster=plyr::mapvalues(counts$cluster, short_names, names(short_names))
coefs$cluster=plyr::mapvalues(coefs$cluster, short_names, names(short_names))

mean_logcpm = function(x) {
  log2(mean(2^x))
}

std_dev = function(x) {
  mn = mean_logcpm(x)
  sdv = sd(x)
  list(y = mn, ymin = mn - sdv, ymax = mn + sdv)
}



make_plot = function(genes_to_plot, clusters_to_plot, center=T) {
  build_df=function(cluster_to_plot,gene_to_plot) {
    dat    = counts %>% dplyr::filter(cluster==cluster_to_plot, gene==gene_to_plot)
    crv    =  coefs %>% dplyr::filter(cluster==cluster_to_plot, gene==gene_to_plot)
    
    
    dat$time=as.numeric(gsub("ZT", "", dat$time))
    
    dat2 = rbind(dat,
                 dat %>%
                   mutate(time =time+24 ))
    
    
    
    samples_df = dat2 %>% dplyr::rename(sample=sample_name)%>%
      distinct(sample, genotype,sex) 
    
    pred_data = expand.grid(
      sample = samples_df$sample,
      time   = seq(0, 42, length.out = 100)
    ) %>%
      left_join(samples_df, by = "sample") %>%
      mutate(pred_lin = pmap_dbl(
        list(genotype, sex,time),
        function(gno, sx, tm) {
          val = crv$Intercept
          
          if (gno == "WT")  val = val + crv$genotype_WT_vs_APP23
          
          if (sx == "M")          val = val + crv$sex_M_vs_F
          
          time_pi = (tm %% 24) * 2*pi/24
          
          tc = crv$t_c
          ts = crv$t_s
          if (gno=="WT") {
            return(val+(tc+crv$genotypeWT.t_c)*cos(time_pi) + (ts+crv$genotypeWT.t_s)*sin(time_pi))
          } else {
          return(val + tc*cos(time_pi) + ts*sin(time_pi))}
        }
      )) %>%
      mutate(pred_log2 = log2(2^pred_lin+1))
    
    pred_avg = pred_data %>%
      group_by(genotype, time) %>%
      summarise(pred_mean = mean_logcpm(pred_log2), .groups = "drop")
   if(center) { mn=mean(pred_avg$pred_mean)
    pred_avg$pred_mean=pred_avg$pred_mean-mn
    dat2$l2expr=dat2$l2expr-mn
   }
    list(to_p=dat2, crv=pred_avg %>% dplyr::mutate(gene=gene_to_plot, cluster=cluster_to_plot))
  }
  combos=expand.grid(genes_to_plot, clusters_to_plot)
  
  dfs=pbmclapply(1:nrow(combos), function(i) {
    build_df(combos$Var2[i], combos$Var1[i])
  })
  
  dat2=lapply(dfs, function(x) {x$to_p})%>% bind_rows() 
  crv=lapply(dfs, function(x) {x$crv})%>% bind_rows()
  dat2= dat2 %>% mutate(gene=factor(gene, levels=genes_to_plot),
                        cluster=factor(cluster, levels=clusters_to_plot))
  crv=crv%>% mutate(gene=factor(gene, levels=genes_to_plot),
                    cluster=factor(cluster, levels=clusters_to_plot))
  
  plt = ggplot() +
    annotate("rect", xmin = 0,  xmax = 12, ymin = -Inf, ymax = Inf,
             fill = "#F6F18F", alpha = 0.1) +
    annotate("rect", xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf,
             fill = "#606161", alpha = 0.1) +
    annotate("rect", xmin = 24, xmax = 36, ymin = -Inf, ymax = Inf,
             fill = "#F6F18F", alpha = 0.1) +
    annotate("rect", xmin = 36, xmax = 42, ymin = -Inf, ymax = Inf,
             fill = "#606161", alpha = 0.1) +
    
    geom_point(
      data = dat2,
      aes(x = time, y = l2expr, col = genotype),
      stat = "summary", fun = mean_logcpm, size = 0.5
    ) +
    stat_summary(
      data = dat2,
      aes(x = time, y = l2expr, col = genotype),
      fun.data = std_dev,
      geom = "errorbar", width = 0.02, linewidth = 0.21
    ) +
    
    geom_line(
      data = crv,
      aes(x = time, y = pred_mean, col = genotype),
      linewidth = 0.5
    ) +
    
    scale_color_manual(name="",
                       values = ggsci::pal_nejm()(2) %>% setNames(c("APP23","WT")),
                       breaks=c("APP23","WT"),labels=c("APP23","NTG")
    )+ 
    facet_grid2(gene~cluster, switch = "y", scales = "free_y")+
    theme_bw(base_size = 6, base_family = "ArialMT") +
    arial_theme +
    theme(
      axis.title.x     = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 0, hjust = 0.5),
      axis.title.y     = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x = element_text(size=6, family = "ArialMT"),
      strip.text.y = element_text(size=8, family = "ArialMT", face="italic"),
      legend.position  = "bottom"
    ) +
    scale_x_continuous(
     # limits=c(0,24),
      breaks = c(0, 12, 24, 36),
      labels = c("0", "12", "", "")
    ) +
    guides(color = guide_legend(keywidth = 0.6, keyheight = 0.3)) +
    scale_y_continuous(breaks = pretty_breaks(n = 3))
 
  return(plt)
}



plt=make_plot(c("Per1","Per2", "Bhlhe40","Cry2"), grep("^L", coefs$cluster, value = T) %>% unique() %>% sort %>% .[1:6])



pdf("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/figures/sin_clock_cortex.pdf",height = 3, width =4)

plot(plt)
dev.off()

int_res=readRDS("interaction_results_clock_g_y.rds") %>% 
  dplyr::filter(gene %in% c("Per1","Per2", "Bhlhe40","Cry2"), grepl("^Cortex|CTX", cluster)) %>% dplyr::filter(fdr<0.1)


int_res_sum = int_res %>% group_by(gene) %>%
  mutate(wtts = t_s + genotypeWT.t_s,
         wttc = t_c + genotypeWT.t_c) %>%
  summarise(
    tc = mean(t_c),
    ts = mean(t_s),
    wtts = mean(wtts),
    wttc = mean(wttc)
  ) %>% ungroup() %>%
  mutate(
    amp_wt = sqrt(wtts ^ 2 + wttc ^ 2),
    phihr_wt = (12 / pi) * atan2(wtts, wttc) %% (2 * pi),
    amp_app = sqrt(ts ^ 2 + tc ^ 2),
    phihr_app = (12 / pi) * atan2(ts, tc) %% (2 * pi)
  ) %>%
  pivot_longer(
    cols      = c(starts_with("amp_"), starts_with("phihr_")),
    names_to  = c(".value", "genotype"),
    names_sep = "_"
  )
int_res=int_res %>% dplyr::select(-app_phi,-wt_phi) %>%
  pivot_longer(
  cols      = c(ends_with("_amp"), ends_with("_phi_hr")),
  names_to  = c( "genotype",".value"),
  names_sep = "_"
)
plts=lapply(c("Per1","Per2", "Bhlhe40","Cry2"), function(gn){
  df=int_res %>% dplyr::filter(gene==gn)
  df_sum=int_res_sum %>% dplyr::filter(gene==gn)
  brks=round(c(0, 0.5*max(df$amp))*10)/10
  
 ggplot() +
  geom_segment(data=df_sum, aes(x = phihr, xend = phihr, y = 0, yend = amp, col = genotype),
               arrow = arrow(length = unit(0.05, "npc")), linewidth = .2) +
  geom_point(data=df, aes(x=phi, amp, col=genotype), size=0.2)+
  coord_polar(start = 0) +
  scale_x_continuous(breaks = 0:3 * 6, limits = c(0, 24)) +
  scale_color_manual(name="",
                     values = ggsci::pal_nejm()(2) %>% setNames(c("app","wt")),
                     breaks=c("app","wt"),labels=c("APP23","NTG")
  )+ 
  theme_minimal() +
   arial_theme+# Removes background, axes, gridlines, titles
  theme(
    strip.text = element_blank(),         # Hides facet labels
    legend.position = "none",             # Hides legend
    panel.border = element_blank(),       # Removes panel border
    panel.spacing = unit(0.2, "lines"),
    axis.title = element_blank(),
    axis.text.x = element_blank()
    
  )+scale_y_continuous(breaks =brks)
 
})

pdf("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/figures/mrv_core_clock.pdf",height = 3, width =1)

plot(wrap_plots(plts, ncol = 1))

dev.off()


int_res=readRDS("interaction_results_readjusted_y.rds") %>% dplyr::filter(fdr<0.05, grepl("sg", cluster))

make_plot(int_res$gene, "DGsg", center = F) +facet_wrap(~gene, scales = "free")
int_res$gene[c(8,9,20,21,22)]
p=make_plot(c("Tet3", "Cxxc5", "Acsl5", "Smim20", "Kdm7a", "Tmem158", "P4ha1", 
            "Hcfc2"), "DGsg", center = F) +facet_wrap(gene~., scales = "free", ncol=2, strip.position = "left")+theme(strip.placement = "outside")+scale_y_continuous(breaks=pretty_breaks(n=3))
pdf("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/figures/sin_dg.pdf", height = 2.8, width = 2.2)
plot(p)
dev.off()


int_res=readRDS("interaction_results_readjusted_y.rds") %>% dplyr::filter(fdr<0.1, grepl("2/3", cluster))
p=make_plot(c("Midn", "Cd4", "Dusp4", "Dusp6", "Otud1", "Egr1", "P4ha1", 
              "Irs2"), "L2/3", center = F) +facet_wrap(gene~., scales = "free", ncol=2, strip.position = "left")+theme(strip.placement = "outside")+scale_y_continuous(breaks=pretty_breaks(n=3))
pdf("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/figures/sin_l23.pdf", height = 2.8, width = 2.2)
plot(p)
dev.off()
int_res=readRDS("interaction_results_readjusted_y.rds") %>% dplyr::filter(fdr<0.05) %>% mutate(dif=app_amp-wt_amp) %>% arrange(dif)
 res_y_aw <- readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/young/res_y_aw.rds") %>% dplyr::filter(padj<.1, genotype=="WT")

 res_y_aw=left_join(res_y_aw %>% dplyr::select(cluster, gene), int_res)%>% na.omit

  p=make_plot(c("Efna2", "Bcar1", "Pik3r1", "Smim20", "Kdm7a", "Tmem158", "P4ha1", 
               "Hcfc2"), "DGsg", center = F) +facet_wrap(gene~., scales = "free", ncol=2, strip.position = "left")+theme(strip.placement = "outside")+scale_y_continuous(breaks=pretty_breaks(n=3))
 pdf("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/figures/sin_dg.pdf", height = 2.8, width = 2.2)
 plot(p)
 dev.off()
 
 
 p=make_plot(sym2entrez$SYMBOL, "L5b", center = F) +facet_wrap(gene~., scales = "free", ncol=2, strip.position = "left")+theme(strip.placement = "outside")+scale_y_continuous(breaks=pretty_breaks(n=3))
 pdf("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/figures/sin_dg.pdf", height = 2.8, width = 2.2)
 plot(p)
 dev.off()
 
 
 
 p=make_plot(c("Bcar1", "Efna3","Efna2"  ), "DGsg", center = F) +facet_wrap(gene~., scales = "free", ncol=3, strip.position = "top")+theme(strip.placement = "outside")+scale_y_continuous(breaks=pretty_breaks(n=3))
 pdf("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/figures/sin_dg_rap1_fdr0.1.pdf", height = 1, width =3.1)
 plot(p)
 dev.off()
 
 
 p=make_plot("Egr2", "L4", center = F) +facet_wrap(gene~., scales = "free", ncol=3, strip.position = "top")+theme(strip.placement = "outside")+scale_y_continuous(breaks=pretty_breaks(n=3))
 #pdf("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/figures/sin_dg_rap1_fdr0.1.pdf", height = 1, width =3.1)
 plot(p)
 dev.off()