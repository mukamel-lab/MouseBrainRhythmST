library(DESeq2)
library(tidyverse)
library(patchwork)
library(ggrepel)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_geno_sex_int")
deseqs_old=readRDS("~/desp1/precast/prec_c25q25g3000/deseq_geno_sex_int/sex_geno_old_wald_full/deseqs.rds")
sex_genes=c("Xist", "Tsix","Eif2s3y", "Ddx3y", 
            "Uty","Kdm5d", "Ddx3x", "Eif2sx", "Eif2s3x")

res=lapply(deseqs_old, function(x) {
  results(x, 
          #contrast =list( c("genotype_APP23_vs_WT","genotypeAPP23.sexF")),
          #name="genotypeAPP23.sexF",
          contrast =list( c("sex_F_vs_M","genotypeAPP23.sexF")),
          tidy = T)
}) %>% 
  bind_rows(.id="cluster") %>%
  dplyr::rename(gene=row) %>%
  filter(!(gene %in% sex_genes))

cluster_color=readRDS("../objects/cluster_color.rds")

res=res %>% group_by(cluster) %>% 
  arrange(padj) %>%
  mutate(top_sig = ifelse(row_number() <= 5, cluster, "non-sig")) %>%
  mutate(top_sig=ifelse(top_sig %in% names(cluster_color)[1:11],
                        top_sig, "non-sig")) %>%
  mutate(sig=ifelse(padj<0.05, cluster, "non-sig"))
           

ns=c("gray65")

names(ns)="non-sig"
cluster_color=c(cluster_color, ns)
res$cluster=gsub("Olfactory Tubercle", "Amygdala",res$cluster)

res=res %>% mutate(padj=ifelse(padj<1e-4, 1e-4, padj))

res=split(res, res$top_sig)

res$`non-sig`=sample_n(res$`non-sig` %>% ungroup(),size = 20000)

res=bind_rows(res, .id="top_sig")

ggplot(res %>% na.omit() %>% arrange(-padj))+
  geom_point(aes(x=log2(baseMean+1), y=log2FoldChange, col=sig, size=ifelse(sig=="non-sig", "ns", "s")))+
  geom_text_repel(data=res %>% filter(top_sig !="non-sig"),
             aes(x=log2(baseMean+1), y=log2FoldChange, label=gene),
             inherit.aes = F,
             max.overlaps = 10)+
  scale_color_manual(name="",values  =cluster_color)+
  scale_size_manual(values=c(.1,1))+theme_bw() +xlab("Normalized Expression") +
  ylab("Effect Size: Genotype + Genotype:Sex")+
  guides(size="none")

ggplot(res %>% na.omit() %>% arrange(-padj))+
  geom_hline(yintercept=4,linetype = "dashed")+
  geom_point(aes(x=log2FoldChange, y=-1*log10(padj), col=sig, size=ifelse(sig=="non-sig", "ns", "s")))+
  geom_text_repel(data=res %>% filter(top_sig !="non-sig"),
                   aes(x=log2FoldChange, y=-1*log10(padj), label=gene),
                   inherit.aes = F,
                   max.overlaps = 10)+
  scale_color_manual(name="",values  =cluster_color, breaks = head(names(cluster_color),-1))+
  scale_size_manual(values=c(.5,1))+theme_bw(base_size = 20) +ylab("-log10(FDR)") +
  xlab("Effect Size: Genotype + Genotype:Sex")+
  guides(size="none")+xlim(c(-1.7,1.7))+ylim(c(NA,4.3))
  

