library(tidyverse)
library(AnnotationDbi)
library(org.Mm.eg.db)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new")

y=readRDS("young/young_gsea_rhth_amp_app_vs_wt_kegg.rds") 
# %>% rowwise()%>%mutate(
#   genes=paste0(as.character(mapIds(
#     org.Mm.eg.db, str_split_1(core_enrichment, "/"),  'SYMBOL','ENTREZID'
#   )),collapse=",")
# )

cnts_y=y %>% group_by(cluster) %>% summarise(n2=sum(p.adjust<0.2), 
                                           n1=sum(p.adjust<0.1),
                                           n5=sum(p.adjust<0.05))
o=readRDS("old/old_gsea_rhth_amp_app_vs_wt_kegg.rds") 
# %>% rowwise()%>%mutate(
#   genes=paste0(as.character(mapIds(
#     org.Mm.eg.db, str_split_1(core_enrichment, "/"),  'SYMBOL','ENTREZID',
#   )),collapse=",")
# )

cnts_o=o %>% group_by(cluster) %>% summarise(n2=sum(p.adjust<0.2), 
                                             n1=sum(p.adjust<0.1),
                                             n5=sum(p.adjust<0.05))

df_plot <- y %>%
  filter(p.adjust < 0.2, cluster=="DGsg") %>%
  mutate(color_group = case_when( 
  
    p.adjust < 0.1 & NES >  0 ~ "padj<0.1, NES>0",
     p.adjust < 0.2 & NES >  0 ~ "padj<0.2, NES>0",
    
    p.adjust < 0.1 & NES <  0 ~ "padj<0.1, NES<0",
    p.adjust < 0.2 & NES <  0 ~ "padj<0.2, NES<0"
  ))

# Plot
ggplot(df_plot, aes(
  x = reorder(Description, NES),
  y = NES,
  fill = color_group
)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    name = "Category",
    values = c(
      "padj<0.1, NES>0"  = "#4783B5",
      "padj<0.2, NES>0"  = "#B0C4DE",
      "padj<0.2, NES<0"  = "#FFCF77",
      "padj<0.1, NES<0"  = "#F78C1E"
    )
  ) +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "Normalized Enrichment Score (NES)",
    title = "Significant KEGG Terms by NES",
    subtitle = "Colored by p.adjust threshold and NES direction"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )


df_plot <- o %>%
  filter(p.adjust < 0.2, cluster=="L5b") %>%
  mutate(color_group = case_when( 
    
    p.adjust < 0.1 & NES >  0 ~ "padj<0.1, NES>0",
    p.adjust < 0.2 & NES >  0 ~ "padj<0.2, NES>0",
    
    p.adjust < 0.1 & NES <  0 ~ "padj<0.1, NES<0",
    p.adjust < 0.2 & NES <  0 ~ "padj<0.2, NES<0"
  ))

# Plot
ggplot(df_plot, aes(
  x = reorder(Description, NES),
  y = NES,
  fill = color_group
)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    name = "Category",
    values = c(
      "padj<0.1, NES>0"  = "#4783B5",
      "padj<0.2, NES>0"  = "#B0C4DE",
      "padj<0.2, NES<0"  = "#FFCF77",
      "padj<0.1, NES<0"  = "#F78C1E"
    )
  ) +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "Normalized Enrichment Score (NES)",
    title = "Significant KEGG Terms by NES",
    subtitle = "Colored by p.adjust threshold and NES direction"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

