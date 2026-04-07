library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(patchwork)

setwd("~/desp1/precast/prec_c25q25g3000/")

old=readRDS("deseq_rhth_int2/sex_int_geno_spec/main_effect/interaction_with_genotype/deseq_interaction_old.rds")
young=readRDS("deseq_rhth_int2/sex_int_geno_spec/main_effect/interaction_with_genotype/deseq_interaction_young.rds")

interact_old=pbmclapply(old, function(dds) {
  dds=dds %>% DESeq(fitType = "local", test="LRT", reduced = ~sex+genotype)
  results(dds, tidy=T) %>% dplyr::rename(gene=row)
  
  
}) 
interact_old=interact_old %>% bind_rows(.id="cluster")
interact_young=pbmclapply(young, function(dds) {
  
  dds=dds %>% DESeq(fitType = "glmGamPoi", test="LRT", reduced = ~sex+genotype)
  
  results(dds, tidy=T) %>% dplyr::rename(gene=row)
  
  
}) 

interact_young=interact_young %>% bind_rows(.id="cluster")

saveRDS(interact_old, "deseq_rhth_int2/sex_int_geno_spec/main_effect/interaction_with_genotype/int_res_old_LRT.rds")
saveRDS(interact_young, "deseq_rhth_int2/sex_int_geno_spec/main_effect/interaction_with_genotype/int_res_young_LRT.rds")

plot_deseq_box = function(dds, gene,
group_x   = "sex",
facet_by  = "genotype",
assay     = c("counts","vst","rlog"),
transform = TRUE) {
assay <- match.arg(assay)
# extract expression
if (assay %in% c("vst","rlog")) {
mat <- assay(dds)[gene, , drop=TRUE]
} else {
mat <- counts(dds, normalized = TRUE)[gene, , drop=TRUE]
if (transform) mat <- log2(mat + 1)
}
# build data.frame
df <- as.data.frame(colData(dds)[ , c(group_x, facet_by), drop=FALSE])
df$expression <- mat
df[[group_x]]  <- factor(df[[group_x]])
df[[facet_by]] <- factor(df[[facet_by]])
# plot
p <- ggplot(df, aes_string(x = group_x, y = "expression", fill = group_x)) +
geom_boxplot(outlier.shape = NA, width = 0.6) +
geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
facet_wrap(as.formula(paste("~", facet_by)), nrow = 1) +
theme_minimal(base_size = 11) +
xlab(group_x) +
ylab(paste0(ifelse(transform, "log2(", ""), assay,
ifelse(transform, "+1)", ""))) +
theme(
strip.background = element_blank(),
strip.text       = element_text(face = "bold"),
legend.position  = "none"
)
return(p)
}
