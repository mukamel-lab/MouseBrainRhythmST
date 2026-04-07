library(tidyverse)


sexdrg <- readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/sex_int_geno_spec/interaction_results_readjusted.rds")
chr_ano <- readRDS("~/desp1/precast/prec_c25q25g3000/figures/sex_rhythmicity/chr_ano.rds")

sexdrg=split(sexdrg, sexdrg$cluster)


res = lapply(sexdrg, function(gns) {
  
  sigs=gns %>% dplyr::filter(fdr<0.1) %>% pull(gene)
  if(length(sigs>=30)){
  
  univ_df <- data.frame(
    gene = gns$gene,
    chr  = plyr::mapvalues(gns$gene, from = chr_ano$gene, to = chr_ano$chr),
    stringsAsFactors = FALSE
  )
  
  hit_df <- data.frame(
    gene = sigs,
    chr  = plyr::mapvalues(sigs, from = chr_ano$gene, to = chr_ano$chr),
    stringsAsFactors = FALSE
  )
  
  univ_counts <- table(univ_df$chr)
  hit_counts  <- table(hit_df$chr)
  all_chr     <- names(univ_counts)
  
  hg <- data.frame(
    chr = all_chr,
    M   = as.integer(univ_counts[all_chr]),                # # genes on chr in universe
    k   = as.integer(hit_counts[all_chr]),                  # # hits on chr
    stringsAsFactors = FALSE
  )
  hg$k[is.na(hg$k)] <- 0
  N <- sum(hg$M)  
  K <- sum(hg$k)    
  
  hg$pval <- mapply(function(k, M) {
    phyper(q = k - 1,  
           m = M,
           n = N - M,
           k = K,
           lower.tail = FALSE)
  }, hg$k, hg$M)
  
  hg$padj <- p.adjust(hg$pval, method = "BH")

  obs <- hg$k

    exp_props <- hg$M / N
  
  chi <- chisq.test(x           = obs,
                    p           = exp_props,
                    simulate.p.value = TRUE,
                    B           = 10000)
  
  list(
    hyper = hg,
    chi   = chi
  )} else { return(NULL)}
})





