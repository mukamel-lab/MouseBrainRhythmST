
deseqs=lapply(agg_exp[c(1,11)], function(cl){
  # full model specified in design
  # reduced model specified in reduced
  try({
    suppressWarnings({
      counts1=cl
      
      md=meta.data %>% dplyr::rename(id=sample)
      maxs=apply(counts1 %>% dplyr::select(-gene) , 1, max)
      counts2=counts1[maxs>10, ]
      cs=colSums(counts2 %>% select(-gene))
      cs=names(cs)[cs>2e5]
      ngens=apply(counts2 %>% dplyr::select(-gene), 2, function(x){sum(x>0)})
      ngens=ngens/nrow(counts2)
      ngens=names(ngens)[ngens>0.9]
      
      md=md %>% filt(id %in% intersect(cs, ngens), genotype=="WT")
      
      counts2[,c("gene", md$id)] %>% 
        filt(gene %in% res_sig$gene) })})})


ctx=deseqs$`Cortex Layer 2/3`
colnames(ctx)[2:34]=paste0("ctx_",colnames(ctx)[2:34])

dgsg=deseqs$`Dentate Gyrus-sg`
colnames(dgsg)[2:33]=paste0("dg_",colnames(dgsg)[2:33])

exp_b=left_join(ctx, dgsg)

md2=data.frame(id=colnames(exp_b)[2:66])

md2=md2 %>% separate(id, c("reg", "sample"), "_", remove = F)

md2=left_join(md2, meta.data)


dsq=  DESeqDataSetFromMatrix(countData=exp_b, 
                             colData=md2, 
                             design=~sex+age+reg+t_s+t_c+reg:(t_s+t_c),
                             tidy = T) %>%
  estimateSizeFactors(type="poscounts") %>%
  estimateDispersions(fitType="glmGamPoi") %>%
  DESeq(test="LRT", reduced=~sex+age+reg+t_s+t_c)

r1=results(dsq, tidy = T)
