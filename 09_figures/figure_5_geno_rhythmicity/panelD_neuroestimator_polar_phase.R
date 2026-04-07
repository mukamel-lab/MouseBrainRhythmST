library(tidyverse)
library(pbmcapply)
library(lmtest)

setwd("/home/agelber/desp1/precast/prec_c25q25g3000/neuroestimator")


meta.data=readRDS("../objects/md_ser_neuroE.rds")

meta.data$t_c=as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s=as.numeric(gsub("ZT", "", meta.data$time))
meta.data= meta.data%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                               t_c=round(cos(2*pi*(t_c)/24)))

md =meta.data %>% group_by(sample, time, t_s,
                           t_c,
                           sex, age, genotype, anot) %>% 
  summarise(neuroE=mean(predicted_activity))

md=split(md, md$anot)

md=lapply(md, function(df) { split(df, paste0(df$genotype,"_", df$age))})
res=pbmclapply(md, function(clst){
  pbmclapply(clst, function(df) {
f1=glm(neuroE ~1+t_s+t_c,data=df)
f2=glm(neuroE ~1,data=df)

full_ll = logLik(f1)
red_ll = logLik(f2)

lr = -2 * (as.numeric(red_ll)-as.numeric(full_ll))

p.val = pchisq(lr, df = 2, lower.tail = FALSE)

data.frame(
  cluster = df$anot[1],
  age = df$age[1],
  genotype = df$genotype[1],
  pval = p.val,
  t_s = coefficients(f1)[2],
  t_c = coefficients(f1)[3],
  intcp = coefficients(f1)[1]
) %>%
  mutate(
    phi = atan2(t_s, t_c) %% (2 * pi),
    amp = sqrt(t_c ^ 2 + t_s ^ 2),
    rel_amp = amp / intcp
  )

  })  %>% bind_rows()}) %>% bind_rows(.id="cluster") %>%
  {df=.;rownames(df)=NULL;df}


res$padj=p.adjust(res$pval, method = "BH")
