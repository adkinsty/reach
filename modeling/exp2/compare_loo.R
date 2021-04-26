# Compares models using leave-one-out (LOO) expected log posterior density (elpd)

library(tidyverse)
library(ggthemes)
library(loo)


pal <- c("#000000", # Black
         "#cc79a7", # Reddish purple
         "#f0e442", # Yellow
         "#009e73", # Bluish green
         "#e69f00", # Orange
         "#d55e00", # Vermillion
         "#56b4e9", # Sky blue
         "#0072b2") # Blue

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

h_loo <- loo(readRDS("modeling/exp2/stan/rds/H.rds"),cores = 4)
meg_loo <- loo(readRDS("modeling/exp2/stan/rds/MEG.rds") ,cores = 4)
mleg_loo <- loo(readRDS("modeling/exp2/stan/rds/MLEG.rds"),cores = 4)
# mswu_loo <- loo(readRDS("modeling/exp2/stan/rds/MSWU.rds") ,cores = 4)
# mswu_g_loo <- loo(readRDS("modeling/exp2/stan/rds/MSWU_gamma.rds") ,cores = 4)
meu_loo <- loo(readRDS("modeling/exp2/stan/rds/MEU.rds") ,cores = 4)
meg_vp_loo <- loo(readRDS("modeling/exp2/stan/rds/MEG_vp.rds") ,cores = 4)


print(h_loo)
print(meg_loo)
print(mleg_loo)
# print(mswu_loo)
# print(mswu_g_loo)
print(meu_loo)
print(loo_compare(h_loo,meg_loo,mleg_loo))
print(loo_compare(list(meg_vp=meg_vp_loo,meg=meg_loo)))
# print(loo_compare(x = list(h=h_loo,meg=meg_loo,mleg=mleg_loo,mswu=mswu_loo,mswu_g=mswu_g_loo,meu=meu_loo)))
print(loo_compare(list(h=h_loo,meg=meg_loo,mleg=mleg_loo,meu=meu_loo)))

tibble(ELPD = c(meu_loo$elpd_loo,h_loo$elpd_loo,meg_loo$elpd_loo,mleg_loo$elpd_loo),
       se = c(meu_loo$se_elpd_loo,h_loo$se_elpd_loo,meg_loo$se_elpd_loo,mleg_loo$se_elpd_loo),
       M = c("MEU","H","MEV","MELV")) %>%
  mutate(M = factor(M, levels=c("MEU","H","MEV","MELV"))) %>%
  ggplot(aes(x=M,y=ELPD,ymax=ELPD+se,ymin=ELPD-se,colour=M)) +
  geom_point() +
  geom_errorbar(width=0) +
  scale_x_discrete("Model") +
  scale_y_continuous("LOO") +
  scale_colour_manual("Model",
                      values=c("#0072b2","#e69f00","#009e73","#cc79a7"),
                      labels=c("MEU", "H","MEV","VP")) +  
  theme_tufte(base_size=10,base_family="sans") +
  theme(axis.line=element_line(size=.25),
        legend.pos="none")
ggsave("visuals/raw/exp2/exp2_LOOIC.pdf",units="in",height=2,width=2)
