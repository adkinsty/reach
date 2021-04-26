# Script for comparing the models in terms of LOO elpd
# Author: Tyler Adkins
# Date: 9/10/20

library(tidyverse)
library(ggthemes)
library(loo)

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

# h_loo <- loo(read_rds("modeling/exp1/stan/rds/H.rds") ,cores = 4)
# write_rds(h_loo,"modeling/exp1/loo/h_loo.rds")
h_loo <- read_rds("modeling/exp1/loo/h_loo.rds")

# meg_loo <- loo(read_rds("modeling/exp1/stan/rds/MEG.rds") ,cores = 4)
# write_rds(meg_loo,"modeling/exp1/loo/meg_loo.rds")
meg_loo <- read_rds("modeling/exp1/loo/meg_loo.rds")

# meg_vp_loo <- loo(read_rds("modeling/exp1/stan/rds/MEG_vp.rds") ,cores = 4)
# write_rds(meg_vp_loo,"modeling/exp1/loo/meg_vp_loo.rds")
meg_vp_loo <- read_rds("modeling/exp1/loo/meg_vp_loo.rds")

# mleg_loo <- loo(read_rds("modeling/exp1/stan/rds/MLEG.rds"),cores = 4)
# write_rds(mleg_loo,"modeling/exp1/loo/mleg_loo.rds")
mleg_loo <- read_rds("modeling/exp1/loo/mleg_loo.rds")

# meu_loo <- loo(read_rds("modeling/exp1/stan/rds/MEU.rds") ,cores = 4)
# write_rds(meu_loo,"modeling/exp1/loo/meu_loo.rds")
meu_loo <- read_rds("modeling/exp1/loo/meu_loo.rds")

# mswu_g_loo <- loo(read_rds("modeling/exp1/stan/rds/MSWU_gamma.rds") ,cores = 4)
# meu_g_loo <- loo(read_rds("modeling/exp1/stan/rds/MEU.rds") ,cores = 4)

print(h_loo)
print(meg_loo)
print(mleg_loo)
# print(mswu_loo)
# print(mswu_g_loo)
print(meu_loo)
print(loo_compare(h_loo,meg_loo,mleg_loo))
print(loo_compare(list(meu=meu_loo,h=h_loo,meg=meg_loo,melg=mleg_loo)))
print(loo_compare(list(meg=meg_loo,meg_vp=meg_vp_loo)))


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
                      labels=c("MEU", "H","MEV","MELV")) +  
  theme_tufte(base_size=10,base_family="sans") +
  theme(axis.line=element_line(size=.25),
        legend.pos="none")
ggsave("visuals/raw/exp1/exp1_LOOIC.pdf",units="in",height=2,width=2)
