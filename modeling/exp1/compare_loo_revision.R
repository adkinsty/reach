# Compares models using leave-one-out (LOO) expected log posterior density (elpd)
# Author: Tyler Adkins
# Date: 2/12/21

library(tidyverse)
library(ggthemes)
library(loo)

wong <- c("#000000", "#e69f00", "#56b4e9", "#009e73", "#f0e442", "#0072b2", "#d55e00", "#cc79a7")

ggplot(data = tibble(x = c(5, 4), y = c(3, 2)), aes(x = x, y = y)) + geom_point(colour = wong[3])

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

models <- c("MEG", "MELV", "MEU", "MSWU", "UH", "H")

meg_loo <- loo(readRDS("/Users/adkinsty/files/offline_reach/exp1/rds/MEG.rds") ,cores = 1)
mleg_loo <- loo(readRDS("/Users/adkinsty/files/offline_reach/exp1/rds/MLEG.rds"),cores = 1)
meu_loo <- loo(readRDS("/Users/adkinsty/files/offline_reach/exp1/rds/MEU.rds") ,cores = 1)
mswu_loo <- loo(readRDS("/Users/adkinsty/files/offline_reach/exp1/rds/MSWU.rds") ,cores = 1)
h_loo <- loo(readRDS("/Users/adkinsty/files/offline_reach/exp1/rds/H.rds"),cores = 1)
uh_loo <- loo(readRDS("/Users/adkinsty/files/offline_reach/exp1/rds/UH.rds"),cores = 1)
#ulh_loo <- loo(readRDS("/Users/adkinsty/files/offline_reach/exp1/rds/ULH.rds"),cores = 4)
#meg_vp_loo <- loo(readRDS("/Users/adkinsty/files/offline_reach/exp1/rds/MEG_vp.rds") ,cores = 4)

comparison <- loo_compare(list(meg = meg_loo,
      mleg = mleg_loo,
      meu = meu_loo,
      mswu = mswu_loo,
      uh = uh_loo,
      h = h_loo))

loo_compare(list(meg = meg_loo,
      mleg = mleg_loo,
      meu = meu_loo,
      mswu = mswu_loo,
      h = h_loo))

loo_compare(list(meg = meg_loo,
      mleg = mleg_loo,
      h = h_loo))


models_plt <- c("MLEG", "MEG", "MSWU", "UH")

pdat <- tibble(ELPD = c(mleg_loo$elpd_loo,
            meg_loo$elpd_loo,
            mswu_loo$elpd_loo,
            uh_loo$elpd_loo),
       se = c(mleg_loo$se_elpd,
            meg_loo$se_elpd,
            mswu_loo$se_elpd,
            uh_loo$se_elpd),
       M = models_plt) %>%
  mutate(M = factor(M, levels=models_plt))

labels <- c("LA", "MEG", "NLVP", "UH")

pdat %>%
  ggplot(aes(x=M,y=ELPD,ymax=ELPD+se,ymin=ELPD-se,colour=M)) +
  geom_point() +
  geom_errorbar(width=0) +
  scale_x_discrete("Model",labels=labels) +
  scale_y_continuous("LOO ELPD") +
  scale_colour_manual("Model", values=wong[c(2,4,3,1)], labels = labels) + 
  theme_tufte(base_size=10,base_family="sans") +
  theme(axis.line=element_line(size=.25),
        legend.pos="none")
ggsave("visuals/raw/exp1/revision/exp1_LOOIC.pdf",units="in",height=2,width=2)
