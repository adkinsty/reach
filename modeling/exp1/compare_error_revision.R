# Script for comparing the models in terms of posterior prediction error
# Author: Tyler Adkins
# Date: 9/10/20

library(tidyverse)
library(ggthemes)

wong <- c("#000000", "#e69f00", "#56b4e9", "#009e73", "#f0e442", "#0072b2", "#d55e00", "#cc79a7")

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

data <- read_csv("data/clean/exp1/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2) %>%
  mutate(lr = ifelse(ratio==0,0,ifelse(ratio==-1,1,2)))

y <- data$x # observation
ymm <- y * 8.55 # mm per radius

data %>% 
  mutate(ratio=factor(loss,levels=c(0,-1,-3,-5,-15)),
         x = x*8.55,
         dist=ifelse(dist==1,"Near","Far")) %>%
  ggplot(aes(x=x,y=ratio)) +
  facet_wrap(dist~.,nrow=2) +
  geom_point(size=.1,alpha=.1,position = position_jitter(width = 0, height = .25)) +
  stat_summary() + 
  scale_x_continuous("Observed endpoint (mm)",breaks=seq(-10,10,2)) +
  scale_y_discrete("Loss ratio",labels=c("0:1","1:1","3:3","5:1","15:3")) +
  geom_vline(xintercept = 0, colour="black",linetype="dashed") +
  coord_cartesian(xlim=c(0,8.55)) +
  theme_tufte(base_size=10,base_family="sans") +
  theme(axis.line=element_line(size=.25),
        legend.pos="none")
ggsave("visuals/raw/exp1/revision/exp1_obs.pdf",units="in",height=2,width=3.5)


models <- c("MEG", "MLEG", "MEU", "MSWU", "UH", "H")

# posterior replications
meg_yrep <- rstan::extract(readRDS("/Users/adkinsty/files/offline_reach/exp1/rds/MEG.rds"), 'yrep')$yrep
mleg_yrep <- rstan::extract(readRDS("/Users/adkinsty/files/offline_reach/exp1/rds/MLEG.rds"), 'yrep')$yrep
meu_yrep <- rstan::extract(readRDS("/Users/adkinsty/files/offline_reach/exp1/rds/MEU.rds"), 'yrep')$yrep
mswu_yrep <- rstan::extract(readRDS("/Users/adkinsty/files/offline_reach/exp1/rds/MSWU.rds"), 'yrep')$yrep
uh_yrep <- rstan::extract(readRDS("/Users/adkinsty/files/offline_reach/exp1/rds/UH.rds"), 'yrep')$yrep
h_yrep <- rstan::extract(readRDS("/Users/adkinsty/files/offline_reach/exp1/rds/H.rds"), 'yrep')$yrep
#ulh_yrep <- rstan::extract(readRDS("modeling/exp1/stan/rds/ULH.rds"), 'yrep')$yrep
#meg_vp_yrep <- rstan::extract(readRDS("modeling/exp1/stan/rds/MEG_vp.rds"), 'yrep')$yrep

mu_meg_yrep <- colMeans(meg_yrep)
mu_mleg_yrep <- colMeans(mleg_yrep)
mu_meu_yrep <- colMeans(meu_yrep)
mu_mswu_yrep <- colMeans(mswu_yrep)
mu_uh_yrep <- colMeans(uh_yrep)
mu_h_yrep <- colMeans(h_yrep)
#mu_ulh_yrep <- colMeans(ulh_yrep)
#mu_meg_vp_yrep <- colMeans(meg_vp_yrep)

yrep_data <- data %>%
  mutate(MEG=mu_meg_yrep,
         MLEG=mu_mleg_yrep,
         MEU=mu_meu_yrep,
         MSWU=mu_mswu_yrep,
         UH=mu_uh_yrep,
         H = mu_h_yrep) %>%
  pivot_longer(cols=c(MEG,MLEG,MSWU,UH),names_to="model",values_to="yrep")

models_plt <- c("MLEG", "MEG", "MSWU", "UH")
labels <- c("LA", "MEG", "NLVP", "UH")

yrep_data %>%
  mutate(ratio=factor(loss,levels=c(0,-1,-3,-5,-15)),
         obs=x*8.55,
         yrep=yrep*8.55,
         dist=ifelse(dist==1,"8.55 mm spread","11.75 mm spread"),
         model = factor(model,levels = models_plt))%>%
  ggplot(aes(x = ratio, y = yrep - obs, fill = model, group = interaction(model, ratio, dist))) + 
  facet_wrap(dist~.,nrow=2) +
  stat_summary(geom="col",position=position_dodge(.75),width=.75) +
  stat_summary(geom="errorbar",fun.data = "mean_cl_boot", position=position_dodge(.75),width=.1, colour="grey")+
  geom_hline(yintercept = 0,linetype="dashed",size=.25) +
  scale_y_continuous("Predictive error (mm)",breaks=seq(-3,2,1)) +
  scale_x_discrete("Loss:Gain",labels=c("0:1","1:1","3:3","5:1","15:3")) +
  scale_fill_manual("Model", values=wong[c(2,4,3,1)], labels=labels) + 
  coord_cartesian(ylim=c(-3.5,1.5)) +
  theme_tufte(base_size=10,base_family="sans") +
  theme(axis.line=element_line(size=.25),
        legend.key.size = unit(.15,"in"),
        legend.position = "top")
ggsave("visuals/raw/exp1/revision/exp1_ppc_mean_error.pdf",
       units="in",height=3,width=3)

