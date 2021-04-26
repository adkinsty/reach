# Script for comparing the models in terms of posterior prediction error
# Author: Tyler Adkins
# Date: 9/10/20

library(tidyverse)
library(ggthemes)

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
  geom_point(size=.1,alpha=.5,position = position_jitter(),colour="black") +
  stat_summary() + 
  scale_x_continuous("Observed endpoint (mm)",breaks=seq(-10,10,2)) +
  scale_y_discrete("Loss ratio",labels=c("0:1","1:1","3:3","5:1","15:3")) +
  geom_vline(xintercept = 0, colour="black",linetype="dashed") +
  coord_cartesian(xlim=c(0,8.55)) +
  theme_tufte(base_size=10,base_family="sans") +
  theme(axis.line=element_line(size=.25),
        legend.pos="none")
ggsave("visuals/raw/exp1/exp1_obs.pdf",units="in",height=2.25,width=4)

# posterior draws

h_yrep <- rstan::extract(readRDS("modeling/exp1/stan/rds/H.rds"), 'yrep')$yrep
meg_yrep <- rstan::extract(readRDS("modeling/exp1/stan/rds/MEG.rds"), 'yrep')$yrep
mleg_yrep <- rstan::extract(readRDS("modeling/exp1/stan/rds/MLEG.rds"), 'yrep')$yrep
meu_yrep <- rstan::extract(readRDS("modeling/exp1/stan/rds/MEU.rds"), 'yrep')$yrep
# mswu_yrep <- rstan::extract(readRDS("modeling/exp1/stan/rds/MSWU.rds"), 'yrep')$yrep
# mswu_g_yrep <- rstan::extract(readRDS("modeling/exp1/stan/rds/MSWU_gamma.rds"), 'yrep')$yrep

mu_h_yrep <- colMeans(h_yrep)
mu_meg_yrep <- colMeans(meg_yrep)
mu_mleg_yrep <- colMeans(mleg_yrep)
mu_meu_yrep <- colMeans(meu_yrep)

# mu_mswu_yrep <- colMeans(mswu_yrep)
# mu_mswu_g_yrep <- colMeans(mswu_g_yrep)

yrep_data <- data %>%
  mutate(H = mu_h_yrep,MEG=mu_meg_yrep,
         MLEG=mu_mleg_yrep,MEU=mu_meu_yrep) %>%
  pivot_longer(cols=c(MEU,H,MEG,MLEG),names_to="model",values_to="yrep") %>%
  mutate(model = factor(model,c("MEU","H","MEG","MLEG")))

yrep_data %>%
  mutate(ratio=factor(loss,levels=c(0,-1,-3,-5,-15)),
         x=x*8.55,
         yrep=yrep*8.55,
         dist=ifelse(dist==1,"8.55 mm spread","11.75 mm spread"))%>%
  ggplot(aes(x=ratio,y=x-yrep,fill=model)) + 
  facet_wrap(dist~.,nrow=2) +
  stat_summary(geom="col",position=position_dodge(.75),width=.75) +
  stat_summary(geom="errorbar",position=position_dodge(.75),width=.1)+
  geom_hline(yintercept = 0,linetype="dashed",size=.25) +
  scale_y_continuous("Predictive error (mm)",breaks=seq(-2,3,1)) +
  scale_x_discrete("Loss ratio",labels=c("0:1","1:1","3:3","5:1","15:3")) +
  scale_fill_manual("",values=c("#0072b2","#e69f00","#009e73","#cc79a7"),
                      labels=c("MEU", "H","MEV","MELV")) + 
  coord_cartesian(ylim=c(-1.5,3)) +
  theme_tufte(base_size=10,base_family="sans") +
  theme(axis.line=element_line(size=.25),
        legend.key.size = unit(.15,"in"),
        legend.position = "top")
ggsave("visuals/raw/exp1/exp1_ppc_mean_error.pdf",
       units="in",height=3,width=3.25)
