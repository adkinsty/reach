# Title     : Linear regression model
# Objective : Interface with stan model
# Created by: adkinsty
# Created on: 6/22/20

library(rstan)
library(tidyverse)
library(shinystan)
library(ggthemes)
library(shotGroups)
library(bayesplot)
library(tidyr)

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2) %>%
  mutate(z_cond = ifelse(ratio==0,-1,ifelse(ratio==-1,0,1)))

la_df <- read_tsv("data/raw/loss_aversion/All_lambdas.csv") %>%
  filter(subNum %in% data$id) %>% arrange(subNum) %>% 
  mutate(id = 1:n(),
         z_la = scale(lambdas)[,1])

input <- list(
  N=nrow(data),                                      # num trials
  K=2,                                               # num trial predictors
  J=length(unique(data$sub_id)),                     # num subjects
  L=1,                                               # num subject predictors
  jj=data$sub_id,                                    # subject for trial
  x=as.matrix(cbind(rep(1,nrow(data)),data$z_cond)),  # trial predictors
  u=as.matrix(cbind(rep(1,nrow(la_df)))), # subject predictors
  y=data$x,                                          # observations
  get_yrep=1,  
  get_log_lik=1,
  prior_only=0
)

model <- stan_model(
  file="modeling/stan/LRc.stan",
  model_name="LRc"
)

mcmc_est <- sampling(
  object=model,
  data=input,
  chains=4,
  iter=1000,
  warmup=900,
  cores=4,
  control=list(max_treedepth=10))
write_rds(mcmc_est,path="modeling/stan/rds/LRc.rds")