# Title     : Multi-normal regression
# Objective : Interface with stan model
# Created by: adkinsty
# Created on: 6/20/20

library(rstan)
library(tidyverse)
library(shinystan)
library(ggthemes)


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
 
setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2) %>%  # reach filter
  filter(rt < 1000) %>%  # speed filter
  filter(id == 524)  # subject filter

x <- data$x
y <- data$y
obs <- as.array(cbind(x,y))
colnames(obs) <- NULL

ratio <- data$ratio

input <- list(N=length(x),
              J=ncol(obs),
              ratio=ratio,
              obs=obs,
              get_obs_rep=1,
              get_log_lik=0,
              prior_only=0)

fit <- stan(
    file="modeling/stan/multi_norm_lin.stan",
    data=input,
    chains=2,
    iter=2000,
    cores=1)

x_rep = rstan::extract(fit, 'obs_rep')$obs_rep[,,1]
y_rep = rstan::extract(fit, 'obs_rep')$obs_rep[,,2]

launch_shinystan(fit)
# 
# yrep = rstan::extract(fit, 'y_rep')$y_rep[idx_rep,]






