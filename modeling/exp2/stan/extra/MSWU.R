# Title     : Experiment 2, maximum subjective weighted utility model
# Objective : Interface with stan model
# Created by: adkinsty
# Updated on: 6/29/20

library(rstan)
library(tidyverse)
library(shinystan)
library(ggthemes)
library(shotGroups)
library(bayesplot)
library(optimx)
library(mvtnorm)
library(raster)

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp2/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

# key variables
jj <- data$sub_id
cc <- data$cond_id
J <- length(unique(jj))
N <- nrow(data)
C <- length(unique(cc))
cond <- sort(unique(data$loss))
y <- data$x

mu_cond <- read_rds("data/clean/exp2/MSWU_aims.rds")

input <- list(N=N,
              J=J,
              C=C,
              cc=cc,
              jj=jj,
              aim=mu_cond,
              y=y,
              get_gq=1,
              prior_only=0)

model <- stan_model(file="modeling/exp2/stan/aim.stan",model_name="MSWU")

mcmc_est <- sampling(
  object=model,
  data=input,
  chains=4,
  iter=5000,
  warmup=2500,
  cores=4)
write_rds(mcmc_est,path="modeling/exp2/stan/rds/MSWU.rds")