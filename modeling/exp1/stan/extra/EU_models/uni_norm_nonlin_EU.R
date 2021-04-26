# Title     : Univariate Normal Nonlinear Expected Utility Model
# Objective : Interface with stan model
# Created by: adkinsty
# Created on: 6/20/20

library(rstan)
library(tidyverse)
library(shinystan)
library(ggthemes)
library(shotGroups)

get_p_gain <- function(x,s=sigma) {
  # get probability of gain, given aim-point
  aim <- c(x, 0)
  p <- pmvnEll(r = 1, x0 = c(0,0), e = diag(2), mu = aim, sigma = s)
  return(p)
}
get_p_loss <- function(x,s=sigma) {
  # get probability of loss, given aim-point
  aim <- c(x, 0)
  p <- pmvnEll(r = 1, x0 = c(-1,0), e = diag(2), mu = aim, sigma = s)
  return(p)
}

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2) %>%  # reach filter
  filter(rt < 1000) %>%  # speed filter
  filter(id > 530) # subject filter

grid <- seq(0,1,.001) # .03125 = approx 1 px

x <- data$x
y <- data$y 
sigma <- array(c(var(x),cov(x,y),cov(x,y),var(y)),dim=c(2,2))

p_gain <- sapply(X = grid, FUN = get_p_gain, simplify = T, USE.NAMES = F)
p_loss <- sapply(X = grid, FUN = get_p_loss, simplify = T, USE.NAMES = F)

obs <- data$x
gain <- data$gain
loss <- data$loss
cond_id <- data$cond_id

input <- list(N=length(x),
              G=length(grid),
              C=length(unique(loss)),
              cond=sort(unique(loss)),
              grid=grid,
              p_gain=p_gain,
              p_loss=p_loss,
              cond_id=data$cond_id,
              gain=gain,
              loss=loss,
              obs=obs,
              get_obs_rep=1,
              get_log_lik=1,
              prior_only=0)

model <- stan_model(file="modeling/stan/uni_norm_nonlin_EU.stan",
  model_name="uni_norm_nonlin_EU",
  allow_undefined = TRUE,
  includes = paste0('\n#include "',file.path(getwd(),'modeling/stan/hpp/get_max_idx.hpp'), '"\n'))
fit <- sampling(
  object=model,
  data=input,
  chains=2,
  iter=2000,
  cores=2,
  control=list(max_treedepth=10))
obs_rep <- rstan::extract(fit, 'obs_rep')$obs_rep
launch_shinystan(fit)

