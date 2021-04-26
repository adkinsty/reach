# Title     : Univariate Max Expected Value Model
# Objective : Interface with stan model
# Created by: adkinsty
# Created on: 6/20/20

library(rstan)
library(tidyverse)
library(shinystan)
library(ggthemes)
library(shotGroups)

options(mc.cores = parallel::detectCores())


setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp2/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2) %>%  # reach filter
  filter(rt < 1000) %>%  # speed filter
  filter(id == 524)  # subject filter

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2) %>%  # reach filter
  filter(rt < 1000) %>%  # speed filter
  filter(id == 524)  # subject filter

grid <- seq(-5,5,.0312) # approx -48 to 48 by 1 px

sigma <- array(c(var(train$x),cov(train$x,train$y),
                 cov(train$x,train$y),var(train$y)),dim=c(2,2))

get_p_gain <- function(x,s=sigma) {
  aim <- c(x, 0)
  p <- pmvnEll(r = 1, x0 = c(0,0), e = diag(2), mu = aim, sigma = s)
  return(p)
}
get_p_loss <- function(x,s=sigma) {
  aim <- c(x, 0)
  p <- pmvnEll(r = 1, x0 = c(-1,0), e = diag(2), mu = aim, sigma = s)
  return(p)
}

p_gain <- sapply(X = grid, FUN = get_p_gain, simplify = T, USE.NAMES = F)
p_loss <- sapply(X = grid, FUN = get_p_loss, simplify = T, USE.NAMES = F)

x <- data$x
y <- data$y 

obs <- data$x
gain <- data$gain
loss <- data$loss

input <- list(N=length(x),
              G=length(grid),
              C=length(unique(loss)),
              cond=sort(unique(loss)),
              prior_sigma=var(train$x),
              grid=grid,
              p_gain=p_gain,
              p_loss=p_loss,
              gain=gain,
              loss=loss,
              obs=obs,
              get_obs_rep=1,
              get_log_lik=1,
              prior_only=0)

fit <- stan(
  file="modeling/stan/max_EV_uni_normal.stan",
  data=input,
  chains=2,
  iter=2000,
  cores=1)

obs_rep = rstan::extract(fit, 'obs_rep')$obs_rep

launch_shinystan(fit)

