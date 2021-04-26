# Title     : Hierarchical Univariate Normal Maximum Expected Utility Model (non-centered param)
# Objective : Interface with stan model
# Created by: adkinsty
# Created on: 6/22/20

library(rstan)
library(tidyverse)
library(shinystan)
library(ggthemes)
library(shotGroups)


options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp2/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2) %>%
  filter(sub_id < 10)

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2) %>%
  mutate(loss = ifelse(loss==0,-.01,loss)) %>%
  filter(sub_id < 10)

sub_id <- data$sub_id
S <- length(unique(sub_id))
grid <- seq(0,1,.001) # approx -48 to 48 by 1 px
G <- length(grid)
p_gain <- array(rep(0,G*S),dim=c(G,S))
p_loss <- array(rep(0,G*S),dim=c(G,S))
for (s in 1:S) {
  tmp <- train %>% filter(sub_id == s & block > 5)
  sigma <- array(c(var(tmp$x),cov(tmp$x,tmp$y),cov(tmp$x,tmp$y),var(tmp$y)),dim=c(2,2))
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
  p_gain[,s] <- sapply(X = grid, FUN = get_p_gain, simplify = T, USE.NAMES = F)
  p_loss[,s] <- sapply(X = grid, FUN = get_p_loss, simplify = T, USE.NAMES = F)
}

x <- data$x
y <- data$y
obs <- x
colnames(obs) <- NULL

gain <- data$gain
loss <- data$loss
cond_id <- data$cond_id

input <- list(N=length(obs),
              G=length(grid),
              C=length(unique(loss)),
              S=length(unique(sub_id)),
              K=2,
              cond=sort(unique(loss)),
              sub=1:S,
              sub_id=sub_id,
              cond_id=cond_id,
              grid=grid,
              p_gain=p_gain,
              p_loss=p_loss,
              gain=gain,
              loss=loss,
              obs=obs,
              get_obs_rep=1,
              get_log_lik=1,
              prior_only=0)

model <- stan_model(file="modeling/stan/scratch_model.stan",
                    model_name="scratch_model",
                    allow_undefined = TRUE,
                    includes = paste0('\n#include "',
                                      file.path(getwd(),
                                                'modeling/stan/hpp/get_max_idx.hpp'), 
                                      '"\n'))
fit <- sampling(
  object=model,
  data=input,
  chains=3,
  iter=3000,
  cores=3)

obs_rep <- rstan::extract(fit, 'obs_rep')$obs_rep
launch_shinystan(fit)

