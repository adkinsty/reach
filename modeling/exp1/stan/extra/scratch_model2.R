# Title     : Hierarchical Univariate Normal Maximum Expected Utility Model (non-centered param)
# Objective : Interface with stan model
# Created by: adkinsty
# Created on: 6/22/20

library(rstan)
library(tidyverse)
library(shinystan)
library(ggthemes)
library(shotGroups)
library(bayesplot)


options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp2/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2) 

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2) 

sub_id <- data$sub_id
S <- length(unique(sub_id))
grid <- seq(0,1,.01) 
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
              cond_id=cond_id,
              sub_id=sub_id,
              gain=gain,
              loss=loss,
              cond=sort(unique(loss)),
              obs=obs,
              grid=grid,
              p_gain=p_gain,
              p_loss=p_loss,
              get_obs_rep=1,
              get_log_lik=0,
              prior_only=0)

# # function form 2 with an argument named `chain_id`
# initf2 <- function(chain_id = 1) {
#   # cat("chain_id =", chain_id, "\n")
#   list(beta_mu = 1 + rnorm(1,0,.05),
#        lambda_mu = 1 + rnorm(1,0,.05),
#        sigma_mu = .5 + rnorm(1,0,.05),
#        beta_sd = .5 + rnorm(1,0,.025),
#        lamdba_sd = .5 + rnorm(1,0,.025),
#        sigma_sd = .5 + rnorm(1,0,.025),
#        beta_re = rnorm(S,0,0.1),
#        lambda_re = rnorm(S,0,0.1),
#        sigma_re = rnorm(S,0,0.1))
# }
# # generate a list of lists to specify initial values
# n_chains <- 4
# init_ll <- lapply(1:n_chains, function(id) initf2(chain_id = id))

model <- stan_model(file="modeling/stan/scratch_model2.stan",
                    model_name="scratch_model2",
                    allow_undefined = TRUE,
                    includes = paste0('\n#include "',
                                      file.path(getwd(),
                                                'modeling/stan/hpp/get_max_idx.hpp'), 
                                      '"\n'))
fit <- sampling(
  object=model,
  data=input,
  # init=init_ll,
  chains=4,
  iter=4000,
  cores=4,
  control=list(max_treedepth=10))

obs_rep <- rstan::extract(fit, 'obs_rep')$obs_rep
launch_shinystan(fit)

color_scheme_set("darkgray")
mcmc_pairs(
  as.matrix(fit),
  pars = c("beta_mu", "lambda_mu", "sigma_mu"),
  regex_pars = "mu_cond",
  off_diag_fun = "hex")
