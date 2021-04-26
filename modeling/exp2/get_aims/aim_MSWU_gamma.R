# Title     : Get aimpoints for MSWU gamma-restricted model
# Objective : inputs to stan model
# Created by: adkinsty
# Updated on: 6/29/20

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

library(tidyverse)
library(optimx)

get_u <- function(v,alpha) {
  if (v < 0){
    u <- -((-v)^alpha)
  } else {
    u <- v
  }
  return(u)
}
get_w <- function(p,gamma) {
  return(exp(-(-log(p))^gamma))
}
get_swu <- function(wp,gain,both,loss) {
  swu <- wp[1,]*gain + wp[2,]*both + wp[3,]*loss
  return(swu)
}
get_mu <- function(p,loss,alpha,gamma){
  hyp <- seq(0,1,.01) # aims
  gain <- 1
  uloss <- get_u(loss,alpha)
  uboth <- uloss + gain
  wp <- get_w(p,gamma)
  swu <- get_swu(wp,gain,uboth,uloss)
  mu <- hyp[which.max(swu)]
  return(mu)
}
get_obs <- function(data,j) {
  tmp = data %>% filter(sub_id == j)
  return(tmp$x)
}
get_loglik <- function(x,mu,sigma) {
  return(sum(-dnorm(x,mu,sigma,log=TRUE)))
}

train <- read_csv("data/clean/exp2/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

# forward model: outcome probabilities by aim-point
fm <- read_rds("data/clean/exp2/forward_model.rds")

# key variables
jj <- data$sub_id
cc <- data$cond_id
J <- length(unique(jj))
N <- nrow(data)
C <- length(unique(cc))
cond <- sort(unique(data$loss))
y <- data$x

param <- matrix(nrow=2, ncol=J) #  MLE parameters
mu_cond <- matrix(nrow=3,ncol=J) # MLE aim-points

# Maximum Likelihood Estimation
for (j in 1:J) {
  tmp <- data %>% filter(sub_id == j)
  obs <- tmp$x
  cond_id <- tmp$cond_id
  sig <- sd(obs)
  objective <- function(par,p=fm[,,j],y=obs,cc=cond_id,s=sig) {
    alpha=par[1]; gamma=par[2]
    cond <- c(-5, -1, 0)
    mc <- c()
    for (c in 1:length(cond)) {
      mc <- append(mc,get_mu(p,cond[c],alpha,gamma))
    }
    mu <- c()
    for (n in 1:length(y)) {
      mu[n] = mc[cc[n]]
    }
    sigma <- rep(s,length(y))
    ll <- get_loglik(y,mu,sigma)
    return(ll)
  }
  print(sprintf("Subject %s",j))
  fit <- optimx(par = c("alpha"=1,"gamma"=1),
                lower=c(0.001,0.001),upper=c(100,1.001),
                fn = objective,method="bobyqa")
  param[1,j] <- fit$alpha
  param[2,j] <- fit$gamma
  mu_cond[1,j] <- get_mu(p=fm[,,j],loss = -5, alpha = fit$alpha, gamma=fit$gamma)
  mu_cond[2,j] <- get_mu(p=fm[,,j],loss = -1, alpha = fit$alpha, gamma=fit$gamma)
  mu_cond[3,j] <- get_mu(p=fm[,,j],loss = 0, alpha = fit$alpha, gamma=fit$gamma)
}

write_rds(param,"data/clean/exp2/MSWU_gamma_param.rds")
write_rds(mu_cond,"data/clean/exp2/MSWU_gamma_aims.rds")