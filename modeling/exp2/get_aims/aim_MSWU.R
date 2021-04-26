# Title     : Get aimpoints for MSWU model
# Objective : inputs to stan model
# Created by: adkinsty
# Updated on: 3/11/21

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

library(tidyverse)
library(optimx)

get_u <- function(v,w) {
  # non-linear value transform
  if (v < 0){
    u <- -((-v)^w)
  } else {
    if (v > 0) {
      u <- v^w
    } else {
      u <- 0
    }
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
get_mu <- function(p,loss,alpha,beta,gamma){
  hyp <- seq(0,1,.01) # aims
  gain <- 1
  uloss <- get_u(loss,alpha)
  ugain <- get_u(gain,beta)
  uboth <- uloss + ugain
  wp <- get_w(p,gamma)
  swu <- get_swu(wp,ugain,uboth,uloss)
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
K <- length(unique(cc))
cond <- sort(unique(data$loss))
y <- data$x

param <- matrix(nrow=3, ncol=J) #  MLE parameters
mu_cond <- matrix(nrow=3,ncol=J) # MLE aim-points

 # Maximum Likelihood Estimation
for (j in 1:J) {
  tmp <- data %>% filter(sub_id == j)
  obs <- tmp$x
  cond_id <- tmp$cond_id
  sig <- sd(obs)
  objective <- function(par,p=fm[,,j],y=obs,kk=cond_id,s=sig) {
    alpha=par[1]; beta=par[2]; gamma=par[3]
    uk <- c(-5, -1, 0)
    mk <- c()
    for (k in 1:length(uk)) {
      mk <- append(mk,get_mu(p,uk[k],alpha,beta,gamma))
    }
    mu <- c()
    for (i in 1:length(y)) {
      mu[i] = mk[kk[i]]
    }
    sigma <- rep(s,length(y))
    ll <- get_loglik(y,mu,sigma)
    return(ll)
  }
  print(sprintf("Subject %s",j))
  fit <- optimx(par = c("alpha"=1,"beta"=1,"gamma"=1),
                lower=c(0.001,0.001,0.001),upper=c(100,100,100),
                fn = objective,method="bobyqa")
  param[1,j] <- fit$alpha
  param[2,j] <- fit$beta
  param[3,j] <- fit$gamma
  mu_cond[1,j] <- get_mu(p=fm[,,j],loss = -5, param[1,j], param[2,j], param[3,j])
  mu_cond[2,j] <- get_mu(p=fm[,,j],loss = -1, param[1,j], param[2,j], param[3,j])
  mu_cond[3,j] <- get_mu(p=fm[,,j],loss = 0,  param[1,j], param[2,j], param[3,j])
}

write_rds(param,"data/clean/exp2/MSWU_param.rds")
write_rds(mu_cond,"data/clean/exp2/MSWU_aims.rds")
