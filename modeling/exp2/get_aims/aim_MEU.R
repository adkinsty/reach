# Title     : Get aimpoints for MEU model
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
get_eu <- function(p,gain,both,loss) {
  eu <- p[1,]*gain + p[2,]*both + p[3,]*loss
  return(eu)
}
get_mu <- function(p,loss,alpha,beta){
  hyp <- seq(0,1,.01) # aims
  gain <- 1
  uloss <- get_u(loss,alpha)
  ugain <- get_u(gain,beta)
  uboth <- uloss + ugain
  eu <- get_eu(p,gain,uboth,uloss)
  mu <- hyp[which.max(eu)]
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
kk <- data$cond_id
J <- length(unique(jj))
K <- length(unique(kk))

param <-  matrix(nrow=2,ncol=J) #  MLE parameters
mu_cond <- matrix(nrow=3,ncol=J) # MLE aim-points

# Maximum Likelihood Estimation
for (j in 1:J) {
  tmp <- data %>% filter(sub_id == j)
  obs <- tmp$x
  cond_id <- tmp$cond_id
  sig <- sd(obs)
  objective <- function(par,p=fm[,,j],y=obs,kk=cond_id,s=sig) {
    alpha <- par[1]; beta <- par
    uk <- c(-5, -1, 0)
    mk <- c()
    for (k in 1:length(uk)) {
      mk <- append(mk,get_mu(p,uk[k],alpha,beta))
    }
    mu <- c()
    for (n in 1:length(y)) {
      mu[n] = mk[kk[n]]
    }
    sigma <- rep(s,length(y))
    ll <- get_loglik(y,mu,sigma)
    return(ll)
  }
  print(sprintf("Subject %s",j))
  fit <- optimx(par = c("alpha"=1,"beta"=1),
                lower=c(0.001,0.001),upper=c(100,100),
                fn = objective,method="bobyqa")
  param[1,j] <- fit$alpha
  param[2,j] <- fit$beta

  mu_cond[1,j] <- get_mu(p=fm[,,j],loss = -5, param[1,j], param[2,j])
  mu_cond[2,j] <- get_mu(p=fm[,,j],loss = -1, param[1,j], param[2,j])
  mu_cond[3,j] <- get_mu(p=fm[,,j],loss = 0,  param[1,j], param[2,j])
}
2 
write_rds(param,"data/clean/exp2/MEU_param.rds")
write_rds(mu_cond,"data/clean/exp2/MEU_aims.rds")
