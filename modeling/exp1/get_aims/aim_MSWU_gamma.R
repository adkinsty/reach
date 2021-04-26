# Title     : Uses MLE to obtain predictions for the MSWU model
# Objective : generates input to stan model
# Created by: adkinsty
# Updated on: 8/7/20

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
  swu <- wp[1,]*gain + wp[2,]*(gain+loss) + wp[3,]*loss
  return(swu)
}
get_mu <- function(p,c,alpha,gamma){
  hyp <- seq(0,1,.01) # aims
  
  # get outcome probs based on separation distance
  pg <- p[1,]; pgl <- p[2,];  pl <- p[3,]; pn <- p[4,]
  
  # get loss value based on cond id
  if (c %in% c(1,6)) {
    l <- -15
  } else if (c %in% c(2, 7)) {
    l <- -5
  } else if (c %in% c(3, 8)) {
    l <- -3
  } else if (c %in% c(4, 9)) {
    l <- -1
  } else {
    l <- 0
  }
  print(l)
  
  # get gain value based on loss 
  if (l %in% c(-15, -3)) {
    g <- 3
  } else {
    g <- 1
  }
  
  uloss <- get_u(l,alpha)
  uboth <- g + uloss
  
  wp <- get_w(p,gamma)
  swu <- get_swu(wp,g,uboth,uloss)
  mu <- hyp[which.max(swu)]
  return(mu)
}
get_loglik <- function(x,mu,sigma) {
  return(sum(-dnorm(x,mu,sigma,log=TRUE)))
}

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp1/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

data <- read_csv("data/clean/exp1/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

fm <- read_rds("data/clean/exp1/forward_model.rds")

# key variables
jj <- data$sub_id
cc <- data$cond_id
J <- length(unique(jj))
N <- nrow(data)
C <- length(unique(cc))
cond <- sort(unique(data$cond_id))
y <- data$x

param <- matrix(nrow=2, ncol=J) #  MLE parameters
mu_cond <- matrix(nrow=C,ncol=J) # MLE aim-points

# Maximum Likelihood Estimation
for (j in 1:J) {
  tmp <- data %>% filter(sub_id == j)
  obs <- tmp$x
  cond_id <- tmp$cond_id
  sig <- sd(obs)
  all_p = fm[,,j,]
  objective <- function(par,P=all_p,y=obs,cc=cond_id,s=sig) {
    alpha=par[1]; gamma=par[2]
    cond <- 1:10
    mc <- c()
    for (i in cond) {
      if (i < 6) {
        p <- P[,,1]
      } else {
        p <- P[,,2]
      }
      mc <- append(mc,get_mu(p,i,alpha,gamma))
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
  for (c in 1:10) {
    if (c < 6) {
      p <- all_p[,,1]
    } else {
      p <- all_p[,,2]
    }
    mu_cond[c,j] <- get_mu(p,c,fit$alpha,fit$gamma)
  }
}

write_rds(param,"data/clean/exp1/MSWU_gamma_param.rds")
write_rds(mu_cond,"data/clean/exp1/MSWU_gamma_aims.rds")
