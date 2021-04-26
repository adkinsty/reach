# Title     : universal optimal heuristic model
# Objective : Interface with stan model
# Created by: adkinsty
# Updated on: 2/18/21

library(tidyverse)
library(optimx)

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp1/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

data <- read_csv("data/clean/exp1/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

fm <- read_rds("data/clean/exp1/forward_model.rds")

get_l <- function(k) {
  # get loss value based on cond id
  if (k %in% c(1,6)) {
    l <- -15
  } else if (k %in% c(2, 7)) {
    l <- -5
  } else if (k %in% c(3, 8)) {
    l <- -3
  } else if (k %in% c(4, 9)) {
    l <- -1
  } else {
    l <- 0
  }
  return(l)
}
get_g <- function(l) {
  # get gain value based on loss 
  if (l %in% c(-15, -3)) {
    g <- 3
  } else {
    g <- 1
  }
  return(g)
}
get_p <- function(k,P) {
  if (k < 6) {
    p <- P[,,1]
  } else {
    p <- P[,,2]
  }
  return(p)
}
get_mu <- function(p,k,intercept,slope){
  hyp <- seq(0,1,.01) # aims
  # get outcome probs based on separation distance
  pg <- p[1,]; pgl <- p[2,];  pl <- p[3,]; pn <- p[4,]
  l <- get_l(k)
  g <- get_g(l)
  ratio <- abs(l / g)
  v <- ifelse(ratio == 0, -.5, ifelse(ratio == 1, 0, .5))
  mu <- intercept + slope * v
  return(mu)
}
get_loglik <- function(x,mu,sigma) {
  return(sum(-dnorm(x,mu,sigma,log=TRUE)))
}

# key variables
N <- nrow(data)
jj <- data$sub_id
J <- length(unique(jj))
cc <- data$cond_id
K <- length(unique(cc))
y <- data$x

# Grid optimization to obtain predicted aim points (means)
hyp <- seq(0,1,.01) 
aim <- matrix(nrow=K,ncol=J)
short_EU <- c()
long_EU <- c()
for (j in 1:J) {
  tmp <- data %>% filter(sub_id == j)
  
  all_p = fm[,,j,]
    
  # calculate win probabilities and get aim with max
  for (k in 1:K) {
    
    if (k < 6) {
      p <- all_p[,,1]
    } else {
      p <- all_p[,,2]
    }
    
    # get outcome probs based on separation distance
    pg <- p[1,]; pgl <- p[2,];  pl <- p[3,]; pn <- p[4,]
    
    # get loss value based on cond id
    if (k %in% c(1,6)) {
      l <- -15
    } else if (k %in% c(2, 7)) {
      l <- -5
    } else if (k %in% c(3, 8)) {
      l <- -3
    } else if (k %in% c(4, 9)) {
      l < -1
    } else {
      l <- 0
    }
    
    # get gain value based on loss
    if (l %in% c(-15, -3)) {
      g <- 3
    } else {
      g <- 1
    }
    
    gl = g + l
    
    EU = pg*g + pgl*gl + pl*l+ pn*0

    if (k < 6) {
      short_EU <- rbind(short_EU, EU)
    } else {
      long_EU <- rbind(long_EU, EU)
    }
    
  }
}

# universal optimal
short_UEU <- colMeans(short_EU)
long_UEU <- colMeans(long_EU)
init_aim <- c(hyp[which.max(short_UEU)], hyp[which.max(long_UEU)])

param <- c() #  MLE parameters
mu_cond <- matrix(nrow=K,ncol=J) # MLE aim-points

# Maximum Likelihood Estimation
for (j in 1:J) {
  tmp <- data %>% filter(sub_id == j)
  obs <- tmp$x
  cond_id <- tmp$cond_id
  sig <- sd(obs)
  all_p = fm[,,j,]
  
  objective <- function(par,P=all_p,y=obs,kk=cond_id,s=sig,intercepts = init_aim) {
    slope=par
    mk <- c()
    for (k in 1:K) {
      p <- get_p(k,P)
      if (k < 6) {
        intercept <- intercepts[1]
      } else {
        intercept <- intercepts[2]
      }
      mk <- append(mk,get_mu(p,k,intercept,slope))
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
  fit <- optimx(par = c("slope"=-.5),
                lower=c(-1),upper=c(1),
                fn = objective,method="bobyqa")
  # results
  param[j] <- fit$slope
  for (k in 1:K) {
    if (k < 6) {
      intercept <- init_aim[1]
    } else {
      intercept <- init_aim[2]
    }
    p <- get_p(k, all_p)
    mu_cond[k,j] <- get_mu(p,k,intercept,fit$slope)
  }
}

write_rds(param,"data/clean/exp1/ULH_param.rds")
write_rds(mu_cond,"data/clean/exp1/ULH_aims.rds")
