# Title     : universal heuristic model
# Objective : get aimpoints for uh model (input to stan model)
# Created by: adkinsty
# Updated on: 2/5/21

library(tidyverse)
library(optimx)

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp2/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

fm <- read_rds("data/clean/exp2/forward_model.rds")

get_mu <- function(loss,intercept,slope){
  gain <- 1
  ratio <- abs(loss / gain)
  v <- ifelse(ratio == 0, 1, 0)
  mu <- intercept + slope * v
  return(mu)
}
get_obs <- function(data,j) {
  tmp = data %>% filter(sub_id == j)
  return(tmp$x)
}
get_loglik <- function(x,mu,sigma) {
  return(sum(-dnorm(x,mu,sigma,log=TRUE)))
}

# key variables
jj <- data$sub_id
kk <- data$cond_id
J <- length(unique(jj))
K <- length(unique(kk))
cond <- sort(unique(data$loss))
y <- data$x

# compute expected value lanscapes
hyp <- seq(0,1,.01) 

short_EU <- c()
long_EU <- c()

for (j in 1:J) {
  
  p <- fm[,,j]
  
  pg <- p[1,]; pgl <- p[2,]; pl <- p[3,]; pn <- p[4,]
  
  for (k in 1:K) {
    g <- 1
    l <- cond[k]
    gl <- g + l
    EU <- pg*g + pgl*gl + pl*l+ pn*0
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

param <- matrix(nrow=3, ncol=J) #  MLE parameters
mu_cond <- matrix(nrow=3,ncol=J) # MLE aim-points

# Maximum Likelihood Estimation
for (j in 1:J) {
  tmp <- data %>% filter(sub_id == j)
  obs <- tmp$x
  cond_id <- tmp$cond_id
  sig <- sd(obs)
  objective <- function(par,y=obs,kk=cond_id,s=sig,intercepts=init_aim) {
    uk <- c(-5, -1, 0)
    mk <- c()
    for (k in 1:length(uk)) {
      if (k < 6) {
        intercept <- intercepts[1]
      } else {
        intercept <- intercepts[0]
      }
      mk <- append(mk,get_mu(uk[k], intercept=intercept, slope=par['slope']))
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
  fit <- optimx(par = c("slope"=-0.001),lower = c(-1), upper = c(1),
                fn = objective, method = "bobyqa")
  param[1,j] <- fit$slope
  param[2,j] <- init_aim[1]
  param[3,j] <- init_aim[2]

  mu_cond[1,j] <- get_mu(loss = -5, intercept = init_aim, slope = fit$slope)
  mu_cond[2,j] <- get_mu(loss = -1, intercept = init_aim, slope = fit$slope)
  mu_cond[3,j] <- get_mu(loss = 0, intercept = init_aim, slope = fit$slope)
}
write_rds(param,"data/clean/exp2/UH_param.rds")
write_rds(mu_cond,"data/clean/exp2/UH_aims.rds")
