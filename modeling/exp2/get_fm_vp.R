# Title     : Experiment 2, Forward model with variable precision
# Objective : get probabilistic consequences of different aimpoints given training variance
# Created by: adkinsty
# Updated on: 9/6/20

library(tidyverse)
library(mvtnorm)
library(raster)

get_sigma <-function(dat,j,k){
  tmp <- dat %>% filter(sub_id == j & cond_id == k)
  sigma <- array(c(var(tmp$x),cov(tmp$x,tmp$y),cov(tmp$x,tmp$y),var(tmp$y)),dim=c(2,2))
  return(sigma)
}
sim_bhv <- function(x,s) {
  n=1e5 # num trials
  sim <- rmvnorm(n,mean=c(x,0),sigma=s)  # simulate dat
  gpos <- cbind(rep(0,n),rep(0,n))  # gain position
  lpos <- cbind(rep(-1,n),rep(0,n)) # loss position
  gdist <- pointDistance(p1=sim,p2=gpos,lonlat = F) # dist from sim to gain
  ldist <- pointDistance(p1=sim,p2=lpos,lonlat = F) # dist from sim to loss
  return(cbind(gdist,ldist))
}
get_fm <- function(dat) {
  # gets outcome probabilites conditional on aimpoints
  # aka fm or forward models
  hyp <- seq(0,1,.01) # aims
  G = length(hyp) # num aims
  O <- 4 # num discrete outcome
  J <- length(unique(dat$id)) # num subjects
  K <- length(unique(dat$cond_id)) # num conditions
  fm <- array(dim=c(O,G,J,K)) # outcome probability container
  for (j in 1:J) {
    for (k in 1:K) {
      sigma <- get_sigma(dat,j,k)
      get_p <- function(x,s=sigma) {
        dist <- sim_bhv(x=x,s=s)
        pg <- mean(dist[,1] < 1 & dist[,2] > 1) # p of gain and no loss
        pl <- mean(dist[,1] > 1 & dist[,2] < 1) # probability of loss and no gain
        pgl <- mean(dist[,1] < 1 & dist[,2] < 1) # p of gain and loss
        pn <- mean(dist[,1] > 1 & dist[,2] > 1) # p of neither gain nor loss
        return(matrix(data = c(pg,pgl,pl,pn),ncol = 1))
      }
      p <- sapply(X = hyp, FUN = get_p, simplify = T, USE.NAMES = F)
      fm[,,j,k] = p
    }
  }
  return(fm)
}

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

test <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

fm <- get_fm(test)

write_rds(fm,"data/clean/exp2/forward_model_vp.rds")

