# Title     : Experiment 1, compute 'forward models' using training variances
# Objective : prepare input for models
# Created by: adkinsty
# Updated on: 9/10/20

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

library(shotGroups)
library(raster)
library(mvtnorm)
library(tidyverse)

get_prior_sigma <-function(test,j,ck){
  tmp <- test %>% filter(sub_id == j & ratio == ck)
  s <- array(c(var(tmp$x),cov(tmp$x,tmp$y),cov(tmp$x,tmp$y),var(tmp$y)),dim=c(2,2))
  return(s)
}
sim_bhv <- function(x,s,l_pos) {
  n=1e5 # num trials
  sim <- rmvnorm(n,mean=c(x,0),sigma=s)  # simulate data
  gpos <- cbind(rep(0,n),rep(0,n))  # gain position
  lpos <- cbind(rep(l_pos,n),rep(0,n)) # loss position
  gdist <- pointDistance(p1=sim,p2=gpos,lonlat = F) # dist from sim to gain
  ldist <- pointDistance(p1=sim,p2=lpos,lonlat = F) # dist from sim to loss
  return(cbind(gdist,ldist))
}
get_fm_vp <- function(test) {
  # gets outcome probabilites conditional on aimpoints
  # aka fm or forward models
  hyp <- seq(0,1,.01) # aims
  G = length(hyp) # num aims
  O <- 4 # num discrete outcome
  J <- length(unique(test$id)) # num subjects
  jj <- test$sub_id
  kk <- test$ratio
  cond <- sort(unique(kk))
  J <- length(unique(jj))
  K <- length(cond)
  fm <- array(dim=c(O,G,J,K,2)) # outcome probability container
  for (j in 1:J) {
    for (k in 1:K) {
      
      sigma <- get_prior_sigma(test,j,cond[k])
      
      get_p <- function(x,s=sigma,l_pos=-1) {
        dist <- sim_bhv(x=x,s=s,l_pos=l_pos)
        pg <- mean(dist[,1] <= 1 & dist[,2] > 1) # p of gain and no loss
        pl <- mean(dist[,1] > 1 & dist[,2] <= 1) # probability of loss and no gain
        pgl <- mean(dist[,1] <= 1 & dist[,2] <= 1) # p of gain and loss
        pn <- mean(dist[,1] > 1 & dist[,2] > 1) # p of neither gain nor loss
        return(matrix(data = c(pg,pgl,pl,pn),ncol = 1))
      }
      p_short <- sapply(X = hyp, FUN = get_p, simplify = T, USE.NAMES = F)
      get_p <- function(x,s=sigma,l_pos=-1.375) {
        dist <- sim_bhv(x=x,s=s,l_pos=l_pos)
        pg <- mean(dist[,1] <= 1 & dist[,2] > 1) # p of gain and no loss
        pl <- mean(dist[,1] > 1 & dist[,2] <= 1) # probability of loss and no gain
        pgl <- mean(dist[,1] <= 1 & dist[,2] <= 1) # p of gain and loss
        pn <- mean(dist[,1] > 1 & dist[,2] > 1) # p of neither gain nor loss
        return(matrix(data = c(pg,pgl,pl,pn),ncol = 1))
      }
      p_long <- sapply(X = hyp, FUN = get_p, simplify = T, USE.NAMES = F)
      
      fm[,,j,k,1] = p_short
      fm[,,j,k,2] = p_long
    }
  }
  return(fm)
}

data <- read_csv("data/clean/exp1/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

fm <- get_fm_vp(data) # forward model: outcome probabilities by aim-point

write_rds(fm,"data/clean/exp1/forward_model_vp.rds")

