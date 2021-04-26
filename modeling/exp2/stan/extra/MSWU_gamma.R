# Title     : Experiment 2, maximum subjective weighted utility model
# Objective : Interface with stan model
# Created by: adkinsty
# Updated on: 6/29/20

library(rstan)
library(tidyverse)
library(shinystan)
library(ggthemes)
library(shotGroups)
library(bayesplot)
library(optimx)
library(mvtnorm)
library(raster)

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp2/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

get_prior_sigma <-function(train=train,j){
  tmp <- train %>% filter(sub_id == j & block > 5)
  sigma <- array(c(var(tmp$x),cov(tmp$x,tmp$y),cov(tmp$x,tmp$y),var(tmp$y)),dim=c(2,2))
  return(sigma)
}
sim_bhv <- function(x,sigma) {
  n=1e5 # num trials
  sim <- rmvnorm(n,mean=c(x,0),sigma=sigma)  # simulate data
  gpos <- cbind(rep(0,n),rep(0,n))  # gain position
  lpos <- cbind(rep(-1,n),rep(0,n)) # loss position
  gdist <- pointDistance(p1=sim,p2=gpos,lonlat = F) # dist from sim to gain
  ldist <- pointDistance(p1=sim,p2=lpos,lonlat = F) # dist from sim to loss
  return(cbind(gdist,ldist))
}
get_fm <- function(train,test) {
  # gets outcome probabilites conditional on aimpoints
  # aka fm or forward models
  hyp <- seq(0,1,.01) # aims
  G = length(hyp) # num aims
  O <- 4 # num discrete outcome
  J <- length(unique(data$id)) # num subjects
  fm <- array(dim=c(O,G,J)) # outcome probability container
  for (j in 1:length(unique(data$id))) {
    sigma <- get_prior_sigma(train,j)
    get_p <- function(x,s=sigma) {
      dist <- sim_bhv(x=x,s=s)
      pg <- mean(dist[,1] < 1 & dist[,2] > 1) # p of gain and no loss
      pl <- mean(dist[,1] > 1 & dist[,2] < 1) # probability of loss and no gain
      pgl <- mean(dist[,1] < 1 & dist[,2] < 1) # p of gain and loss
      pn <- mean(dist[,1] > 1 & dist[,2] > 1) # p of neither gain nor loss
      return(matrix(data = c(pg,pgl,pl,pn),ncol = 1))
    }
    p <- sapply(X = hyp, FUN = get_p, simplify = T, USE.NAMES = F)
    fm[,,j] = p
  }
  return(fm)
}
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

fm <- get_fm(train,data) # forward model: outcome probabilities by aim-point

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
  fit <- optimx(par = c("alpha"=.5,"gamma"=.5),
                lower=c(0.001,0.001),upper=c(100,1),
                fn = objective,method="bobyqa")
  param[1,j] <- fit$alpha
  param[2,j] <- fit$gamma
  mu_cond[1,j] <- get_mu(p=fm[,,j],loss = -5, alpha = fit$alpha, gamma=fit$gamma)
  mu_cond[2,j] <- get_mu(p=fm[,,j],loss = -1, alpha = fit$alpha, gamma=fit$gamma)
  mu_cond[3,j] <- get_mu(p=fm[,,j],loss = 0, alpha = fit$alpha, gamma=fit$gamma)
}

write.csv(param,"data/clean/exp2/MSWU_gamma_param.csv")
write.csv(mu_cond,"data/clean/exp2/MSWU_gamma_aims.csv")


input <- list(N=N,
              J=J,
              C=C,
              cc=cc,
              jj=jj,
              aim=mu_cond,
              y=y,
              get_gq=1,
              prior_only=0)

model <- stan_model(file="modeling/exp2/stan/aim.stan",model_name="MSWU")

mcmc_est <- sampling(
  object=model,
  data=input,
  chains=4,
  iter=5000,
  warmup=2500,
  cores=4)
write_rds(mcmc_est,path="modeling/exp2/stan/rds/MSWU_gamma.rds")




la_df <- read_tsv("data/raw/loss_aversion/All_lambdas.csv") %>%
  filter(subNum %in% data$id) %>% arrange(subNum) %>% 
  mutate(id = subNum, la = lambdas)



# parameters
la <- la_df$lambdas
la <- ifelse(la < -10 | la > -1 | is.na(la), NA,abs(la))
s <- c()
for (j in 1:J) {
  tmp <- train %>% filter(sub_id == j)
  s[j] <- sd(tmp$x)
}

tibble(est=c(param[1,],param[2,]),
       par=c(rep("a",J),rep("g",J))) %>%
  ggplot(aes(x=factor(par,levels=c("s","a","g")),y=est)) +
  geom_violin(draw_quantiles = TRUE) +
  geom_point(alpha=.5,position=position_jitter(.05),colour="#0065cc") +
  stat_summary(geom="point",fun.y = "median",size=2) +
  geom_hline(yintercept = 1, linetype="dashed") +
  scale_x_discrete(
    "Parameter",
    labels=c(expression(alpha),
             expression(gamma))) +
  scale_y_continuous("Estimate") +
  theme_tufte(base_size = 12, base_family = "sans") +
  theme(axis.line = element_line(size=.25))

ggsave("visuals/raw/exp2_mswu_gamma_par_est.pdf",units = "in",height = 2.4,width = 3.5)







  
