library(rstan)
library(tidyverse)
library(shinystan)
library(ggthemes)
library(bayesplot)
library(raster)
library(mvtnorm)

get_p <- function(x) {
  dist <- sim_dist(x)
  pg <- mean(dist[,1] < 1 & dist[,2] > 1) # p of gain and no loss
  pl <- mean(dist[,1] > 1 & dist[,2] < 1) # probability of loss and no gain
  pgl <- mean(dist[,1] < 1 & dist[,2] < 1) # p of gain and loss
  pn <- mean(dist[,1] > 1 & dist[,2] > 1) # p of neither gain nor loss
  return(matrix(data = c(pg,pgl,pl,pn),ncol = 1))
}

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

losses <- c(-5)
alphas <- seq(0,1,.1)
sigmas <- c(.3)
lambdas <- c(1)
gammas <- seq(1,2,.1)

aim <- c(); los <- c(); alph <- c(); sig <- c(); lam <- c(); gam <- c()
hyp <- seq(0,1,.01) 

for (j in 1:length(sigmas)) {
  sigma <- array(c(sigmas[j],0,0,sigmas[j]),dim=c(2,2))
  sim_dist <- function(x,s=sigma,n=1e5) {
    sim <- rmvnorm(n,mean=c(x,0),sigma=s)  # simulate data
    gpos <- cbind(rep(0,n),rep(0,n))  # gain position
    lpos <- cbind(rep(-1,n),rep(0,n)) # loss position
    gdist <- pointDistance(p1=sim,p2=gpos,lonlat = F) # dist from sim to gain
    ldist <- pointDistance(p1=sim,p2=lpos,lonlat = F) # dist from sim to loss
    return(cbind(gdist,ldist))
  }
  op <- sapply(X = hyp, FUN = get_p, simplify = T, USE.NAMES = F)
  for (h in 1:length(gammas)) {
    p <- exp(-(-log(op))^gammas[h])
    for (i in 1:length(losses)) {
      for (k in 1:length(alphas)) {
        for (q in 1:length(lambdas)) {
          g <- 1
          l <- -lambdas[q]*(-losses[i])^alphas[k]
          gl <- -lambdas[q]*(-(g+losses[i]))^alphas[k]
          EU <- p[1,]*g + p[2,]*(g+l) + p[3,]*l+ p[4,]*0
          hyp[which.max(EU)]
          
          aim <- append(aim, hyp[which.max(EU)])
          los <- append(los, losses[i])
          alph <- append(alph, alphas[k])
          sig <- append(sig, sigmas[j])
          lam <- append(lam, lambdas[q])
          gam <- append(gam, gammas[h])
        }
      }
    }
  }
}

# low_v <- tibble(x=factor(los),y=aim,b=bet)
sim <- tibble(x=factor(los),y=aim,a=alph,la=lam,sd=sig,g=gam)

sim %>%  
  filter(a == 1) %>%
  ggplot(aes(y,x,group=sd,colour=factor(sd))) + 
  scale_y_discrete("Loss / Gain") +
  scale_x_continuous(expression(mu),limits=c(0,1)) +
  scale_colour_brewer(expression(sigma),palette = "Reds") +
  geom_col(width=.01,fill="white",position=position_dodge(.5)) +
  geom_point(size=3,position=position_dodge(.5)) + 
  theme_tufte(base_family = "sans",base_size = 12) +
  theme(legend.position = "top") 
ggsave("visuals/raw/sigma_shift.pdf",units="in",height=3,width=3)

sim %>%  
  filter(g==1) %>%
  ggplot(aes(y,x,group=a,colour=a)) + 
  scale_y_discrete("Loss / Gain") +
  scale_x_continuous(expression(mu),limits=c(0,1)) +
  # scale_colour_brewer(expression(alpha),palette = "Blues") +
  geom_col(width=.01,fill="white",position=position_dodge(.5)) +
  geom_point(size=3,position=position_dodge(.5)) + 
  theme_tufte(base_family = "sans",base_size = 12) +
  theme(legend.position = "top") 
ggsave("visuals/raw/alpha_shift.pdf",units="in",height=3,width=3)

sim %>%  
  filter(b==1 & sd ==.25 & g==1) %>%
  ggplot(aes(y,x,group=la,colour=factor(la))) + 
  scale_y_discrete("Loss / Gain") +
  scale_x_continuous(expression(mu),limits=c(0,1)) +
  scale_colour_brewer(expression(lambda),palette = "Greens") +
  geom_col(width=.01,fill="white",position=position_dodge(.5)) +
  geom_point(size=3,position=position_dodge(.5)) + 
  theme_tufte(base_family = "sans",base_size = 12) +
  theme(legend.position = "top") 
ggsave("visuals/raw/lambda_shift.pdf",units="in",height=3,width=3)

sim %>%  
  filter(b==1 & sd ==.25 & la==1) %>%
  ggplot(aes(y,x,group=g,colour=factor(g))) + 
  scale_y_discrete("Loss / Gain") +
  scale_x_continuous(expression(mu),limits=c(0,1)) +
  scale_colour_brewer(expression(gamma),palette = "Purples") +
  geom_col(width=.01,fill="white",position=position_dodge(.5)) +
  geom_point(size=3,position=position_dodge(.5)) + 
  theme_tufte(base_family = "sans",base_size = 12) +
  theme(legend.position = "top") 
ggsave("visuals/raw/gamma_shift.pdf",units="in",height=3,width=3)












