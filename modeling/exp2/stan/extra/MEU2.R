# Title     : maximum expected gain model
# Objective : Interface with stan model
# Created by: adkinsty
# Created on: 6/22/20

library(rstan)
library(tidyverse)
library(shinystan)
library(ggthemes)
library(shotGroups)
library(bayesplot)

get_pg <- function(x) {
  dist <- sim_dist(x)
  pg <- mean(dist[,1] < 1 & dist[,2] > 1) # p of gain and no loss
  return(pg)
}
get_pl <- function(x) {
  dist <- sim_dist(x)
  pl <- mean(dist[,1] > 1 & dist[,2] < 1) # probability of loss and no gain
  return(pl)
}
get_pgl <- function(x) {
  dist <- sim_dist(x)
  pgl <- mean(dist[,1] < 1 & dist[,2] < 1) # p of gain and loss
  return(pgl)
}
get_pn <- function(x) {
  dist <- sim_dist(x)
  pn <- mean(dist[,1] > 1 & dist[,2] > 1) # p of neither gain nor loss
  return(pn)
}
u <- function(x) {
  if (x>=0) {
    return(x)
  } else {
    return(-((-x)^.95))
  }
}


options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp2/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

# key variables
jj <- data$sub_id
cc <- data$cond_id
J <- length(unique(jj))
N <- nrow(data)
C <- length(unique(cc))
cond <- sort(unique(data$loss))
y <- data$x

# Get optimal aim points using training variances
hyp <- seq(0,1,.005) # approx -48 to 48 by 1 px
G <- length(grid)
aim <- matrix(nrow=3,ncol=J)
for (j in 1:J) {
  print(j)
  tmp <- train %>% filter(sub_id == j & block > 5)
  sigma <- array(c(var(tmp$x),cov(tmp$x,tmp$y),cov(tmp$x,tmp$y),var(tmp$y)),dim=c(2,2))
  
  sim_dist <- function(x,s=sigma,n=1e4) {
    sim <- rmvnorm(n,mean=c(x,0),sigma=s)  # simulate data
    gpos <- cbind(rep(0,n),rep(0,n))  # gain position
    lpos <- cbind(rep(-1,n),rep(0,n)) # loss position
    gdist <- pointDistance(p1=sim,p2=gpos,lonlat = F) # dist from sim to gain
    ldist <- pointDistance(p1=sim,p2=lpos,lonlat = F) # dist from sim to loss
    return(cbind(gdist,ldist))
  }

  pg <- sapply(X = hyp, FUN = get_pg, simplify = T, USE.NAMES = F)
  pl <- sapply(X = hyp, FUN = get_pl, simplify = T, USE.NAMES = F)
  pgl <- sapply(X = hyp, FUN = get_pgl, simplify = T, USE.NAMES = F)
  pn <- sapply(X = hyp, FUN = get_pn, simplify = T, USE.NAMES = F)
  
  for (c in 1:C) {
    ug = u(100)
    ul = u(cond[c]*100)
    ugl = ug+ul
    EU = pg*ug + pgl*ugl + pl*ul+ pn*0
    aim[c,j] = hyp[which.max(EU)]
  }
}

input <- list(N=N,
              J=J,
              C=C,
              cc=cc,
              jj=jj,
              aim=aim,
              y=y,
              get_yrep=1,
              get_log_lik=1,
              prior_only=0)

model <- stan_model(file="modeling/stan/MEU.stan",
                    model_name="MEU")

map_est = optimizing(object=model,data=input,as_vector=FALSE)
map_est$par
mcmc_est <- sampling(
  object=model,
  data=input,
  chains=4,
  iter=1500,
  warmup=1000,
  cores=4,
  control=list(max_treedepth=10))
write_rds(mcmc_est,path="modeling/stan/rds/MEU.rds")

print(mcmc_est,pars = c("sigma_m","sigma_s","sigma"),probs = c(.025,.975))
stan_plot(mcmc_est,pars = c("sigma_m","sigma_s","sigma"))

yrep <- rstan::extract(mcmc_est, 'yrep')$yrep
pred <- as_tibble(t(yrep)) %>% cbind(data) %>%
  pivot_longer(cols=starts_with("V"),
               names_to="draw",
               values_to="yrep")
color_scheme_set("red")
ppc_dens_overlay(y,yrep[sample(nrow(yrep), 25), ]) + theme(legend.position="top")
ppc_ecdf_overlay(y, yrep[sample(nrow(yrep), 25), ]) + theme(legend.position="top")
ppc_stat(y,yrep,stat = "mean") + theme(legend.position = "top")
ppc_stat(y,yrep,stat = "sd") + theme(legend.position = "top")
ppc_stat(y,yrep,stat = "max") + theme(legend.position = "top")
ppc_stat(y,yrep,stat = "min") + theme(legend.position = "top")
ppc_violin_grouped(y,yrep[sample(nrow(yrep), 25), ],group=data$ratio,
                   y_draw = "points",y_alpha = .2,y_jitter = .2,y_size = .2)

loo_est <- loo(mcmc_est, save_psis = TRUE, cores = 4)
psis <- loo_est$psis_object
lw <- weights(psis)
ppc_loo_pit_overlay(y, yrep, lw = lw) + theme(legend.position = "top")
ppc_loo_pit_qq(y, yrep, lw = lw) + theme(legend.position = "top")

keep_obs <- sample(nrow(yrep), 50)
ppc_loo_intervals(y, yrep, psis_object = psis, subset = keep_obs,order = "median")
ppc_loo_ribbon(y, yrep, psis_object = psis, subset = keep_obs)

# mu_obs vs mu_rep scatter
pred %>% 
  group_by(draw,id) %>% 
  summarise(mu_rep=mean(yrep),
            mu_obs=mean(x)) %>% ungroup() %>%
  ggplot() +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
  geom_abline(slope=1,intercept=0,colour="grey",linetype="dashed") +
  stat_bin_hex(aes(x=mu_rep,y=mu_obs),binwidth = .05) +
  scale_fill_viridis_b(option="A") +
  theme_tufte(base_family = "sans",base_size=15) +
  theme(axis.line = element_line(size=.25),
        legend.position = "none") 

# sd_obs vs sd_rep scatter
pred %>% 
  group_by(draw,id) %>% 
  summarise(sd_rep=sd(yrep),
            sd_obs=sd(x)) %>%
  ggplot() +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
  geom_abline(slope=1,intercept=0,colour="grey",linetype="dashed") +
  stat_bin_hex(aes(x=sd_rep,y=sd_obs),binwidth = .05) +
  scale_fill_viridis_b(option="A") +
  theme_tufte(base_family = "sans",base_size=15) +
  theme(axis.line = element_line(size=.25),
        legend.position = "none") 

# mu_obs vs mu_rep scatter by ratio
pred %>% 
  group_by(draw,id,ratio) %>% 
  summarise(mu_rep=mean(yrep),
            mu_obs=mean(x)) %>%
  ggplot() +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed") +
  stat_bin_hex(aes(x=mu_rep,y=mu_obs),binwidth = .05) +
  scale_fill_viridis_b(option="A") +
  facet_wrap(.~ratio,nrow=3) +
  theme_tufte(base_family = "sans",base_size=15) +
  theme(axis.line = element_line(size=.25),
        legend.position = "none") 

# sd_obs vs sd_rep scatter by ratio
pred %>% 
  group_by(draw,id,ratio) %>% 
  summarise(sd_rep=sd(yrep),
            sd_obs=sd(x)) %>%
  ggplot() +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed") +
  stat_bin_hex(aes(x=sd_rep,y=sd_obs),binwidth = .05,alpha=.8) +
  scale_fill_viridis_b(option="A") +
  facet_wrap(.~ratio,nrow=3) +
  theme_tufte(base_family = "sans",base_size=15) +
  theme(axis.line = element_line(size=.25),
        legend.position = "none") 
