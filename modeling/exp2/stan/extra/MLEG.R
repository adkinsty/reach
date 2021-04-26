# Title     : maximum loss-averse expected gain model
# Objective : Interface with stan model
# Created by: adkinsty
# Updated on: 6/22/20

library(rstan)
library(tidyverse)
library(shinystan)
library(ggthemes)
library(raster)
library(mvtnorm)
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

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp2/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

la_df <- read_tsv("data/raw/loss_aversion/All_lambdas.csv") %>%
  filter(subNum %in% data$id) %>% arrange(subNum) %>% 
  mutate(id = subNum, la = lambdas)

# key variables
jj <- data$sub_id
cc <- data$cond_id
J <- length(unique(jj))
subs <- sort(unique(data$id))
N <- nrow(data)
C <- length(unique(cc))
cond <- sort(unique(data$loss))
y <- data$x

read_rds("data/clean/exp2/MLEG_aims.csv")

input <- list(N=N,
              J=J,
              C=C,
              cc=cc,
              jj=jj,
              aim=aim,
              y=y,
              get_gq=1,
              prior_only=0)

model <- stan_model(file="modeling/exp2/stan/aim.stan",model_name="MLEG")

mcmc_est <- sampling(
  object=model,
  data=input,
  chains=4,
  iter=5000,
  warmup=2500,
  cores=4)
write_rds(mcmc_est,path="modeling/exp2/stan/rds/MLEG.rds")









# 
# 
# print(mcmc_est,pars = c("sigma_m","sigma_s","sigma"),probs = c(.025,.975))
# stan_plot(mcmc_est,pars = c("sigma_m","sigma_s","sigma"))
# 
# yrep <- rstan::extract(mcmc_est, 'yrep')$yrep
# pred <- as.tibble(t(yrep)) %>% cbind(data) %>%
#   pivot_longer(cols=starts_with("V"),
#                names_to="draw",
#                values_to="yrep")
# 
# ppc_dens_overlay(y,yrep[sample(nrow(yrep), 25), ]) + theme(legend.position="top")
# ppc_ecdf_overlay(y, yrep[sample(nrow(yrep), 25), ]) + theme(legend.position="top")
# ppc_stat(y,yrep,stat = "mean") + theme(legend.position = "top")
# ppc_stat(y,yrep,stat = "sd") + theme(legend.position = "top")
# ppc_stat(y,yrep,stat = "max") + theme(legend.position = "top")
# ppc_stat(y,yrep,stat = "min") + theme(legend.position = "top")
# ppc_violin_grouped(y,yrep[sample(nrow(yrep), 25), ],group=data$ratio,
#                    y_draw = "points",y_alpha = .1)
# 
# # mu_obs vs mu_rep scatter
# pred %>% 
#   group_by(draw,id) %>% 
#   summarise(mu_rep=mean(yrep),
#             mu_obs=mean(x)) %>% ungroup() %>%
#   ggplot() +
#   coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
#   geom_abline(slope=1,intercept=0,colour="grey",linetype="dashed") +
#   stat_bin_hex(aes(x=mu_rep,y=mu_obs),binwidth = .05) +
#   scale_fill_viridis_b(option="A") +
#   theme_tufte(base_family = "sans",base_size=15) +
#   theme(axis.line = element_line(size=.25),
#         legend.position = "none") 
# 
# # sd_obs vs sd_rep scatter
# pred %>% 
#   group_by(draw,id) %>% 
#   summarise(sd_rep=sd(yrep),
#             sd_obs=sd(x)) %>%
#   ggplot() +
#   coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
#   geom_abline(slope=1,intercept=0,colour="grey",linetype="dashed") +
#   stat_bin_hex(aes(x=sd_rep,y=sd_obs),binwidth = .05) +
#   scale_fill_viridis_b(option="A") +
#   theme_tufte(base_family = "sans",base_size=15) +
#   theme(axis.line = element_line(size=.25),
#         legend.position = "none") 
# 
# # mu_obs vs mu_rep scatter by ratio
# pred %>% 
#   group_by(draw,id,ratio) %>% 
#   summarise(mu_rep=mean(yrep),
#             mu_obs=mean(x)) %>%
#   ggplot() +
#   coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
#   geom_abline(slope=1,intercept=0,colour="black",linetype="dashed") +
#   stat_bin_hex(aes(x=mu_rep,y=mu_obs),binwidth = .05) +
#   scale_fill_viridis_b(option="A") +
#   facet_wrap(.~ratio,nrow=3) +
#   theme_tufte(base_family = "sans",base_size=15) +
#   theme(axis.line = element_line(size=.25),
#         legend.position = "none") 
