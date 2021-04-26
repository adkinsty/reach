# Title     : heuristic
# Objective : Interface with stan model
# Created by: adkinsty
# Created on: 6/25/20

library(rstan)
library(tidyverse)
library(shinystan)
library(ggthemes)
library(shotGroups)
library(bayesplot)

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

aim <- matrix(nrow=3,ncol=J)
for (j in 1:J) {
  # calculate win probabilities and get aim with max
  for (c in 1:C) {
    if (cond[c]==0) {
      aim[c,j] = 0
    } else {
      aim[c,j] = .5
    }
  }
}
write.csv(aim,"data/clean/exp2/H_aims.csv")


input <- list(N=N,J=J,C=C,
              cc=cc,jj=jj,
              aim=aim,
              y=y,
              get_yrep=1,
              get_log_lik=1,
              prior_only=0)

model <- stan_model(file="modeling/exp2/stan/aim.stan",model_name="H")

map_est = optimizing(object=model,data=input,as_vector=FALSE)
map_est$par$sigma_m; map_est$par$sigma_s; map_est$par$sigma

# mcmc_est <- readRDS("modeling/stan/rds/H.rds")
mcmc_est <- sampling(
  object=model,
  data=input,
  chains=4,
  iter=5000,
  warmup=2500,
  cores=4)
write_rds(mcmc_est,path="modeling/exp2/stan/rds/H.rds")










# 
# 
# print(mcmc_est,pars = c("sigma_m","sigma_s","sigma"),probs = c(.025,.975))
# stan_plot(mcmc_est,pars = c("sigma_m","sigma_s","sigma"))
# 
# yrep <- rstan::extract(mcmc_est, 'yrep')$yrep
# pred <- as_tibble(t(yrep)) %>% cbind(data) %>%
#   pivot_longer(cols=starts_with("V"),names_to="draw",values_to="yrep")
# color_scheme_set("red")
# ppc_dens_overlay(y,yrep[sample(nrow(yrep), 25), ]) + theme(legend.position="top")
# ppc_ecdf_overlay(y, yrep[sample(nrow(yrep), 25), ]) + theme(legend.position="top")
# ppc_stat(y,yrep,stat = "mean") + theme(legend.position = "top")
# ppc_stat(y,yrep,stat = "sd") + theme(legend.position = "top")
# ppc_stat(y,yrep,stat = "max") + theme(legend.position = "top")
# ppc_stat(y,yrep,stat = "min") + theme(legend.position = "top")
# ppc_violin_grouped(y,yrep[sample(nrow(yrep), 25), ],group=data$ratio,
#                    y_draw = "points",y_alpha = .2,y_jitter = .2,y_size = .2)
# 
# loo_est <- loo(mcmc_est, save_psis = TRUE, cores = 4)
# psis <- loo_est$psis_object
# lw <- weights(psis)
# ppc_loo_pit_overlay(y, yrep, lw = lw) + theme(legend.position = "top")
# ppc_loo_pit_qq(y, yrep, lw = lw) + theme(legend.position = "top")
# 
# keep_obs <- sample(nrow(yrep), 50)
# ppc_loo_intervals(y, yrep, psis_object = psis, subset = keep_obs,order = "median")
# ppc_loo_ribbon(y, yrep, psis_object = psis, subset = keep_obs)
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
# 
# # sd_obs vs sd_rep scatter by ratio
# pred %>% 
#   group_by(draw,id,ratio) %>% 
#   summarise(sd_rep=sd(yrep),
#             sd_obs=sd(x)) %>%
#   ggplot() +
#   coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
#   geom_abline(slope=1,intercept=0,colour="black",linetype="dashed") +
#   stat_bin_hex(aes(x=sd_rep,y=sd_obs),binwidth = .05,alpha=.8) +
#   scale_fill_viridis_b(option="A") +
#   facet_wrap(.~ratio,nrow=3) +
#   theme_tufte(base_family = "sans",base_size=15) +
#   theme(axis.line = element_line(size=.25),
#         legend.position = "none") 
