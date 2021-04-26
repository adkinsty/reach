# Title     : Linear regression model
# Objective : Interface with stan model
# Created by: adkinsty
# Created on: 6/22/20

library(rstan)
library(tidyverse)
library(shinystan)
library(ggthemes)
library(shotGroups)
library(bayesplot)
library(tidyr)

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2) %>%
  mutate(z_cond = ifelse(ratio==0,-1,ifelse(ratio==-1,0,1)))

la_df <- read_tsv("data/raw/loss_aversion/All_lambdas.csv") %>%
  filter(subNum %in% data$id) %>% arrange(subNum) %>% 
  mutate(id = 1:n(),
         z_la = scale(lambdas)[,1])

input <- list(
  N=nrow(data),                                      # num trials
  K=2,                                               # num trial predictors
  J=length(unique(data$sub_id)),                     # num subjects
  L=1,                                               # num subject predictors
  jj=data$sub_id,                                    # subject for trial
  x=as.matrix(cbind(rep(1,nrow(data)),data$z_cond)),  # trial predictors
  u=as.matrix(cbind(rep(1,nrow(la_df)))), # subject predictors
  y=data$x,                                          # observations
  get_yrep=1,  
  get_log_lik=1,
  prior_only=0
)

model <- stan_model(
  file="modeling/stan/LR.stan",
  model_name="LR"
)


mcmc_est <- sampling(
  object=model,
  data=input,
  chains=4,
  iter=5000,
  warmup=2500,
  cores=4,
  control=list(max_treedepth=10))
write_rds(mcmc_est,path="modeling/stan/rds/LR.rds")


s <- rstan::summary(mcmc_est,pars=c("beta","sigma"),prob=c(.25,.975))

stan_plot(mcmc_est,pars = c("sigma_m","sigma_s","sigma")) + vline_at(1)

print(mcmc_est,pars = c("gamma","beta"),probs = c(.025,.05,.10,.90,.95,.975))
stan_plot(mcmc_est,pars = c("gamma","beta")) + vline_0()

yrep <- rstan::extract(mcmc_est, 'yrep')$yrep
pred <- as.tibble(t(yrep)) %>% cbind(data) %>%
  pivot_longer(cols=starts_with("V"),
               names_to="draw",
               values_to="yrep")


# mu_obs vs mu_rep histogram
pred %>% 
  group_by(draw) %>% 
  summarise(mu_rep=mean(yrep)) %>%
  ggplot() +
  geom_histogram(aes(x=mu_rep),binwidth=.003,colour="white",fill="grey") +
  geom_density(aes(x=mu_rep),colour="grey",size=1) +
  geom_vline(xintercept=mean(data$x),size=1,linetype="solid",colour="green3") + 
  theme_tufte(base_family = "sans",base_size=15) +
  theme(axis.line = element_line(size=.25))


# sd_obs vs sd_rep histogram
pred %>% 
  group_by(draw) %>% 
  summarise(sd_rep=sd(yrep)) %>%
  ggplot() +
  geom_histogram(aes(x=sd_rep),binwidth = .003,colour="white",fill="grey") +
  geom_density(aes(x=sd_rep),colour="grey",size=1) +
  geom_vline(xintercept=sd(data$x),size=1,linetype="solid",colour="green3") + 
  theme_tufte(base_family = "sans",base_size=15) +
  theme(axis.line = element_line(size=.25))


# mu_obs vs mu_rep by ratio histogram
pred %>% 
  group_by(draw,ratio) %>% 
  summarise(mu_rep=mean(yrep)) %>%
  ggplot() +
  geom_histogram(aes(x=mu_rep),binwidth = .01,colour="white",fill="grey") +
  geom_vline(data=data%>%group_by(ratio)%>%summarise(mu=mean(x)),
             aes(xintercept=mu),colour="green3",size=1) +
  facet_wrap(.~ratio,nrow=3,scales="fixed") +
  theme_tufte(base_family = "sans",base_size=15) +
  theme(axis.line = element_line(size=.25))


# mu_obs vs mu_rep scatter
pred %>% 
  group_by(draw,id) %>% 
  summarise(mu_rep=mean(yrep),
            mu_obs=mean(x)) %>%
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
