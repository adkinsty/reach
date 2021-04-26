# Title     : Hierarchical Univariate Normal Maximum Expected Utility Model (non-centered param)
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
  x=as.matrix(cbind(rep(1,nrow(data)),data$ratio)),  # trial predictors
  u=as.matrix(cbind(rep(1,nrow(la_df)))), # subject predictors
  y=data$x,                                          # observations
  get_yrep=1,                     
  get_log_lik=0,
  prior_only=1
)

model <- stan_model(
  file="modeling/stan/scratch_model4.stan",
  model_name="scratch_model4"
)

fit <- sampling(
  object=model,
  data=input,
  chains=3,
  iter=1100,
  warmup=1000,
  cores=3,
  control=list(max_treedepth=10))

y <- scale(data$x)[,1]
yrep <- rstan::extract(fit, 'yrep')$yrep

# launch_shinystan(fit)

pred <- as.tibble(t(yrep)) %>% cbind(data) %>%
  pivot_longer(cols=starts_with("V"),
               names_to="draw",
               values_to="yrep")

pred %>% 
  group_by(draw) %>% 
  summarise(mu_rep=mean(yrep)) %>%
  ggplot() +
  # coord_cartesian(xlim=c(.36,.41)) +
  geom_histogram(aes(x=mu_rep),binwidth=.003,colour="white",fill="grey") +
  geom_density(aes(x=mu_rep),colour="grey",size=1) +
  geom_vline(xintercept=mean(data$x),size=1,linetype="solid",colour="green3") + 
  theme_tufte(base_family = "sans",base_size=15) +
  theme(axis.line = element_line(size=.25))


pred %>% 
  group_by(draw) %>% 
  summarise(sd_rep=sd(yrep)) %>%
  ggplot() +
  geom_histogram(aes(x=sd_rep),binwidth = .003,colour="white",fill="grey") +
  geom_density(aes(x=sd_rep),colour="grey",size=1) +
  geom_vline(xintercept=sd(data$x),size=1,linetype="solid",colour="green3") + 
  theme_tufte(base_family = "sans",base_size=15) +
  theme(axis.line = element_line(size=.25))


pred %>% 
  group_by(draw,ratio) %>% 
  summarise(mu_rep=mean(yrep)) %>%
  ggplot() +
  geom_histogram(aes(x=mu_rep),binwidth = .01,colour="white",fill="grey") +
  geom_vline(data=data%>%group_by(ratio)%>%summarise(mu=mean(x)),
             aes(xintercept=mu),linetype="dotted",size=1) +
  facet_wrap(.~ratio,nrow=1,scales="fixed") +
  theme_tufte(base_family = "sans",base_size=15) +
  theme(axis.line = element_line(size=.25))


pred %>% 
  group_by(draw,id) %>% 
  summarise(mu_rep=mean(yrep),
            mu_obs=mean(x)) %>%
  ggplot() +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
  geom_abline(slope=1,intercept=0,colour="grey",linetype="dashed") +
  geom_point(aes(x=mu_rep,y=mu_obs),alpha=.008) +
  theme_tufte(base_family = "sans",base_size=15) +
  theme(axis.line = element_line(size=.25))

pred %>% 
  group_by(draw,id) %>% 
  summarise(sd_rep=sd(yrep),
            sd_obs=sd(x)) %>%
  ggplot() +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
  geom_abline(slope=1,intercept=0,colour="grey",linetype="dashed") +
  geom_point(aes(x=sd_rep,y=sd_obs),alpha=.008) +
  theme_tufte(base_family = "sans",base_size=15) +
  theme(axis.line = element_line(size=.25))

pred %>% 
  group_by(draw,id,ratio) %>% 
  summarise(mu_rep=mean(yrep),
            mu_obs=mean(x)) %>%
  ggplot() +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
  geom_abline(slope=1,intercept=0,colour="grey",linetype="dashed") +
  geom_point(aes(x=mu_rep,y=mu_obs,colour=factor(ratio)),alpha=.01) +
  scale_colour_colorblind() +
  theme_tufte(base_family = "sans",base_size=15) +
  theme(axis.line = element_line(size=.25),
        legend.position = "none")



pred %>% 
  group_by(draw,id,ratio) %>% 
  summarise(mu_rep=mean(yrep),
            mu_obs=mean(x))