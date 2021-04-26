# Runs bayesian models

library(rstan)
library(tidyverse)

#options(mc.cores = parallel::detectCores())
#parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
options(mc.cores = 1)

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp1/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

data <- read_csv("data/clean/exp1/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

# key variables
jj <- data$sub_id
cc <- data$cond_id
J <- length(unique(jj))
N <- nrow(data)
C <- length(unique(cc))
cond <- sort(unique(data$loss))
y <- data$x

input <- list(N=N,J=J,C=C,cc=cc,jj=jj,aim=NA,y=y,get_gq=1,prior_only=0)

model <- stan_model(file="modeling/exp1/stan/stan_model.stan",model_name="model")

input$aim <- read_rds("data/clean/exp1/H_aims.rds") 
est <- sampling(object=model,data=input,chains=4)
write_rds(est,path="/Users/adkinsty/files/offline_reach/exp1/rds/H.rds")

input$aim <- read_rds("data/clean/exp1/MEG_aims.rds")
est <- sampling(model,data=input)
write_rds(est,path="/Users/adkinsty/files/offline_reach/exp1/stan/rds_revision/MEG.rds")

input$aim <- read_rds("data/clean/exp1/MLEG_aims.rds")
est <- sampling(model,data=input)
write_rds(est,path="/Users/adkinsty/files/offline_reach/exp1/stan/rds_revision/MLEG.rds")

input$aim <- read_rds("data/clean/exp1/MEU_aims.rds")
est <-sampling(model,data=input)
write_rds(est,path="/Users/adkinsty/files/offline_reach/exp1/stan/rds_revision/MEU.rds")

input$aim <- read_rds("data/clean/exp1/UH_aims.rds")
est <-sampling(model,data=input)
write_rds(est,path="/Users/adkinsty/files/offline_reach/exp1/stan/rds_revision/UH.rds")

input$aim <- read_rds("data/clean/exp1/MSWU_aims.rds")
est <- sampling(model,data=input)
write_rds(est,path="/Users/adkinsty/files/reach_offline/exp1/stan/rds_revision/MSWU.rds")

# input$aim <- read_rds("data/clean/exp1/MEG_vp_aims.rds")
# est <-sampling(model,data=input)
# write_rds(est,path="modeling/exp1/stan/rds_revision/MEG_vp.rds")

# input$aim <- read_rds("data/clean/exp1/MSWU_gamma_aims.rds")
# est <-sampling(model,data=input)
# write_rds(est,path="modeling/exp1/stan/rds_revision/MSWU_gamma.rds")