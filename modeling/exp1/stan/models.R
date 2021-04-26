# Runs bayesian models

library(rstan)
library(tidyverse)

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

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
est_H <- sampling(object=model,data=input,chains=4)
write_rds(est_H,path="modeling/exp1/stan/rds/H.rds")

input$aim <- read_rds("data/clean/exp1/MEG_aims.rds")
est_MEG <- sampling(model,data=input)
write_rds(est_MEG,path="modeling/exp1/stan/rds/MEG.rds")

input$aim <- read_rds("data/clean/exp1/MEG_vp_aims.rds")
est_MEG_vp <-sampling(model,data=input)
write_rds(est_MEG_vp,path="modeling/exp1/stan/rds/MEG_vp.rds")

input$aim <- read_rds("data/clean/exp1/MLEG_aims.rds")
est_MLEG <- sampling(model,data=input)
write_rds(est_MLEG,path="modeling/exp1/stan/rds/MLEG.rds")

input$aim <- read_rds("data/clean/exp1/MEU_aims.rds")
est_MEU <-sampling(model,data=input)
write_rds(est_MEU,path="modeling/exp1/stan/rds/MEU.rds")

input$aim <- read_rds("data/clean/exp1/UH_aims.rds")
est_UH <-sampling(model,data=input)
write_rds(est_UH,path="modeling/exp1/stan/rds/UH.rds")

input$aim <- read_rds("data/clean/exp1/LH_aims.rds")
est_LH <-sampling(model,data=input)
write_rds(est_LH,path="modeling/exp1/stan/rds/LH.rds")

input$aim <- read_rds("data/clean/exp1/MSWU_aims.rds")
est_MSWU <- sampling(model,data=input)
write_rds(est_MSWU,path="modeling/exp1/stan/rds/MSWU.rds")

# input$aim <- read_rds("data/clean/exp1/MSWU_gamma_aims.rds")
# est_MSWUg <-sampling(model,data=input)
# write_rds(est_MSWUg,path="modeling/exp1/stan/rds/MSWU_gamma.rds")