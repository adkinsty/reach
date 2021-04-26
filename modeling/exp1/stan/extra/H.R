# Title     : Maximum Expected Win Model
# Objective : Interface with stan model
# Created by: adkinsty
# Created on: 6/30/20

library(rstan)
library(tidyverse)
library(shinystan)
library(ggthemes)
library(shotGroups)
library(bayesplot)

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp1/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

data <- read_csv("data/clean/exp1/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

y <- data$x

# key variables
jj <- data$sub_id
cc <- data$cond_id
dd <- data$dist_id
con <- sort(unique(cc))
dis <- unique(dd)
C <- length(con)
D <- length(dis)
J <- length(unique(jj))
N <- nrow(data)

aim <- matrix(nrow=C,ncol=J)
for (j in 1:J) {
  # calculate win probabilities and get aim with max
  for (c in 1:C) {
    if (c %in% c(5,10)) {
      aim[c,j] = 0
    } else if (c %in% 1:4){
      aim[c,j] = .5
    } else {
      aim[c,j] = (1 - 1.375) + (1.375/2)
    }
  }
}
write.csv(aim,"data/clean/exp1/H_aims.csv")


input <- list(N=N,J=J,C=C,
              cc=cc,jj=jj,
              aim=aim,
              y=y,
              get_yrep=1,
              get_log_lik=1,
              prior_only=0)

model <- stan_model(file="modeling/exp1/stan/aim.stan",
                    model_name="H")

# mcmc_est <- readRDS("modeling/exp1/stan/rds/H.rds")
mcmc_est <- sampling(
  object=model,
  data=input,
  chains=4,
  iter=5000,
  warmup=2500,
  cores=4)
write_rds(mcmc_est,path="modeling/exp1/stan/rds/H.rds")