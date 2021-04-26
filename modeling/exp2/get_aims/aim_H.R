# Title     : get heuristic aimpoints
# Objective : Interface with stan model
# Created by: adkinsty
# Created on: 6/25/20

library(tidyverse)

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp2/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

# key variables
jj <- data$sub_id
kk <- data$cond_id
J <- length(unique(jj))
K <- length(unique(kk))
cond <- sort(unique(data$loss))
y <- data$x

aim <- matrix(nrow=3,ncol=J)
for (j in 1:J) {
  # calculate win probabilities and get aim with max
  for (k in 1:K) {
    if (cond[k]==0) {
      aim[k,j] = 0
    } else {
      aim[k,j] = .5
    }
  }
}
write_rds(aim,"data/clean/exp2/H_aims.rds")
