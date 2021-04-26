# Title     : Get aimpoints for heuristic model
# Objective : inputs to stan interface
# Created by: adkinsty
# Created on: 9/10/20

library(tidyverse)

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
kk <- data$cond_id
K <- length(unique(kk))
J <- length(unique(jj))

aim <- matrix(nrow=K,ncol=J)
for (j in 1:J) {
  # calculate win probabilities and get aim with max
  for (k in 1:K) {
    if (k %in% c(5,10)) {
      aim[k,j] = 0
    } else if (k %in% 1:4){
      aim[k,j] = .5
    } else {
      aim[k,j] = (1 - 1.375) + (1.375/2)
    }
  }
}
write_rds(aim,"data/clean/exp1/H_aims.rds")
