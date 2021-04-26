# Title     : maximum expected gain model with variable precision 
# Objective : Interface with stan model
# Created by: adkinsty
# Updated on: 6/29/20

library(tidyverse)

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp2/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

fm <- read_rds("data/clean/exp2/forward_model_vp.rds")

# key variables
jj <- data$sub_id
kk <- data$cond_id
J <- length(unique(jj))
N <- nrow(data)
K <- length(unique(kk))
cond <- sort(unique(data$loss))
y <- data$x

# Grid optimization to obtain predicted aim points (means)
hyp <- seq(0,1,.01) 
G <- length(hyp)
aim <- matrix(nrow = 3,ncol = J)
for (j in 1:J) {
  
  for (k in 1:K) {
    
    p = fm[,,j,k]
    
    pg <- p[1,]; pgl <- p[2,]; pl <- p[3,]; pn <- p[4,]
    
    g = 1
    l = cond[k]
    gl = g + l
    EU = pg*g + pgl*gl + pl*l + pn*0
    aim[k,j] = hyp[which.max(EU)]
  }
}
write_rds(aim,"data/clean/exp2/MEG_vp_aims.rds")
