# Title     : maximum expected gain model
# Objective : Interface with stan model
# Created by: adkinsty
# Updated on: 3/11/21

library(tidyverse)

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp2/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

fm <- read_rds("data/clean/exp2/forward_model.rds")

# key variables
jj <- data$sub_id
kk <- data$cond_id
J <- length(unique(jj))
K <- length(unique(kk))
cond <- sort(unique(data$loss))
y <- data$x

# Grid optimization to obtain predicted aim points (means)
hyp <- seq(0,1,.01) 
aim <- matrix(nrow=3,ncol=J)
for (j in 1:J) {
  
  p = fm[,,j]
  
  pg <- p[1,]; pgl <- p[2,]; pl <- p[3,]; pn <- p[4,]
  
  for (k in 1:K) {
    g = 1
    l = cond[k]
    gl = g + l
    EU = pg*g + pgl*gl + pl*l+ pn*0
    aim[k,j] = hyp[which.max(EU)]
  }
}
write_rds(aim,"data/clean/exp2/MEG_aims.rds")
