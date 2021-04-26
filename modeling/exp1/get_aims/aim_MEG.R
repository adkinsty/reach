# Title     : maximum expected gain model
# Objective : Interface with stan model
# Created by: adkinsty
# Updated on: 3/11/21

library(tidyverse)
library(optimx)

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp1/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

data <- read_csv("data/clean/exp1/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

fm <- read_rds("data/clean/exp1/forward_model.rds")

# key variables
jj <- data$sub_id
kk <- data$cond_id
J <- length(unique(jj))
subs <- sort(unique(data$id))
K <- length(unique(kk))

# Grid optimization to obtain predicted aim points (means)
hyp <- seq(0,1,.01) 
aim <- matrix(nrow=K,ncol=J)
for (j in 1:J) {
  
  print(sprintf("Subject %s",j))

  tmp <- data %>% filter(sub_id == j)
  
  all_p = fm[,,j,]
    
  # calculate win probabilities and get aim with max
  for (k in 1:K) {
    
    if (k < 6) {
      p <- all_p[,,1]
    } else {
      p <- all_p[,,2]
    }
    
    # get outcome probs based on separation distance
    pg <- p[1,]; pgl <- p[2,];  pl <- p[3,]; pn <- p[4,]
    
    # get loss value based on cond id
    if (k %in% c(1,6)) {
      l <- -15
    } else if (k %in% c(2, 7)) {
      l <- -5
    } else if (k %in% c(3, 8)) {
      l <- -3
    } else if (k %in% c(4, 9)) {
      l < -1
    } else {
      l <- 0
    }
    
    # get gain value based on loss
    if (l %in% c(-15, -3)) {
      g <- 3
    } else {
      g <- 1
    }
    
    gl = g + l
    
    EU = pg*g + pgl*gl + pl*l+ pn*0
    aim[k,j] = hyp[which.max(EU)]
    
  }
}
write_rds(aim,"data/clean/exp1/MEG_aims.rds")
