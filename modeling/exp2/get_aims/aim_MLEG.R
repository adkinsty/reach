# Title     : maximum loss-averse expected gain model
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

la_df <- read_tsv("data/raw/loss_aversion/All_lambdas.csv") %>%
  filter(subNum %in% data$id) %>% arrange(subNum) %>% 
  mutate(id = subNum, la = lambdas)

# key variables
jj <- data$sub_id
kk <- data$cond_id
J <- length(unique(jj))
K <- length(unique(kk))
cond <- sort(unique(data$loss))
subs <- sort(unique(data$id))
y <- data$x

# Grid optimization to obtain predicted aim points (means)
hyp <- seq(0,1,.01) 
aim <- matrix(nrow=3,ncol=J)
for (j in 1:J) {
  
  p = fm[,,j]
  
  pg <- p[1,]; pgl <- p[2,]; pl <- p[3,]; pn <- p[4,]

  raw_la <- abs(as.numeric(la_df[la_df$'id'==subs[j],'la']))
  lambda <- ifelse(is.na(raw_la), 1, raw_la)
  
  for (k in 1:K) {
    g = 1
    l = cond[k]
    uloss <- l * lambda
    gl = g + uloss
    EU = pg*g + pgl*gl + pl*uloss+ pn*0
    aim[k,j] = hyp[which.max(EU)]
  }
}
write_rds(aim,"data/clean/exp2/MLEG_aims.rds")
