library(tidyverse)
library(reticulate)

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

use_python("/Users/adkinsty/opt/miniconda3/envs/reach_env/bin/python3.8")
source_python("preproc/help.py")
library("future.apply")

get_id <- function(x) {return(substr(x,15,17))}

# TRAIN DATA
train_files <- Sys.glob(paths = "data/raw/exp1/*[1-9]_reach_train1_output.csv")
train_ids <- sapply(train_files, get_id, USE.NAMES = F)
raw_train <- tibble()
for (i in seq_along(train_ids)) {
    raw_train <- read_csv(train_files[i]) %>%
        mutate(id = rep(train_ids[i],n())) %>%
        mutate(sub_id = rep(i,n())) %>%
        bind_rows(raw_train)
}
clean_train <- raw_train %>%
  filter(pokeTime < 1) %>%
  mutate(radius = 64/2,
         reach_x = unlist(sapply(pokePos, get_x_pos, simplify = T, USE.NAMES=F)),
         reach_y = unlist(sapply(pokePos, get_y_pos, simplify = T, USE.NAMES=F)),
         gain_x = unlist(sapply(aimDotPos, get_x_pos, simplify = T, USE.NAMES=F)),
         gain_y = unlist(sapply(aimDotPos, get_y_pos, simplify = T, USE.NAMES=F)),
         std_reach_x = (reach_x - gain_x) / radius,
         std_reach_y = (reach_y - gain_y) / radius) %>%
  group_by(id) %>%
  mutate(bias_x = mean(std_reach_x,na.rm=T),
         bias_y = mean(std_reach_y,na.rm=T)) %>% ungroup() %>%
  mutate(x = std_reach_x - bias_x,
         y = std_reach_y - bias_y,
         rt = pokeTime) %>%
  dplyr::select(c(x,y,rt,radius,id,sub_id))
write_csv(clean_train, "data/clean/exp1/train_data.csv")


# TEST DATA
test_files <- Sys.glob(paths = "data/raw/exp1/*[1-9]_reach_test_output.csv")
test_ids <- sapply(test_files, get_id, USE.NAMES = F)
raw_test <- tibble()
for (i in seq_along(test_ids)) {
    raw_test <- read_csv(test_files[i]) %>%
        mutate(id = rep(test_ids[i],n())) %>%
        mutate(sub_id = rep(i,n())) %>%
        bind_rows(raw_test)
}
clean_test <- raw_test %>%
  filter(pokeTime < 1) %>%
  mutate(radius = targSize/2,
         reach_x = unlist(sapply(pokePos, get_x_pos, simplify = T, USE.NAMES=F)),
         reach_y = unlist(sapply(pokePos, get_y_pos, simplify = T, USE.NAMES=F)),
         gain_x = unlist(sapply(targetPos, get_x_pos, simplify = T, USE.NAMES=F)),
         gain_y = unlist(sapply(targetPos, get_y_pos, simplify = T, USE.NAMES=F)),
         loss_x = unlist(sapply(penaltyPos, get_x_pos, simplify = T, USE.NAMES=F)),
         loss_y = unlist(sapply(penaltyPos, get_y_pos, simplify = T, USE.NAMES=F)),
         std_reach_x = (reach_x - gain_x) / radius,
         std_reach_y = (reach_y - gain_y) / radius) %>%
  group_by(id) %>%
  mutate(bias_x = mean(std_reach_x,na.rm=T),
         bias_y = mean(std_reach_y,na.rm=T)) %>% ungroup() %>%
  mutate(reflect_x = ifelse(loss_x > gain_x, -std_reach_x, std_reach_x),
         x = reflect_x - bias_x,
         y = std_reach_y - bias_y,
         rt = pokeTime,
         loss = penaltyVal,
         gain = rewardVal,
         ratio = loss / gain,
         cond_id = ifelse(
           loss==-15,1,ifelse(loss==-5,2,ifelse(loss==-3,3,ifelse(loss==-1,4,5)))),
         dist = abs(sepDist / radius),
         dist_id = ifelse(dist==1,1,0),
         cond_id = ifelse(dist_id==1,cond_id,cond_id+5)) %>%
  dplyr::select(c(x,y,rt,gain,loss,ratio,dist,radius,id,block,sub_id,cond_id,dist_id))
write_csv(clean_test, "data/clean/exp1/test_data.csv")
