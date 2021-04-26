library(tidyverse)
library(reticulate)

use_python("/Users/adkinsty/opt/miniconda3/envs/reach_env/bin/python3.8")
source_python("preproc/help.py")

get_id <- function(x) {return(substr(x,15,17))}

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

# TRAIN DATA
train_files <- Sys.glob(paths = "data/raw/exp2/*[1-9]_reach_train_output.csv")
train_ids <- sapply(train_files, get_id, USE.NAMES = F)
raw_train <- tibble()
for (i in seq_along(train_ids)) {
    raw_train <- read_csv(train_files[i]) %>%
        mutate(id = rep(train_ids[i],n())) %>%
        mutate(sub_id = rep(i,n())) %>%
        bind_rows(raw_train)
}
clean_train <- raw_train %>%
  filter(!too_slow & !too_soon) %>%
  mutate(radius = 64/2,
         reach_x = poke_x,
         reach_y = poke_y,
         gain_x = targ_x,
         gain_y = targ_y,
         loss_x = pen_x,
         loss_y = targ_y,
         std_reach_x = (reach_x - gain_x) / radius,
         std_reach_y = (reach_y - gain_y) / radius) %>%
  group_by(id) %>%
  mutate(bias_x = mean(std_reach_x,na.rm=T),
         bias_y = mean(std_reach_y,na.rm=T)) %>% ungroup() %>%
  mutate(reflect_x = ifelse(loss_x > gain_x, -std_reach_x, std_reach_x),
         x = reflect_x - bias_x,
         y = std_reach_y - bias_y,
         rt = trial_time,
         loss = pen_val/100,
         gain = rep(1,n()),
         ratio = loss / gain,
         cond_id = ifelse(loss==-5,1,ifelse(loss==-1,2,3)),
         mult = ifelse(gain == 3, 1, 0),
         dist = radius/radius) %>%
  select(c(x,y,rt,gain,loss,ratio,cond_id,mult,dist,radius,id,sub_id,block))
write_csv(clean_train, "data/clean/exp2/train_data.csv")


# TEST DATA
test_files <- Sys.glob(paths = "data/raw/exp2/*[1-9]_reach_test_output.csv")
test_ids <- sapply(test_files, get_id, USE.NAMES = F)
raw_test <- tibble()
for (i in seq_along(test_ids)) {
    raw_test <- read_csv(test_files[i]) %>%
        mutate(id = rep(test_ids[i],n())) %>%
      mutate(sub_id = rep(i,n())) %>%
        bind_rows(raw_test)
}
clean_test <- raw_test %>%
  filter(!too_slow & !too_soon) %>%
  mutate(radius = 64/2,
         reach_x = poke_x,
         reach_y = poke_y,
         gain_x = targ_x,
         gain_y = targ_y,
         loss_x = pen_x,
         loss_y = targ_y,
         std_reach_x = (reach_x - gain_x) / radius,
         std_reach_y = (reach_y - gain_y) / radius) %>%
  group_by(id) %>%
  mutate(bias_x = mean(std_reach_x,na.rm=T),
         bias_y = mean(std_reach_y,na.rm=T)) %>% ungroup() %>%
  mutate(reflect_x = ifelse(loss_x > gain_x, -std_reach_x, std_reach_x),
         x = reflect_x - bias_x,
         y = std_reach_y - bias_y,
         rt = trial_time,
         loss = pen_val/100,
         gain = rep(1,n()),
         ratio = loss / gain,
         cond_id = ifelse(loss==-5,1,ifelse(loss==-1,2,3)),
         mult = ifelse(gain == 3, 1, 0),
         dist = radius/radius) %>%
  dplyr::select(c(x,y,rt,gain,loss,ratio,cond_id,mult,dist,radius,id,sub_id,block))
write_csv(clean_test, "data/clean/exp2/test_data.csv")