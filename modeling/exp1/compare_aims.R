# Script for comparing the models in terms of predicted vs observed aimpoints
# Author: Tyler Adkins
# Date: 9/10/20

library(tidyverse)
library(ggthemes)

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

data <- read_csv("data/clean/exp1/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2) %>%
  mutate(lr = ifelse(ratio==0,0,ifelse(ratio==-1,1,2)))

y <- data$x # observation
ymm <- y * 8.55 # mm per radius

data %>% 
  mutate(ratio=factor(loss,levels=c(0,-1,-3,-5,-15)),
         x = x*8.55,
         dist=ifelse(dist==1,"Near","Far")) %>%
  ggplot(aes(x=x,y=ratio)) +
  facet_wrap(dist~.,nrow=2) +
  geom_point(size=.1,alpha=.5,position = position_jitter(),colour="black") +
  stat_summary() + 
  scale_x_continuous("Observed endpoint (mm)",breaks=seq(-10,10,2)) +
  scale_y_discrete("Loss ratio",labels=c("0:1","1:1","3:3","5:1","15:3")) +
  geom_vline(xintercept = 0, colour="black",linetype="dashed") +
  coord_cartesian(xlim=c(0,8.55)) +
  theme_tufte(base_size=10,base_family="sans") +
  theme(axis.line=element_line(size=.25),
        legend.pos="none")
ggsave("visuals/raw/exp1/exp1_obs.pdf",units="in",height=2.25,width=4)

# model aimpoints
meg_aim <- read_rds("data/clean/exp1/MEG_aims.rds")
mleg_aim <- read_rds("data/clean/exp1/MLEG_aims.rds")
meu_aim <- read_rds("data/clean/exp1/MEU_aims.rds")
h_aim <- read_rds("data/clean/exp1/H_aims.rds")

aims_array <- array(dim=c(10,22,4))
aims_array[,,1] <- meu_aim
aims_array[,,2] <- h_aim
aims_array[,,3] <- meg_aim
aims_array[,,4] <- mleg_aim
J <- ncol(meg_aim) # num subjects
K <- 4 # num models
aims <- tibble()
for (k in 1:K) {
  tmp <- as.data.frame(t(aims_array[,,k]))
  colnames(tmp) <- 1:10
  aims <- tmp %>%
    mutate(sub_id =  1:n(),
           model = rep(k,n())) %>%
    pivot_longer(cols=1:10,
                 names_to="cond_id",
                 values_to="x") %>%
    mutate(cond_id = as.integer(cond_id)) %>%
    rbind(aims)
}

obs_pred <- data %>% group_by(sub_id,cond_id) %>%
  summarise(obs = mean(x)) %>% ungroup() %>%
  merge(aims) %>%
  mutate(pred=x*8.55,
         obs=obs*8.55,
         rat = ifelse(cond_id==10 | cond_id==5,"zero",
                      ifelse(cond_id == 1 |
                             cond_id == 2 |
                             cond_id == 6 |
                             cond_id == 7, "large","small")),
         sep = ifelse(cond_id > 5, "1.4R","1R"),
         model=ifelse(model==1,"MEU",
                      ifelse(model==2,"H",
                              ifelse(model==3,"MEV","MELV"))),
         model=factor(model,levels=c("MEU","H","MEV","MELV"))) %>% ungroup()

obs_pred %>%
  # filter(rat!="zero") %>%
  ggplot(aes(x=pred,y=obs,group=sub_id,colour=model)) +
  geom_point(alpha=.5) +
  abline_01(linetype="dashed",alpha=.5) +
  coord_fixed(xlim=c(-.55,8.55),ylim=c(-.55,8.55)) +
  facet_wrap(model~.,nrow=2) +
  scale_y_continuous("Observed aimpoint (mm)",breaks=seq(0,8,2)) +
  scale_x_continuous("Predicted aimpoint (mm)",breaks=seq(0,8,2)) +
  scale_shape(guide=F) +
  scale_colour_manual("Model",
                      values=c("#0072b2","#e69f00","#009e73","#cc79a7"),
                      labels=c("MEU", "H","MEV","MELV")) + 
  theme_tufte(base_size = 12, base_family = "sans") +
  theme(legend.position = "none",
        legend.key.size = unit(.15,"in"),
        axis.line = element_line(size=.25))
ggsave("visuals/raw/exp1/exp1_obs_pred.pdf",units = "in",height = 3,width = 3)






