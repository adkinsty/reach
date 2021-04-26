# Script for comparing the models in terms of predicted vs observed aimpoints
# Author: Tyler Adkins
# Date: 9/10/20

library(tidyverse)
library(ggthemes)
library(bayesplot)

wong <- c("#000000", "#e69f00", "#56b4e9", "#009e73", "#f0e442", "#0072b2", "#d55e00", "#cc79a7")

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
  geom_point(size=.1,alpha=.1,position = position_jitter(width = 0),colour="black") +
  stat_summary() + 
  scale_x_continuous("Observed endpoint (mm)",breaks=seq(-10,10,2)) +
  scale_y_discrete("Loss:Gain",labels=c("0:1","1:1","3:3","5:1","15:3")) +
  geom_vline(xintercept = 0, colour="black",linetype="dashed") +
  coord_cartesian(xlim=c(0,8.55)) +
  theme_tufte(base_size=10,base_family="sans") +
  theme(axis.line=element_line(size=.25),
        legend.pos="none")
ggsave("visuals/raw/exp1/exp1_obs.pdf",units="in",height=3,width=4)

# model aimpoints
meg_aim <- read_rds("data/clean/exp1/MEG_aims.rds")
mleg_aim <- read_rds("data/clean/exp1/MLEG_aims.rds")
meu_aim <- read_rds("data/clean/exp1/MEU_aims.rds")
mswu_aim <- read_rds("data/clean/exp1/MSWU_aims.rds")
uh_aim <- read_rds("data/clean/exp1/UH_aims.rds")
ulh_aim <- read_rds("data/clean/exp1/ULH_aims.rds")
h_aim <- read_rds("data/clean/exp1/H_aims.rds")
meg_vp_aim <- read_rds("data/clean/exp1/MEG_vp_aims.rds")

K <- 6 # num models
models <- c("MEG", "MLEG", "MEU", "MSWU", "UH", "H")

aims_array <- array(dim=c(10,22,K))
aims_array[,,1] <- meg_aim
aims_array[,,2] <- mleg_aim
aims_array[,,3] <- meu_aim
aims_array[,,4] <- mswu_aim
aims_array[,,5] <- uh_aim
aims_array[,,6] <- h_aim

J <- ncol(meg_aim) # num subjects

aims <- tibble()
for (k in 1:K) {
  tmp <- as.data.frame(t(aims_array[,,k]))
  colnames(tmp) <- 1:10
  aims <- tmp %>%
    mutate(sub_id =  1:n(),
           m = rep(models[k],n())) %>%
    pivot_longer(cols=1:10,
                 names_to="cond_id",
                 values_to="x") %>%
    mutate(cond_id = as.integer(cond_id)) %>%
    rbind(aims)
}

models_plt <- c("MLEG", "MEG", "MSWU", "UH")
obs_pred <- data %>% group_by(sub_id,cond_id) %>%
  summarise(obs = mean(x)) %>% ungroup() %>%
  merge(aims) %>%
  filter(m %in% models_plt) %>%
  mutate(pred=x*8.55,
         obs=obs*8.55,
         rat = ifelse(cond_id==10 | cond_id==5,"zero",
                      ifelse(cond_id == 1 |
                             cond_id == 2 |
                             cond_id == 6 |
                             cond_id == 7, "large","small")),
         sep = ifelse(cond_id > 5, "1.4R","1R"),
         model=factor(m,levels=models)) %>% ungroup()

labels <- c("LA", "MEG", "NLVP", "UH")

obs_pred %>%
  # filter(rat!="zero") %>%
  ggplot(aes(x=pred,y=obs,group=sub_id,colour=model)) +
  geom_point(alpha=.25) +
  abline_01(linetype="dashed",alpha=.5) +
  coord_fixed(xlim=c(-.55,8.55),ylim=c(-.55,8.55)) +
  facet_wrap(model~.,nrow=2, labeller = labeller(model = c("MEG"="MEG","MLEG"="LA", "MSWU"="NLVP","UH"="UH"))) +
  scale_y_continuous("Observed aimpoint (mm)",breaks=seq(0,8,2)) +
  scale_x_continuous("Predicted aimpoint (mm)",breaks=seq(0,8,2)) +
  scale_shape(guide=F) +
  scale_colour_manual("Model", values=wong[c(4,2,3,1)]) + 
  theme_tufte(base_size = 10, base_family = "sans") +
  theme(legend.position = "none",
        legend.key.size = unit(.15,"in"),
        axis.line = element_line(size=.25))
ggsave("visuals/raw/exp1/revision/exp1_obs_pred.pdf",units = "in",height = 4,width = 4)






