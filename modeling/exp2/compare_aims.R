# Compares models using predicted vs observed aimpoint covariance

library(tidyverse)
library(ggthemes)
library(bayesplot)

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2) %>%
  mutate(lr = ifelse(ratio==0,0,ifelse(ratio==-1,1,2)))

y <- data$x # observation
ymm <- y * 8.55 # mm per radius

# model aimpoints
meg_aim <- read_rds("data/clean/exp2/MEG_aims.rds")
mleg_aim <- read_rds("data/clean/exp2/MLEG_aims.rds")
mswu_aim <- read_rds("data/clean/exp2/MSWU_aims.rds")
mswu_g_aim <- read_rds("data/clean/exp2/MSWU_gamma_aims.rds")
meu_aim <- read_rds("data/clean/exp2/MEU_aims.rds")
h_aim <- read_rds("data/clean/exp2/H_aims.rds")

aims_array <- array(dim=c(3,15,4))
aims_array[,,1] <- meu_aim
aims_array[,,2] <- h_aim
aims_array[,,3] <- meg_aim
aims_array[,,4] <- mleg_aim
# aims_array[,,5] <- mswu_aim
# aims_array[,,6] <- mswu_g_aim


J <- ncol(meg_aim) # num subjects
K <- 4 # num models
aims <- tibble()
for (k in 1:K) {
  tmp <- as_tibble(t(aims_array[,,k]))
  colnames(tmp) <- 1:3
  aims <- tmp %>%
    mutate(sub_id =  1:n(),
           model = rep(k,n())) %>%
    pivot_longer(cols=1:3,
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
         rat = ifelse(cond_id==3,"zero",ifelse(cond_id==2,"small","large")),
         sep = ifelse(cond_id > 5, "1.4R","1R"),
         m=ifelse(model==1,"MEU",
                      ifelse(model==2,"H",
                             ifelse(model==3,"MEV","MELV"))),
         m=factor(m,levels=c("MEU","H","MEV","MELV"))) %>% ungroup()

obs_pred %>%
  ggplot(aes(x=pred,y=obs,group=sub_id,colour=m)) +
  geom_point(alpha=.5) +
  abline_01(linetype="dashed",alpha=.5) +
  coord_fixed(xlim=c(-.55,8.55),ylim=c(-.55,8.55)) +
  facet_wrap(m~.,nrow=2) +
  scale_y_continuous("Observed aimpoint (mm)",breaks=seq(0,8,2)) +
  scale_x_continuous("Predicted aimpoint (mm)",breaks=seq(0,8,2)) +
  scale_colour_manual("Model",
                      values=c("#0072b2","#e69f00","#009e73","#cc79a7"),
                      labels=c("MEU", "H","MEV","MELV")) + 
  theme_tufte(base_size = 12, base_family = "sans") +
  theme(legend.position = "none",
        legend.key.size = unit(.15,"in"),
        axis.line = element_line(size=.25))
ggsave("visuals/raw/exp2/exp2_obs_pred.pdf",units = "in",height = 3,width = 3)
