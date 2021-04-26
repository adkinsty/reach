# Compares models using predicted vs observed aimpoint covariance

library(tidyverse)
library(ggthemes)
library(bayesplot)

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

wong <- c("#000000", "#e69f00", "#56b4e9", "#009e73", "#f0e442", "#0072b2", "#d55e00", "#cc79a7")

data <- read_csv("data/clean/exp2/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2) %>%
  mutate(lr = ifelse(ratio==0,0,ifelse(ratio==-1,1,2)))

y <- data$x # observation
ymm <- y * 8.55 # mm per radius

# model aimpoints
meg_aim <- read_rds("data/clean/exp2/MEG_aims.rds")
mleg_aim <- read_rds("data/clean/exp2/MLEG_aims.rds")
meu_aim <- read_rds("data/clean/exp2/MEU_aims.rds")
mswu_aim <- read_rds("data/clean/exp2/MSWU_aims.rds")
uh_aim <- read_rds("data/clean/exp2/UH_aims.rds")
h_aim <- read_rds("data/clean/exp2/H_aims.rds")
#ulh_aim <- read_rds("data/clean/exp2/ULH_aims.rds")
#meg_vp_aim <- read_rds("data/clean/exp2/MEG_vp_aims.rds")

K <- 6 # num models
models <- c("MEG", "MLEG", "MEU", "MSWU", "UH", "H")

aims_array <- array(dim=c(3,15,K))
aims_array[,,1] <- meg_aim
aims_array[,,2] <- mleg_aim
aims_array[,,3] <- meu_aim
aims_array[,,4] <- mswu_aim
aims_array[,,5] <- uh_aim
aims_array[,,6] <- h_aim

J <- ncol(meg_aim) # num subjects
aims <- tibble()
for (k in 1:K) {
  tmp <- as_tibble(t(aims_array[,,k]))
  colnames(tmp) <- 1:3
  aims <- tmp %>%
    mutate(sub_id =  1:n(),
           m = models[k],
           k = rep(k,n())) %>%
    pivot_longer(cols=1:3,
                 names_to="cond_id",
                 values_to="x") %>%
    mutate(cond_id = as.integer(cond_id)) %>%
    rbind(aims)
}

models_plt <- c("MLEG", "MEG", "MSWU", "UH")
labels <- c("LA","MEG","NLVP","UH")

obs_pred <- data %>% group_by(sub_id,cond_id) %>%
  summarise(obs = mean(x)) %>% ungroup() %>%
  merge(aims) %>%
  filter(m %in% models_plt) %>%
  mutate(pred=x*8.55,
    obs=obs*8.55,
    rat = ifelse(cond_id==3,"zero",ifelse(cond_id==2,"small","large")),
    sep = ifelse(cond_id > 5, "1.4R","1R"),
    m = factor(m,levels=models_plt)) %>% ungroup()

obs_pred %>%
  ggplot(aes(x=pred,y=obs,group=sub_id,colour=m)) +
  geom_point(alpha=.25) +
  abline_01(linetype="dashed",alpha=.5) +
  coord_fixed(xlim=c(-.55,8.55),ylim=c(-.55,8.55)) +
  facet_wrap(m~.,nrow=2, labeller = labeller(m = c("MEG"="MEG","MLEG"="LA", "MSWU"="NLVP","UH"="UH"))) +
  scale_y_continuous("Observed aimpoint (mm)",breaks=seq(0,8,2)) +
  scale_x_continuous("Predicted aimpoint (mm)",breaks=seq(0,8,2)) +
  scale_colour_manual("Model", values=wong[c(2,4,3,1)], labels=labels) + 
  theme_tufte(base_size = 10, base_family = "sans") +
  theme(legend.position = "none",
        legend.key.size = unit(.15,"in"),
        axis.line = element_line(size=.25))
ggsave("visuals/raw/exp2/revision/exp2_obs_pred.pdf",units = "in",height = 3,width = 3)
