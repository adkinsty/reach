library(brms)
library(tidyverse)

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")

train <- read_csv("data/clean/exp1/train_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

data <- read_csv("data/clean/exp1/test_data.csv") %>%
  filter(x > -2 & x < 2 & y > -2 & y < 2)

la_df <- read_tsv("data/raw/loss_aversion/All_lambdas.csv") %>%
  filter(subNum %in% data$id) %>% arrange(subNum) %>% 
  mutate(id = subNum, la = lambdas) %>%
  dplyr::select(id, la)

# key variables
jj <- data$sub_id
cc <- data$cond_id
J <- length(unique(jj))
subs <- sort(unique(data$id))
N <- nrow(data)
C <- length(unique(cc))
cond <- sort(unique(data$loss))
y <- data$x

tmp <- train %>% 
  group_by(sub_id) %>%
  summarise(sd = sd(x)) %>% 
  merge(data) %>%
  mutate(la = NA)
for (j in 1:J){
  lambda <- abs(as.numeric(la_df[la_df$id==subs[j],'la']$la))
  tmp <- tmp %>% 
    mutate(la = ifelse(sub_id != j, la, 
                       ifelse(is.na(lambda), 1.5,
                              ifelse(lambda > 5, 5,
                                     ifelse(lambda < 1, 1, lambda)))))
}
input <- tmp %>% 
  filter(ratio != 0) %>%
  mutate(la = ifelse(is.na(la),1,la),
         zx = scale(x), 
         zla = scale(rank(la)),
         zsd = scale(rank(sd)),
         rat = ifelse(ratio==-5,.5,-.5),
         mult = ifelse(cond_id == 1 | cond_id == 3 | cond_id == 6 | cond_id == 7, .5, -.5),
         dist = ifelse(dist==1,.5,-.5))

m <- brm(bf(zx ~ 1 + mult + dist + rat + zla + zsd + (1 |sub_id), 
            sigma ~ 1 + (1 | sub_id)),
         data=input,
         prior=c(prior(normal(0,1),class='b'),
                 prior(normal(0,1),class='Intercept'),
                 prior(normal(0,1),class='sd'),
                 prior(normal(1,1),class='Intercept',dpar="sigma")),
         iter=5000,
         warmup=2500,
         chains=4,
         cores=4)

input0 <- tmp %>% 
  filter(loss == 0)

m0 <- brm(bf(x ~ 1 + (1 |sub_id), 
             sigma ~ 1 + (1 | sub_id)),
          data=input,
          prior=c(prior(normal(0,1),class='Intercept'),
                  prior(normal(0,1),class='sd'),
                  prior(normal(1,1),class='Intercept',dpar="sigma")),
          iter=5000,
          warmup=2500,
          chains=4,
          cores=4)



input <- tmp %>% 
  mutate(la = ifelse(is.na(la),1,la),
         zx = scale(x), 
         zla = scale(rank(la)),
         zsd = scale(rank(sd)),
         rat = ifelse(ratio==-5,.5,ifelse(ratio==-1,0,-.5)))



mo <- brm(bf(zx ~ 1 + rat + (1 |sub_id), 
             sigma ~ 1 + (1 | sub_id)),
          data=input,
          prior=c(prior(normal(0,1),class='b'),
                  prior(normal(0,1),class='Intercept'),
                  prior(normal(0,1),class='sd'),
                  prior(normal(1,1),class='Intercept',dpar="sigma")),
          iter=5000,
          warmup=2500,
          chains=4,
          cores=4)

mv <- brm(bf(zx ~ 1 + rat + (1 |sub_id), 
             sigma ~ 1 + rat + (1 | sub_id)),
          data=input,
          prior=c(prior(normal(0,1),class='b'),
                  prior(normal(0,1),class='Intercept'),
                  prior(normal(0,1),class='sd'),
                  prior(normal(1,1),class='Intercept',dpar="sigma")),
          iter=5000,
          warmup=2500,
          chains=4,
          cores=4)

mo_yrep <- posterior_predict(mo)
mv_yrep <- posterior_predict(mv)
mo_loo <- loo(mo, cores = 6, reloo = TRUE)
mv_loo <- loo(mv, cores = 6, reloo = TRUE)

loo_compare(mo_loo,mv_loo)

y  <- input$zx
# Density
ppc_dens_overlay(y,mv_yrep[sample(nrow(mv_yrep), 25), ]) +
  coord_cartesian(ylim=c(0,1.1),xlim=c(-1,2)) +
  scale_x_continuous("Endpoint") +
  theme_tufte(base_size=15,base_family = "sans") +
  theme(legend.position = "none",
        axis.line.x = element_line(size=.25),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
ppc_loo_pit_overlay(y, h_yrep, lw = weights(h_loo$psis_object),alpha=.5) +
  coord_cartesian(xlim=c(0,1)) +
  theme_tufte(base_size=15,base_family = "sans") +
  theme(legend.position = "top",
        axis.line.x = element_line(size=.25),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())



# sd by loss ratio
data %>% group_by(sub_id,ratio) %>% 
  summarise(var = var(x*8.55)) %>%
  ggplot(aes(x=factor(ratio,levels=c(0,-1,-5)),y=var,group=sub_id)) +
  geom_line(alpha=.5,position = position_dodge(.05),colour="grey") +
  geom_point(size=.5,position = position_dodge(.05),colour="grey") +
  stat_summary(aes(group=1),fun.data = "mean_cl_boot") +
  scale_y_continuous(expression("Movement noise (mm"^"2"*")")) +
  scale_x_discrete("Loss ratio", labels=c("0:1", "1:1", "5:1")) + 
  theme_tufte(base_size=10,base_family = "sans") +
  theme(legend.position = "top",
        axis.line = element_line(size=.25))
ggsave("visuals/raw/exp1/exp1_noise_reduction.pdf",units="in",height=2,width=2)  





data %>%
  mutate(x = x*8.55) %>%
  group_by(ratio,id) %>%
  summarise(y=mean(x)) %>% ungroup() %>%
  mutate(gmu = mean(y)) %>%
  group_by(id) %>%
  mutate(smu = mean(y)) %>% ungroup() %>%
  mutate(ynew = y - smu + gmu) %>%
  group_by(ratio) %>%
  summarise(mean = mean(y),
            se_b = sd(y)/sqrt(n()),
            se_w = sd(ynew)/sqrt(n()))

data %>%
  mutate(x = x*8.55) %>%
  group_by(ratio,id) %>%
  summarise(y=sd(x)) %>% ungroup() %>%
  mutate(gmu = mean(y)) %>%
  group_by(id) %>%
  mutate(smu = mean(y)) %>% ungroup() %>%
  mutate(ynew = y - smu + gmu) %>%
  group_by(ratio) %>%
  summarise(mean = mean(y),
            se_b = sd(y)/sqrt(n()),
            se_w = sd(ynew)/sqrt(n()))
  
  
data %>%
  mutate(x = x*8.55) %>%
  group_by(dist,id) %>%
  summarise(y=mean(x)) %>% ungroup() %>%
  mutate(gmu = mean(y)) %>%
  group_by(id) %>%
  mutate(smu = mean(y)) %>% ungroup() %>%
  mutate(ynew = y - smu + gmu) %>%
  group_by(dist) %>%
  summarise(mean = mean(y),
            se_b = sd(y)/sqrt(n()),
            se_w = sd(ynew)/sqrt(n()))
