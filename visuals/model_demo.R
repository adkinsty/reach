# Title     : Max Expected Value Model Demo
# Objective : Illustrate the building blocks of the EV-maximizing model
# Created by: adkinsty
# Created on: 6/20/20


library(tidyverse)
library(ggthemes)
library(ggsci)
library(mvtnorm)
library(ggforce)
library(shotGroups)
library(functional)


setwd("/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/")


# FAKE 2D REACH DISTRIBUTION
aim = 0
reach = rmvnorm(n=2000, mean=c(aim,0), sigma=diag(2)/4)
data = tibble(x = reach[,1],
              y = reach[,2])
data %>% ggplot() +
  geom_circle(data=tibble(x=0,y=0,r=1),aes(x0=x,y0=y,r=r),
              size=1,colour="lightskyblue",fill="lightskyblue",alpha=.95) +
  geom_circle(data=tibble(x=-1,y=0,r=1),aes(x0=x,y0=y,r=r),
              size=1,colour="orange2",fill="orange2",alpha=.75) +
  stat_density_2d(aes(x=x,y=y,fill = after_stat(level)),
                  geom = "polygon",adjust=2,alpha=.4) +
  geom_point(x=aim,y=0,alpha=.95,size=2) +
  scale_fill_gradient(guide = F,low="grey",high="black") +
  coord_cartesian(xlim=c(-2,2),ylim=c(-2,2)) +
  theme_tufte(base_family="sans",base_size=12) +
  theme(axis.line=element_line(size=.25),
        axis.ticks = element_line(size=.25),
        legend.position = "none")
ggsave("visuals/raw/fake_2d_aim.pdf",units="in",height=3,width=3,dpi=150)

# FAKE 1D REACH DISTRIBUTION
x_dens = with(density(data$x,adjust=2),tibble(x,y))
x_dens %>% 
  mutate(gain = ifelse(x<1&x>-1,x,0),
         loss = ifelse(x<0&x>-2,x,0)) %>%
  ggplot(aes(x,y)) +
  geom_line(aes(colour=y)) +
  geom_area(aes(x = gain),fill="lightskyblue",alpha=.95) +
  geom_area(aes(x = loss),fill="orange2",alpha=.75) +
  coord_cartesian(xlim=c(-2,2),ylim=c(0,.75)) +
  scale_colour_gradient(guide = F,low="grey",high="black") +
  theme_tufte(base_family="sans",base_size=12) +
  theme(axis.line.x=element_line(size=.25),
        axis.ticks.x = element_line(size=.25),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
ggsave("visuals/raw/fake_1d_aim.pdf",units="in",height=3,width=3,dpi=150)




# FAKE OUTCOME PROB BY AIM
get_p_gain = function(x,s=5) {
  aim = c(x,0)
  sig = diag(2)/s
  eig = diag(2)
  p = pmvnEll(r = 1, x0 = c(0,0), e = eig, mu = aim, sigma = sig)
  return(p)
}
get_p_loss= function(x,s=5) {
  aim = c(x,0)
  sig = diag(2)/s
  eig = diag(2)
  p = pmvnEll(r = 1, x0 = c(-1,0), e = eig, mu = aim, sigma = sig)
  return(p)
}
aims = seq(-1,1,0.01)
p_gain = sapply(X=aims, FUN=get_p_gain)
p_loss = sapply(X=aims, FUN=get_p_loss)
tibble(x=rep(aims,2),
  y=c(p_gain,p_loss),
  o=c(rep("gain",length(p_gain)),rep("loss",length(p_loss)))) %>%
  ggplot(aes(x=x,y=y,colour=o,group=o)) +
  geom_line(size=2) +
  scale_x_continuous("Aim") +
  scale_y_continuous("Outcome\nProbability") +
  scale_colour_manual("Outcome",values=c("lightskyblue","orange2")) +
  scale_fill_gradient_tableau(palette = "Gray") +
  coord_cartesian(xlim=c(-1,1),ylim=c(0,1)) +
  theme_tufte(base_family="sans",base_size=12) +
  theme(axis.line =element_line(size=.25),
        axis.ticks.x = element_line(size=.25),
        legend.position = "top")
ggsave("visuals/raw/fake_p_outcome.pdf",units="in",height=3,width=3,dpi=150)


# FAKE EXPECTED GAIN BY AIM
get_E_gain = function(x,s=5,l=1) {
  p_gain = get_p_gain(x,s)
  p_loss = get_p_loss(x,s)
  E_gain = 1*p_gain - l*p_loss
  return(E_gain)
}
aims = seq(-1,1,0.01)
E_gain = sapply(X=aims, FUN=get_E_gain)
optimal = aims[which.max(E_gain)]
tibble(x=aims,y=E_gain) %>%
  ggplot(aes(x=x,y=y,colour=y)) +
  geom_line(size=2) +
  geom_point(x=optimal,y=max(E_gain),size=2,inherit.aes = F) +
  scale_x_continuous("Aim") +
  scale_y_continuous("Expected\nGain ($)") +
  scale_colour_steps("$",low="orange2",high="lightskyblue",n.breaks=3) +
  coord_cartesian(xlim=c(-1,1),ylim=c(-1,1)) +
  theme_tufte(base_family="sans",base_size=12) +
  theme(axis.line =element_line(size=.25),
        axis.ticks.x = element_line(size=.25),
        legend.position = "top",
        legend.box.spacing = unit(0,"in"),
        legend.text = element_text(size=9))
ggsave("visuals/raw/fake_E_gain.pdf",units="in",height=3,width=3,dpi=150)


# FAKE EXPECTED GAIN BY AIM AND SIGMA
s_data = tibble()
aims = seq(-1,1,0.01)
for (s in 1:5) {
  E_gain_fun = Curry(FUN = get_E_gain,s=s) 
  s_data = s_data %>%
    bind_rows(
      tibble(
        x=aims,
        y=sapply(X=aims,FUN=E_gain_fun),
        s=rep(s,length(aims)),
        max_E_gain = rep(max(y),length(aims)),
        max_E_gain_x = rep(aims[which.max(y)],length(aims))
      )
    )
}
s_data_summary = s_data %>%
  group_by(s) %>%
  summarise(x=mean(max_E_gain_x),
            y=mean(max_E_gain))
s_data %>%
  ggplot(aes(x=x,y=y,colour=y,alpha=s,group=s)) +
  geom_line(size=2) +
  geom_point(data=s_data_summary,aes(x=x,y=y),size=2,inherit.aes = F) +
  scale_x_continuous("Aim") +
  scale_y_continuous("Expected\nGain ($)") +
  scale_colour_steps("$",low="orange2",high="lightskyblue",n.breaks=3) +
  scale_alpha_continuous(guide=F) +
  coord_cartesian(xlim=c(-1,1),ylim=c(-1,1)) +
  theme_tufte(base_family="sans",base_size=12) +
  theme(axis.line =element_line(size=.25),
        axis.ticks.x = element_line(size=.25),
        legend.position = "top",
        legend.box.spacing = unit(0,"in"),
        legend.text = element_text(size=9))
ggsave("visuals/raw/fake_sigma_E_gain.pdf",units="in",height=3,width=3,dpi=150)



# FAKE EXPECTED GAIN BY AIM AND P:R RATIO
l_data = tibble()
for (loss in seq(from=0,to=2,by=.5)) {
  E_gain_fun = Curry(FUN = get_E_gain,l=loss) 
  l_data = l_data %>%
    bind_rows(
      tibble(
        x=aims,
        y=sapply(X=aims,FUN=E_gain_fun),
        l=rep(loss,length(aims)),
        max_E_gain = rep(max(y),length(aims)),
        max_E_gain_x = rep(aims[which.max(y)],length(aims))
      )
    )
}

l_data_summary = l_data %>%
  group_by(l) %>%
  summarise(x=mean(max_E_gain_x),
            y=mean(max_E_gain))
l_data %>%
  ggplot(aes(x=x,y=y,colour=y,alpha=l,group=l)) +
  geom_line(size=2) +
  geom_point(data=l_data_summary,aes(x=x,y=y),size=2,inherit.aes = F) +
  scale_x_continuous("Aim") +
  scale_y_continuous("Expected\nGain ($)") +
  scale_colour_steps("$",low="orange2",high="lightskyblue",n.breaks=4) +
  scale_alpha_continuous(guide=F) +
  # coord_cartesian(xlim=c(-1,1),ylim=c(-1,1)) +
  theme_tufte(base_family="sans",base_size=12) +
  theme(axis.line =element_line(size=.25),
        axis.ticks.x = element_line(size=.25),
        legend.position = "top",
        legend.box.spacing = unit(0,"in"),
        legend.text = element_text(size=9))
ggsave("visuals/raw/fake_ratio_E_gain.pdf",units="in",height=3,width=3,dpi=150)
