library(mvtnorm)
library(raster)
library(tidyverse)

par1 <- read_rds("data/clean/exp1/MSWU_param.rds")
par2 <- read_rds("data/clean/exp2/MSWU_param.rds")

tibble(a1 = par1) %>%
  summarise_all(.funs=c("mu"=mean,"med"=median,"min"=min,"max"=max,"sd"=sd))
tibble(a2 = par2) %>%
  summarise_all(.funs=c("mu"=mean,"med"=median,"min"=min,"max"=max,"sd"=sd))


# PLOTS

for (plt in c("exp1","exp2")) {
  
  if (plt =="demo")  {
    v <- seq(-5,1,.1)
    v_breaks <- c(-5,-1,0,0,1)
    v_lim <- c(-5,1)
    alphas <- seq(0.5,2,.1)
    betas <- seq(0.5,2,.1)
    gammas <- seq(0.5,2,.1)
  } else if (plt == "exp1") {
    v <- seq(-15,3,.1)
    #v_breaks <- c(-15,-5,-3,-1,1,3,5)
    v_lim <- c(-15,5)
    alphas <- par1[1,]
    betas <- par1[2,]
    gammas <- par1[3,]
  } else {
    v <- seq(-5,1,.1)
    #v_breaks <- c(-5,-1,0,1,3)
    v_lim <- c(-5,3)
    alphas <- par2[1,]
    betas <- par2[2,]
    gammas <- par2[3,]
  }
  
  adata <- tibble()
  n <- length(v)
  for (i in 1:length(alphas)) {
    alpha <- alphas[i]; beta <- betas[i]
    get_u <- function(v,a=alpha,b=beta) {
      # non-linear value transform
      if (v < 0){
          u <- -((-v)^a)
      } else {
        if (v > 0) {
          u <- v^b
        } else {
          u <- 0
        }
      }
      return(u)
    }
    u <- sapply(v, get_u)
    adata <- tibble(a = rep(alpha,n),b = rep(beta,n),v = v,u = u) %>% rbind(adata)
  }
  adata %>%
    ggplot(aes(x=v,y=u,group=a,colour=a)) +
    geom_abline(size=.25,linetype="dashed") +
    geom_vline(size=.25,xintercept=0,linetype="dashed") +
    geom_hline(size=.25,yintercept=0,linetype="dashed") +
    geom_line(size=.25) +
    scale_colour_viridis_c(expression(alpha),option="C") +
    scale_x_continuous("Value") +
    scale_y_continuous("Utility") +
    coord_fixed(xlim=v_lim,ylim=v_lim) +
    theme_tufte(base_size=10,base_family="sans") +
    theme(legend.position = "top",
          axis.line = element_line(size=.25),
          legend.title = element_text(vjust = 1),
          legend.key.height = unit(x = .1,units = "in"),
          legend.margin = margin(r=15))
  ggsave(sprintf("visuals/raw/%s/revision/u_function.pdf",plt),units="in",height=2,width=2)
  
  
  gdata <- tibble()
  p <- seq(0,1,.01)
  n <- length(p)
  for (g in gammas) {
    get_w <- function(x,gamma=g) {
      return(exp(-(-log(x))^gamma))
    }
    w <- sapply(p, get_w)
    gdata <- tibble(g = rep(g,n),p = p,w = w) %>% rbind(gdata)
  }
  gdata %>%
    ggplot(aes(x=p,y=w,group=g,colour=g)) +
    geom_abline(size=.25,linetype="dashed") +
    geom_line(size=.25) +
    scale_colour_viridis_c(expression(gamma),option = "C") +
    scale_x_continuous("Probability",breaks=c(0,.25,.5,.75,1), labels=c("0",".25",".5",".75","1")) +
    scale_y_continuous("Weight",breaks=c(0,.25,.5,.75,1), labels=c("0",".25",".5",".75","1")) +
    coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
    theme_tufte(base_size=10,base_family="sans") +
    theme(legend.position = "top",
          axis.line = element_line(size=.25),
          legend.title = element_text(vjust = 1),
          legend.key.height = unit(x = .1,units = "in"),
          legend.margin = margin(r=15))
  ggsave(sprintf("visuals/raw/%s/revision/w_functions.pdf",plt),units="in",height=2,width=2)
}




