library(tidyverse)

u <- function(v,a=1,b=.7,l=.75) {
  if (v >= 0) {
    return(v^a)
  } else {
    return(-l*(-v)^b)
  }
}

v = runif(1000,-5,1) 
u = sapply(v,u)

df = tibble(u,v)

df %>% ggplot(aes(v,u)) +
  geom_line() +
  geom_hline(yintercept = 0,
             linetype = "dotted") +
  scale_x_continuous(breaks=-5:1)
