
# 
# 
# 
# 
# # TEST
# 
# 
# tmp_beta <- 0.8315077
# tmp_gamma <- 1.9993164
# j = 1
# c = 10
# cond <- 1:10
# mc <- c()
# tmp <- data %>% filter(sub_id == j)
# obs <- tmp$x
# cond_id <- tmp$cond_id
# sig <- sd(obs)
# all_p <- fm[,,j,]
# P <- all_p
# p <- P[,,2]
# hyp <- seq(0,1,.01) # aims
# # get outcome probs based on separation distance
# pg <- p[1,]; pgl <- p[2,]; pl <- p[3,]; pn <- p[4,]
# l <- 0
# g <- 1
# uloss <- get_u(l,tmp_beta)
# uboth <- g + uloss
# 
# wp <- get_w(p,tmp_gamma)
# swu <- get_swu(wp,g,uboth,uloss)
# print(swu)
# mu <- hyp[which.max(swu)]
# 
# 
# 
# hyp <- seq(0,1,.01)
# N <- length(hyp)
# tmp <- tibble(x = rep(hyp,4),
#        y = c(fm[1,,1,2],fm[2,,1,2],fm[3,,1,2],fm[4,,1,2]),
#        o = c(rep("gain",N),rep("both",N),rep("loss",N),rep("miss",N)))
# tmp %>% ggplot(aes(x=x,y=y,colour=o,group=o)) +
#   geom_line(size=1.5) +
#   scale_y_continuous("Objective Probability") +
#   scale_colour_brewer("Outcome",palette = "Set1") +
#   coord_cartesian(xlim=c(0,1),ylim=c(0,1))
# 
# 
# hyp <- seq(0,1,.01)
# N <- length(hyp)
# tmp <- tibble(x = rep(hyp,4),
#               y = c(get_w(fm[1,,1,2],tmp_gamma),
#                     get_w(fm[2,,1,2],tmp_gamma),
#                     get_w(fm[3,,1,2],tmp_gamma),
#                     get_w(fm[4,,1,2],tmp_gamma)),
#               o = c(rep("gain",N),rep("both",N),rep("loss",N),rep("miss",N)))
# 
# tmp2 <- tibble(x = hyp,
#                y = get_w(fm[1,,1,2],tmp_gamma) + get_w(fm[2,,1,2],tmp_gamma))
# 
# tmp %>% ggplot(aes(x=x,y=y,colour=o,group=o)) +
#   geom_line(size=1.5) +
#   geom_line(data=tmp2,aes(x=x,y=y),colour="black",linetype="dashed",inherit.aes = F) +
#   annotate(geom="point",
#            x=hyp[which.max(get_w(fm[1,,1,2],tmp_gamma) + get_w(fm[2,,1,2],tmp_gamma))],
#            y=max(get_w(fm[1,,1,2],tmp_gamma) + get_w(fm[2,,1,2],tmp_gamma)),
#            colour="red",size=3) +
#   scale_y_continuous("Subjective Probability") +
#   scale_colour_brewer("Outcome",palette = "Set1") +
#   coord_cartesian(xlim=c(0,1),ylim=c(0,1))  
# 
# hyp <- seq(0,1,.01)
# N <- length(hyp)
# 
# tmp %>% ggplot(aes(x=x,y=y)) +
# 
#   scale_y_continuous("p(gain) + p(both)") +
#   scale_colour_brewer("Outcome",palette = "Set1") +
#   coord_cartesian(xlim=c(0,1),ylim=c(0,1))  
# 
# 
# 
# 
# 
# 
# 
# 
# 
