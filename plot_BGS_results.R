rm(list = ls())

library(ggplot2)
library(reshape2)

bgs <- read.csv("~/work/Mike_BGS/run1_temp/FirstTry.csv")

pi_0 <- mean(bgs[bgs$Nes==0,]$window_1)

bgs_m <- melt(bgs, id = c("Nes", "rep", "N", "r"))
bgs_m$B <- bgs_m$value/pi_0

bgs_m <- bgs_m[bgs_m$Nes!=0,]

bgs_m$labs <- factor(bgs_m$variable,
                     levels = c("window_1", "window_2"),
                     labels = c(expression(italic(r) > 0),
                                expression(italic(r)*" = 0")))

ggplot(data = bgs_m)+
  stat_summary(aes( x = Nes*2, y = B))+
  scale_x_log10(expression(italic( "2N"[e]*"s") ) )+
  scale_y_continuous(expression(italic(B)))+
  coord_cartesian(ylim = c(0.8, 1.1))+
  facet_grid(~labs,
             labeller = label_parsed)+
  theme_bw()
  

  