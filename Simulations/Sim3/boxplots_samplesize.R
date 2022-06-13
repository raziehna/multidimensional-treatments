
rm(list = ls(all=TRUE))
cat("\f")

library(mvtnorm) # dmvnorm
library(MASS) # mvrnorm
library(reshape2)

setwd("?")

# set.seed(0) 

n = c(100, 200, 500, 1000, 5000) # number of data points
p = 6 # number of covariates
case = 2 # (1)distributional assumption, (2)no distributional assumption

file_IPW_100 = paste0("diff_IPW_n", 100, "_p", p, "_case", case, "_dmvnorm_GLM.csv")
file_AIPW_100 = paste0("diff_AIPW_n", 100, "_p", p, "_case", case, "_dmvnorm_GLM.csv")

file_IPW_200 = paste0("diff_IPW_n", 200, "_p", p, "_case", case, "_dmvnorm_GLM.csv")
file_AIPW_200 = paste0("diff_AIPW_n", 200, "_p", p, "_case", case, "_dmvnorm_GLM.csv")

file_IPW_500 = paste0("diff_IPW_n", 500, "_p", p, "_case", case, "_dmvnorm_GLM.csv")
file_AIPW_500 = paste0("diff_AIPW_n", 500, "_p", p, "_case", case, "_dmvnorm_GLM.csv")

file_IPW_1000 = paste0("diff_IPW_n", 1000, "_p", p, "_case", case, "_dmvnorm_GLM.csv")
file_AIPW_1000 = paste0("diff_AIPW_n", 1000, "_p", p, "_case", case, "_dmvnorm_GLM.csv")

file_IPW_5000 = paste0("diff_IPW_n", 5000, "_p", p, "_case", case, "_dmvnorm_GLM.csv")
file_AIPW_5000 = paste0("diff_AIPW_n", 5000, "_p", p, "_case", case, "_dmvnorm_GLM.csv")

IPW_100 = read.csv(file_IPW_100, header = TRUE)[, 3]
AIPW_100 = read.csv(file_AIPW_100, header = TRUE)[, 3]

IPW_200 = read.csv(file_IPW_200, header = TRUE)[, 3]
AIPW_200 = read.csv(file_AIPW_200, header = TRUE)[, 3]

IPW_500 = read.csv(file_IPW_500, header = TRUE)[, 3]
AIPW_500 = read.csv(file_AIPW_500, header = TRUE)[, 3]

IPW_1000 = read.csv(file_IPW_1000, header = TRUE)[, 3]
AIPW_1000 = read.csv(file_AIPW_1000, header = TRUE)[, 3]

IPW_5000 = read.csv(file_IPW_5000, header = TRUE)[, 3]
AIPW_5000 = read.csv(file_AIPW_5000, header = TRUE)[, 3]


dat = data.frame(rbind(IPW_100, AIPW_100, 
                       IPW_200, AIPW_200,
                       IPW_500, AIPW_500, 
                       IPW_1000, AIPW_1000,
                       IPW_5000, AIPW_5000
                       # c(IPW_5000, rep(NA, 50)), 
                       # c(AIPW_5000, rep(NA, 50))
                       ))
dat$type = rep(c("1IPW", "2AIPW"), 5)
dat$case = c("n = 0100", "n = 0100", "n = 0200",  "n = 0200",  "n = 0500", 
             "n = 0500", "n = 1000", "n = 1000", "n = 5000", "n = 5000")

stacked.dat = melt(dat, id = c('type', 'case'))
stacked.dat = stacked.dat[, -3] 


png(filename="boxplots_n.png",
    units="in",
    width=6,
    height=8,
    pointsize=12,
    res=200)
boxplots.triple = boxplot(value~type + case, 
                          data = stacked.dat, at = c(4, 6, 
                                                     10, 12, 
                                                     16, 18, 
                                                     22, 24, 
                                                     28, 30), 
                          outline=FALSE, xaxt='n', yaxt='n', width=rep(1, 10),
                          ylim = c(-0.05, 2), col = c('gray', 'gray'), 
                          xlab = "sample size (n)\n")
axis(side=1, at=c(5, 11, 17, 23, 29), line=0.5, lwd=0,
     labels=c('0.1k\n', '0.2k\n', '0.5k\n', '1.0k\n', '5.0k\n'))
axis(side=2, at=seq(0, 2, by = 0.2)) 
title('Boxplots of Euclidean Distances \n (case 2 with p = 6)')


# draw the rectangles for Los Angeles, specifying the left, bottom, right, and top positions
# of the box plots, the density of the lines, and the angles of the lines
# rect(c(1.4, 6.4), boxplots.triple$stats[2, c(2, 5)],
#      c(2.2, 7.2), boxplots.triple$stats[4, c(2, 5)],
#      density=12, angle=45)
# add city labels near the box plots
# text(c(4, 6, 
#        10, 12, 
#        16, 18, 
#        22, 24, 
#        28, 30),
#      c(3, 0.75, 
#        2.45, 0.5, 
#        1.9, 0.4, 
#        1.7, 0.25, 
#        1.35, 0.15),
#      rep(c('IPW','AIPW'), 5), 
#      font = 1)


x = rep(8, 21)
y = seq(0, 2, by = 0.1)
lines(x, y, type="l", lty=3)

x = rep(14, 21)
y = seq(0, 2, by = 0.1)
lines(x, y, type="l", lty=3)

x = rep(20, 21)
y = seq(0, 2, by = 0.1)
lines(x, y, type="l", lty=3)

x = rep(26, 21)
y = seq(0, 2, by = 0.1)
lines(x, y, type="l", lty=3)

dev.off()

