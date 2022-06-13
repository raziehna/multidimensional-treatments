
rm(list = ls(all=TRUE))
cat("\f")

library(mvtnorm) # dmvnorm
library(MASS) # mvrnorm
library(reshape2) # melt

setwd("?")

set.seed(0) 

file_missAY = "missAY_diff_AIPW_n200_p6_case2_dmvnorm_GLM.csv"
file_missA = "missA_mu12_diff_AIPW_n200_p6_case2_dmvnorm_GLM.csv"
file_missY = "missY_diff_AIPW_n200_p6_case2_dmvnorm_GLM.csv"
file_nomiss = "no_diff_AIPW_n200_p6_case2_dmvnorm_GLM.csv"

missAY = read.csv(file_missAY, header = TRUE)[, 3]
missA = read.csv(file_missA, header = TRUE)[, 3]
missY = read.csv(file_missY, header = TRUE)[, 3]
nomiss = read.csv(file_nomiss, header = TRUE)[, 3]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# P = 6 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

dat = data.frame(rbind(missAY, missA, missY, nomiss))
dat$type = c("0missAY", "1missA", "2missY", "3nomiss")
dat$case = c("Case 02_p06", "Case 02_p06", "Case 02_p06", "Case 02_p06")

stacked.dat = melt(dat, id = c('type', 'case'))
stacked.dat = stacked.dat[, -3]

png(filename="boxplots_p6.png", 
    units="in", 
    width=6,
    height=8,
    pointsize=12,
    res=200)
boxplots.triple = boxplot(value~type + case, 
                          data = stacked.dat, 
                          las = 2,
                          names = c('Miss_AY','Miss_A', 'Miss_Y', 'Miss_no'),
                          at = c(1:4), 
                          outline=FALSE, #xaxt='n',
                          yaxt='n', width=rep(1, 4),
                          ylim = c(0, 1.05), 
                          xlab="", ylab="",
                          col = c('gray', 'gray', 'gray', 'gray'))

axis(side=2, at=seq(0, 2, by = 0.25)) 
title('Boxplots of Euclidean Distances \n (n = 200, p = 6, Case 2)')

dev.off()



