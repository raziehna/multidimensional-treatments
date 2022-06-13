
rm(list = ls(all=TRUE))
cat("\f")

library(mvtnorm) # dmvnorm
library(MASS) # mvrnorm
library(reshape2) # melt

setwd("?")
source("data_sim1.R")

set.seed(0) 

# n = 200 # number of data points
# p = 6 # number of covariates
# case = 2 # (1)distributional assumption, (2)no distributional assumption

file_noIPW_6_1 = paste0("diff_noIPW_n", 200, "_p", 6, "_case", 1, "_dmvnorm_GLM.csv")
file_IPW_6_1 = paste0("diff_IPW_n", 200, "_p", 6, "_case", 1, "_dmvnorm_GLM.csv")
file_AIPW_6_1 = paste0("diff_AIPW_n", 200, "_p", 6, "_case", 1, "_dmvnorm_GLM.csv")
file_PCA_6_1 = "PCA_6_1.csv"

file_noIPW_6_2 = paste0("diff_noIPW_n", 200, "_p", 6, "_case", 2, "_dmvnorm_GLM.csv")
file_IPW_6_2 = paste0("diff_IPW_n", 200, "_p", 6, "_case", 2, "_dmvnorm_GLM.csv")
file_AIPW_6_2 = paste0("diff_AIPW_n", 200, "_p", 6, "_case", 2, "_dmvnorm_GLM.csv")
file_PCA_6_2 = "PCA_6_2.csv"

file_noIPW_12_1 = paste0("diff_noIPW_n", 200, "_p", 12, "_case", 1, "_dmvnorm_GLM.csv")
file_IPW_12_1 = paste0("diff_IPW_n", 200, "_p", 12, "_case", 1, "_dmvnorm_GLM.csv")
file_AIPW_12_1 = paste0("diff_AIPW_n", 200, "_p", 12, "_case", 1, "_dmvnorm_GLM.csv")
file_PCA_12_1 = "PCA_12_1.csv"

file_noIPW_12_2 = paste0("diff_noIPW_n", 200, "_p", 12, "_case", 2, "_dmvnorm_GLM.csv")
file_IPW_12_2 = paste0("diff_IPW_n", 200, "_p", 12, "_case", 2, "_dmvnorm_GLM.csv")
file_AIPW_12_2 = paste0("diff_AIPW_n", 200, "_p", 12, "_case", 2, "_dmvnorm_GLM.csv")
file_PCA_12_2 = "PCA_12_2.csv"

noIPW_6_1 = read.csv(file_noIPW_6_1, header = TRUE)[, 3]
IPW_6_1 = read.csv(file_IPW_6_1, header = TRUE)[, 3]
AIPW_6_1 = read.csv(file_AIPW_6_1, header = TRUE)[, 3]
PCA_6_1 = read.csv(file_PCA_6_1, header = TRUE)[, 1]

noIPW_6_2 = read.csv(file_noIPW_6_2, header = TRUE)[, 3]
IPW_6_2 = read.csv(file_IPW_6_2, header = TRUE)[, 3]
AIPW_6_2 = read.csv(file_AIPW_6_2, header = TRUE)[, 3]
PCA_6_2 = read.csv(file_PCA_6_2, header = TRUE)[, 1]

noIPW_12_1 = read.csv(file_noIPW_12_1, header = TRUE)[, 3]
IPW_12_1 = read.csv(file_IPW_12_1, header = TRUE)[, 3]
AIPW_12_1 = read.csv(file_AIPW_12_1, header = TRUE)[, 3]
PCA_12_1 = read.csv(file_PCA_12_1, header = TRUE)[, 1]

noIPW_12_2 = read.csv(file_noIPW_12_2, header = TRUE)[, 3]
IPW_12_2 = read.csv(file_IPW_12_2, header = TRUE)[, 3]
AIPW_12_2 = read.csv(file_AIPW_12_2, header = TRUE)[, 3]
PCA_12_2 = read.csv(file_PCA_12_2, header = TRUE)[, 1]


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# P = 6 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

dat = data.frame(rbind(noIPW_6_1, IPW_6_1, AIPW_6_1, PCA_6_1, 
                       noIPW_6_2, IPW_6_2, AIPW_6_2, PCA_6_2))
dat$type = c("0Regres", "1IPW", "2AIPW", "3PCA", 
             "0Regres", "1IPW", "2AIPW", "3PCA")
dat$case = c("Case 01_p06", "Case 01_p06", "Case 01_p06", "Case 01_p06",
             "Case 02_p06", "Case 02_p06", "Case 02_p06", "Case 02_p06")

stacked.dat = melt(dat, id = c('type', 'case'))
stacked.dat = stacked.dat[, -3]

# png('boxplots_p6.png')
png(filename="boxplots_p6.png", 
    units="in", 
    width=6,
    height=8,
    pointsize=12,
    res=200)
boxplots.triple = boxplot(value~type + case, 
                          data = stacked.dat, 
                          las = 2,
                          names = c('Reg','IPW', 'AIPW', 'PCA',
                                    'Reg','IPW', 'AIPW', 'PCA'),
                          at = c(1:4, 8:11), 
                          outline=FALSE, #xaxt='n',
                          yaxt='n', width=rep(1, 8),
                          ylim = c(0, 2.05), col = c('gray', 'gray', 'gray', 'gray'))

# axis(side=1, at=c(2.5, 9.5), line=2, lwd=0,
#      labels=c('\n Case 1 \n ', '\n Case 2 \n'))
# axis(side=1, at=c(2.5, 5, 8.5, 16.5, 20, 22.5), line=0.5, lwd=0, 
#      labels=c('\n Case 1', 'p = 6', '\n Case 2', 
#               '\n Case 1', 'p = 12', '\n Case 2')
# )

axis(side=2, at=seq(0, 2, by = 0.5)) 
title('Boxplots of Euclidean Distances \n (n = 200, p = 6)')


# draw the rectangles for Los Angeles, specifying the left, bottom, right, and top positions
# of the box plots, the density of the lines, and the angles of the lines
# rect(c(1.4, 6.4), boxplots.triple$stats[2, c(2, 5)],
#      c(2.2, 7.2), boxplots.triple$stats[4, c(2, 5)],
#      density=12, angle=45)
# add city labels near the box plots

text(c(2.5, 9.5),
     c(2.05, 2.05),
     c('Case (1)', 'Case (2)'),
     font = 1)


y = seq(0, 2.1, by = 0.1)
x = rep(6, length(y))
lines(x, y, type="l", lty=2)

dev.off()


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# P = 12
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


dat = data.frame(rbind(noIPW_12_1, IPW_12_1, AIPW_12_1, PCA_12_1, 
                       noIPW_12_2, IPW_12_2, AIPW_12_2, PCA_12_2))
dat$type = c("0Regres", "1IPW", "2AIPW", "3PCA", 
             "0Regres", "1IPW", "2AIPW", "3PCA")
dat$case = c("Case 1_p12", "Case 1_p12", "Case 1_p12", "Case 1_p12",
             "Case 2_p12", "Case 2_p12", "Case 2_p12", "Case 2_p12")

stacked.dat = melt(dat, id = c('type', 'case'))
stacked.dat = stacked.dat[, -3]

# png('boxplots_p12.jpeg')
png(filename="boxplots_p12.png", 
    units="in", 
    width=6,
    height=8,
    pointsize=12,
    res=200)
boxplots.triple = boxplot(value~type + case, 
                          data = stacked.dat, 
                          las = 2,
                          names = c('Reg','IPW', 'AIPW', 'PCA',
                                    'Reg','IPW', 'AIPW', 'PCA'),
                          at = c(1:4, 8:11), 
                          outline=FALSE, # xaxt='n', 
                          yaxt='n', width=rep(1, 8),
                          ylim = c(0, 2.05), col = c('gray', 'gray', 'gray', 'gray'))

# axis(side=1, at=c(2.5, 9.5), line=0.5, lwd=0,
#      labels=c('Case 1 \n ', 'Case 2 \n '))
# axis(side=1, at=c(2.5, 5, 8.5, 16.5, 20, 22.5), line=0.5, lwd=0, 
#      labels=c('\n Case 1', 'p = 6', '\n Case 2', 
#               '\n Case 1', 'p = 12', '\n Case 2')
# )

axis(side=2, at=seq(0, 2, by = 0.5)) 
title('Boxplots of Euclidean Distances \n (n = 200, p = 12)') 

# draw the rectangles for Los Angeles, specifying the left, bottom, right, and top positions
# of the box plots, the density of the lines, and the angles of the lines
# rect(c(1.4, 6.4), boxplots.triple$stats[2, c(2, 5)],
#      c(2.2, 7.2), boxplots.triple$stats[4, c(2, 5)],
#      density=12, angle=45)
# add city labels near the box plots

text(c(2.5, 9.5),
     c(2.05, 2.05),
     c('Case (1)', 'Case (2)'),
     font = 1)

# text(c(1:4, 8:11), 
#      c(9.7, 4.8, 0.7, 2.7, 
#        4.6, 4.3, 0.5, 2.4),
#      c('Reg','IPW', 'AIPW', 'PCA',
#        'Reg','IPW', 'AIPW', 'PCA'), 
#      font = 1)


y = seq(0, 2.1, by = 0.1)
x = rep(6, length(y))
lines(x, y, type="l", lty=2)

dev.off()


