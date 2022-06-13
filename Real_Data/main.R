
rm(list = ls(all=TRUE))
cat("\f")

set.seed(0) 

# args <- commandArgs(TRUE)
# cat("\n Rscript <main.R> <p> <B> \n \n")
# 
# p <- as.numeric(args[[1]])
# B <- as.numeric(args[[2]])

p = 5
B = 50 # number of bootstraps
max_dim = 1

library(MASS)
library(BayesTree)
library(mvtnorm)
library(zoo)
library(np) 
library(mice)

# setwd("?")

source("analysis.R")
source("derivatives.R")
source("propensity.R")
source("regression.R") 

# +++++++++++++++++++++++++++++++
# Data Preparation  
# +++++++++++++++++++++++++++++++

# prepare data
dat = read.csv("radiation.csv", stringsAsFactors = TRUE)
C = dat[, 1:ncol(C)]
A = dat[, (ncol(C)+1):(ncol(C)+p)]
Y = dat[, ncol(dat)]

# +++++++++++++++++++++++++++++++
# Initialization  
# +++++++++++++++++++++++++++++++

n = nrow(dat)
p = ncol(A)

weight = TRUE
augment = FALSE 
bandw = FALSE 

method_A = "bart"
method_Y = "bart"

samp_size = 20
threshold = 0.1 
delta = 0.01
mu = 0.01 
lambda = 1 

if (weight == FALSE){
 file_name = paste0("diff_noIPW_n", n, ".csv")
}else{
 if (augment == FALSE){
  file_name = paste0("diff_IPW_n", n, ".csv")
 }else{
  file_name = paste0("diff_AIPW_n", n, ".csv")
 }
}

# +++++++++++++++++++++++++++++++++++
# Run the analysis - original sample
# +++++++++++++++++++++++++++++++++++

# p(a)/p(a|c) -> propensity.R
Sigma = cov(dat)
obj_A = get_fA_fAC(C, A, Y, dat, method_A, samp_size) # bart vs dmvnorm 
fA_fAC = obj_A$fA_fAC
A_sample = obj_A$A_sample

# E[Y | A, C] -> regression.R
obj_Y = get_regression(C, A, Y, dat, method_Y, samp_size, A_sample) # bart vs GLM vs Kernels
Y_ac = obj_Y$Y_ac
Yrep = obj_Y$Yrep 

est_beta_0 = list()
for (d in 1:max_dim){
 cat("\n ( ****** original sample with d #", d, "****** ) \n")

 # alpha(A)
 alpha = list()
 for (i in 1:nrow(A)){
  alp = rep(1, p)
  for(j in 1:d) alp = cbind(alp, (A[i, ])^j)
  alpha[[i]] = alp[, -1]
 }

 # bandwidth
 h = (n^(-1/(d+4)))

 # initial beta
 beta = matrix(round(rnorm(p*d, 0, 1), 3), nrow=p)

 # run the analysis
 results_0 = run_analysis(beta, A_sample, threshold, delta, mu, lambda, A, C, Y, Y_ac, Yrep,
                          Sigma, alpha, fA_fAC, weight, augment, h, bandw)

 est_beta_0[[d]] = matrix(results_0$est_beta, nrow = p)
}

for (d in 1:max_dim){
 write.csv(est_beta_0[[d]],
           file=paste0("beta_IPW_d", d, "_p", p, ".csv"),
           row.names=FALSE)
}

# est_beta_0 = list()
# for (d in 1:max_dim){
#  est_beta_0[[d]] = as.matrix(read.csv(paste0("beta_IPW_d", d, "_p", p, ".csv")), nrow = p)
# }


# +++++++++++++++++++++++++++++++++++
# Run the analysis - bootstrap sample
# +++++++++++++++++++++++++++++++++++
 
collection_dim = matrix(0, B, max_dim)

for (b in 1:B){
 start_time <- Sys.time()
 cat("\n ( ****** Bootstrap #", b, "****** ) \n")
 
 idx = sample(n, n, replace = TRUE)
 dat_boots = dat[idx, ]
 row.names(dat_boots) = 1:n
 C = as.matrix(dat_boots[, 1:ncol(C)])
 A = as.matrix(dat_boots[, (ncol(C)+1):(ncol(C) + p)])
 Y = dat_boots$Y

 Sigma = cov(dat_boots)
 obj_A = get_fA_fAC(C, A, Y, dat_boots, method_A, samp_size) # bart vs dmvnorm
 fA_fAC = obj_A$fA_fAC
 A_sample = obj_A$A_sample

 obj_Y = get_regression(C, A, Y, dat_boots, method_Y, samp_size, A_sample) # bart vs GLM
 y_ac = obj_Y$Y_ac
 Yrep = obj_Y$Yrep

 for (d in 1:max_dim){
  cat("\n ( *** b#", b, ", d#", d, " *** ) \n")

  # alpha(A)
  alpha = list()
  for (i in 1:nrow(A)){
   alp = rep(1, p)
   for(j in 1:d) alp = cbind(alp, (A[i, ])^j)
   alpha[[i]] = alp[, -1]
  }

  # bandwidth
  h = (n^(-1/(d+4)))

  # initial beta
  beta = matrix(round(rnorm(p*d, 0, 1), 3), nrow=p)

  results = run_analysis(beta, A_sample, threshold, delta, mu, lambda, A, C, Y, Y_ac, Yrep,
                         Sigma, alpha, fA_fAC, weight, augment, h, bandw)
  est_beta = results$est_beta
  est_beta = matrix(est_beta, nrow = p)

  # Compute u, v -> u, v = (A*b1, A*b2) or A*b1, A^2*b2)
  u = A%*%est_beta_0[[d]]
  v = A%*%est_beta

  # Compute r2(u, v)
  u.var = var(u)
  v.var = var(v)
  uv.cov = cov(u, v)
  vu.cov = cov(v, u)
  u.eig = eigen(u.var)
  if (d == 1 ) u.diag = sqrt(u.eig$values)*diag(1) else u.diag = diag(sqrt(u.eig$values))
  u.sqrt = u.eig$vectors %*% u.diag %*% solve(u.eig$vectors)
  u.sqrt.inv = solve(u.sqrt)
  v.inv = solve(v.var)
  L = u.sqrt.inv%*%uv.cov%*%v.inv%*%vu.cov%*%u.sqrt.inv

  eigV = eigen(L)$values
  cat("eigV ", eigV, "\n")
  id_eig = which(eigV == 0)
  if (sum(id_eig) != 0 ){
   if ( length(id_eig) == d ) eigV = 0 else eigV = eigV[-id_eig]
  }
  collection_dim[b, d] = mean(eigV)
 }
 
 cat("\n")
 print(collection_dim[1:b, ])
 cat("\n")
 
 end_time <- Sys.time()
 cat(end_time - start_time, "\n")
}

write.csv(collection_dim, file=file_name, row.names=FALSE)

