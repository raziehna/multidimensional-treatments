
# Causal Sufficient Dimension Reduction - pick the dimension

rm(list = ls(all=TRUE))
cat("\f")

args <- commandArgs(TRUE)
cat("\n Rscript <main.R> <n> <p> <case> <B> \n \n")

n <- as.numeric(args[[1]])
p <- as.numeric(args[[2]])
case <- as.numeric(args[[3]])
B <- as.numeric(args[[4]])

library(mvtnorm) # dmvnorm
library(MASS) # mvrnorm, ginv
library(Matrix) # norm
library(BayesTree)
library(zoo) # coredata ->  to replicate
library(parallel)

# setwd("?")

set.seed(0)

source("data_sim2.R")  
source("analysis_sim2.R")
source("derivatives_sim2.R")
source("propensity_sim2.R")
source("regression_sim2.R")

# +++++++++++++++++++++++++++++++
# Initializations 
# +++++++++++++++++++++++++++++++

# TRUE: E[{p(A)/p(A|C}*U(beta)] = 0
# FALSE: E[U(beta)] = 0
weight = TRUE

# TRUE: E[{p(A)/p(A|C}*U(beta) - sth + sth]
# FALSE: E[{p(A)/p(A|C}*U(beta)] = 0
augment = TRUE  # if weight = FALSE, augment won't matter anymore

# bart vs dmvnorm
method_A = "dmvnorm"

# bart vs GLM
method_Y = "GLM"

# sample A from p(A)
samp_size = 100

# Some Newton-Raphson args 
threshold = 1 # used for stopping criterion 
delta = 0.001 # used in computing the derivatives
mu = 0.01 # used in parameter update (H = H + mu*I_p)
lambda = 1 # step size in update rule: beta = beta - lambda*H^(-1)*J {or lambda/iter}

# Store the Euclid differences 
if (weight == FALSE){
 file_name = paste0("dim_noIPW_n", n, "_p", p, "_case", case, 
                    "_", method_A, "_", method_Y, ".csv")
}else{
 if (augment == FALSE){
  file_name = paste0("dim_IPW_n", n, "_p", p, "_case", case, 
                     "_", method_A, "_", method_Y, ".csv")
 }else{
  file_name = paste0("dim_AIPW_n", n, "_p", p, "_case", case, 
                     "_", method_A, "_", method_Y, ".csv")
 }
}

# +++++++++++++++++++++++++++++++++++
# Run the analysis - original sample
# +++++++++++++++++++++++++++++++++++

# Generate data
output = generate_data(n, p, case)
dat_orig = output$dat 
C = output$C
A = output$A
Y = output$Y
true_beta = output$true_beta  # true A coeff in E[Y | A, C]

# p(a)/p(a|c) -> propensity_sim3.R
Sigma = cov(dat_orig)
obj_A = get_fA_fAC(n, p, case, C, A, Y, dat_orig, method_A, samp_size) # bart vs dmvnorm 
fA_fAC = obj_A$fA_fAC
A_sample = obj_A$A_sample

# E[Y | A, C] -> regression_sim3.R 
obj_Y = get_regression(n, p, C, A, Y, dat, method_Y, samp_size, A_sample) # bart vs GLM
y_ac = obj_Y$y_ac
Yrep = obj_Y$Yrep 

est_beta_0 = list()
for (d in 1:(p-1)){
 # d = 5
 cat("\n ( ****** original sample with d #", d, "****** ) \n")
 
 # alpha(A)
 alpha = list() 
 # for (i in 1:nrow(A)) alpha[[i]] = replicate(d, A[i, ])
 for (i in 1:nrow(A)){
  alp = rep(1, p)
  for(j in 1:d) alp = cbind(alp, (A[i, ])^j) 
  alpha[[i]] = alp[, -1]
 } 
 
 # initial beta
 beta = matrix(round(rnorm(p*d, 0, 1), 3), nrow=p)
 # beta = matrix(rep(1, p*d), nrow=p)
 
 # run the analysis
 results_0 = run_analysis(beta, A_sample, threshold, delta, mu, lambda, A, C, Y, 
                          y_ac, Yrep, Sigma, alpha, fA_fAC, weight, augment)
 # cat("\n", results_0$est_beta, "\n")
 
 est_beta_0[[d]] = matrix(results_0$est_beta, nrow = p)  
}


# +++++++++++++++++++++++++++++++++++
# Run the analysis - bootstrap sample
# +++++++++++++++++++++++++++++++++++
collection_dim = matrix(0, B, p-1)
for (b in 1:B){
 cat("\n ( ****** Bootstrap #", b, "****** ) \n")
 # sample rows to construct the bootstrap sample
 idx = sample(n, n, replace = TRUE)
 dat = dat_orig[idx, ]
 row.names(dat) = 1:n
 C = as.matrix(dat[, 1:4])
 A = as.matrix(dat[, 5:10])
 Y = dat$Y

 # p(a)/p(a|c) -> propensity_sim3.R
 Sigma = cov(dat)
 obj_A = get_fA_fAC(n, p, case, C, A, Y, dat, method_A, samp_size) # bart vs dmvnorm
 fA_fAC = obj_A$fA_fAC
 A_sample = obj_A$A_sample

 # E[Y | A, C] -> regression_sim3.R
 obj_Y = get_regression(n, p, C, A, Y, dat, method_Y, samp_size, A_sample) # bart vs GLM
 y_ac = obj_Y$y_ac
 Yrep = obj_Y$Yrep

 for (d in 1:(p-1)){
   # d = 5

   cat("\n ( *** b#", b, ", d#", d, " *** ) \n")

   # alpha(A)
   alpha = list() 
   # for (i in 1:nrow(A)) alpha[[i]] = replicate(d, A[i, ])
   for (i in 1:nrow(A)){
     alp = rep(1, p)
     for(j in 1:d) alp = cbind(alp, (A[i, ])^j) 
     alpha[[i]] = alp[, -1]
   } 
   
   # initial beta
   beta = matrix(round(rnorm(p*d, 0, 1), 3), nrow=p)
   # beta = matrix(rep(1, p*d), nrow=p)

   results = run_analysis(beta, A_sample, threshold, delta, mu, lambda, A, C, Y,
                          y_ac, Yrep, Sigma, alpha, fA_fAC, weight, augment)
   est_beta = results$est_beta
   est_beta = matrix(est_beta, nrow = p)

   # Compute u, v
   u = matrix(0, n, d)
   v = matrix(0, n, d)
   for (i in 1:d){
     u[, i] = A^i%*%est_beta_0[[d]][, i]
     v[, i] = A^i%*%est_beta[, i]
   }
   
   # Compute r2(u, v)
   u.var = var(u)
   v.var = var(v)
   uv.cov = cov(u, v)
   vu.cov = cov(v, u)
   u.eig = eigen(u.var)
   if (d ==1 ) u.diag = sqrt(u.eig$values)*diag(1) else u.diag = diag(sqrt(u.eig$values))
   u.sqrt = u.eig$vectors %*% u.diag %*% solve(u.eig$vectors)
   u.sqrt.inv = solve(u.sqrt)
   v.inv = solve(v.var)
   L = u.sqrt.inv%*%uv.cov%*%v.inv%*%vu.cov%*%u.sqrt.inv
   
   eigV = eigen(L)$values
   cat("eigV ", eigV, "\n")
   id_eig = which(round(eigV, 2) == 0)
   if (length(id_eig) !=0 ) eigV = eigV[-id_eig]
   collection_dim[b, d] = mean(eigV)
 }
 cat("\n")
 print(collection_dim[1:b, ])
 cat("\n")
}

write.csv(collection_dim, file=file_name, row.names=FALSE)


