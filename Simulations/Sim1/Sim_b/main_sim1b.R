

rm(list = ls(all=TRUE))
cat("\f")

set.seed(0)

library(mvtnorm) # dmvnorm
library(MASS) # mvrnorm
library(Matrix) # norm
library(BayesTree)
library(zoo) # coredata ->  to replicate
library(parallel)
library(np)
# detectCores()

args <- commandArgs(TRUE)
cat("\n Rscript <main.R> <n> <p> <d> <case> <weight> <augment> <iterations> \n \n")
# if (length(args) != 3) stop("usage: Rscript script.R arg1 arg2 arg3 ...")

n <- as.numeric(args[[1]])
p <- as.numeric(args[[2]])
d <- as.numeric(args[[3]])
case <- as.numeric(args[[4]])
weight <- args[[5]]
augment <- args[[6]]
iterations <- as.numeric(args[[7]])

# n = 100 # number of data points
# p = 6 # number of covariates
# d = 2 # number of basis for central space (beta: p*d)
# case = 2 # (1)distributional assumption, (2)no distributional assumption
# 
# # TRUE: E[{p(A)/p(A|C}*U(beta)] = 0
# # FALSE: E[U(beta)] = 0
# weight = TRUE
# 
# # TRUE: E[{p(A)/p(A|C}*U(beta) - sth + sth]
# # FALSE: E[{p(A)/p(A|C}*U(beta)] = 0
# augment = TRUE  # if weight = FALSE, augment won't matter anymore
# 
# # Number of for loop iterations
# iterations = 1

setwd("?")

source("data_sim1b.R")
source("analysis_sim1b.R")
source("derivatives_sim1b.R")
source("propensity_sim1b.R")
source("regression_sim1b.R")


# +++++++++++++++++++++++++++++++
# Initializations 
# +++++++++++++++++++++++++++++++

# Some Newton-Raphson args 
threshold = 0.5 # used for stopping criterion 
delta = 0.001 # used in computing the derivatives
mu = 0.001 # used in parameter update (H = H + mu*I_p)
lambda = 1 # step size in update rule: beta = beta - lambda*H^(-1)*J {or lambda/iter}

# bart vs dmvnorm
method_A = "dmvnorm"

# Kernels bart vs GLM
method_Y = "GLM"

# samp_size A from p(A)
samp_size = 100


# Store the Euclid differences 
if (weight == FALSE){
 file_name = paste0("diff_noIPW_n", n, "_p", p, "_case", case, 
                    "_", method_A, "_", method_Y, ".csv")
}else{
 if (augment == FALSE){
  file_name = paste0("diff_IPW_n", n, "_p", p, "_case", case, 
                     "_", method_A, "_", method_Y, ".csv")
 }else{
  file_name = paste0("diff_AIPW_n", n, "_p", p, "_case", case, 
                     "_", method_A, "_", method_Y, ".csv")
 }
}

f_pca = paste0("PCA_", p, "_", case, ".csv") 


# +++++++++++++++++++++++++++++++
# Run the analysis 
# +++++++++++++++++++++++++++++++
collection_beta = matrix(0, p*d, iterations) # collect all the beta

# f <- function(i){
# for (i in 1:iterations){
  # cat("( ****** iteration #", i, "****** ) \n")

  # ++++++++++++++
  # Generate data
  # ++++++++++++++
  output = generate_data(n, p, case)
  dat = output$dat
  A = output$A
  C = output$C
  Y = output$Y
  true_beta = output$true_beta  # true A coeff in E[Y | A, C]

  # +++++++++++++++
  # Pick functions
  # +++++++++++++++
  # alpha(A)
  alpha = list()
  for (i in 1:nrow(A)){
   alpha[[i]] = 0.01*A[i, ]%*%t(A[i, ])
   # alpha[[i]] = cbind(A[i, ], A[i, ])
  }

  # p(a)/p(a|c) -> propensity_sim2.R
  Sigma = cov(dat)
  obj_A = get_fA_fAC(n, p, case, C, A, Y, dat, method_A, samp_size) # bart vs dmvnorm
  fA_fAC = obj_A$fA_fAC
  A_sample = obj_A$A_sample

  # E[Y | A, C]
  obj_Y = get_regression(n, p, C, A, Y, dat, method_Y, samp_size, A_sample, true_beta) # bart vs Kernels vs GLM
  y_ac = obj_Y$y_ac
  Yrep = obj_Y$Yrep

  # +++++++++++++++++
  # Run the analysis
  # +++++++++++++++++
  beta = matrix(round(rnorm(p*d, 0, 1), 3), nrow=p)
  # beta = matrix(cbind(rep(1, p), rep(-1, p)), nrow=p)/sqrt(p)

  bandw = F
  # h = n^(-1/(d+4))
  h = sd(A)*((0.75*n)^(-1/(d+4)))
  h = 5
  
  results = run_analysis(beta, A_sample, threshold, delta, mu, lambda, A, C, Y, 
                         y_ac, Yrep, Sigma, alpha, fA_fAC, weight, augment, h, bandw)

  (est_beta = results$est_beta)
  # collection_beta[, i] = as.vector(est_beta)
# }
  
# write.csv(collection_beta, file=file_name, row.names=FALSE)
  
# ++++++++++++++++++++++++++++++++++
# Compare ...
# ++++++++++++++++++++++++++++++++++
 
# Yq = A%*%true_beta[, 1] + (A%*%true_beta[, 2])^2
Yq = A%*%true_beta[, 1] + A^2%*%true_beta[, 2]
diff = c()
for (r in seq(0.1, 50, 0.1)){
 Yq.hat_truth = estimate_kernel(Y, A, A, true_beta, r, bandw=F)
 head(data.frame(Yq, Yq.hat_truth), 20)
 diff = c(diff, mean((Yq - Yq.hat_truth)^2))
}
min(diff)



Yq.hat = estimate_kernel(Y, A, A, est_beta, h , bandw)
datq = data.frame(Yq, Yq.hat)
head(datq, 20)
sqrt(sum((Yq - Yq.hat)^2))
  
  
  
  
  
