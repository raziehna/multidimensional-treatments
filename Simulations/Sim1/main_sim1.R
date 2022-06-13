
rm(list = ls(all=TRUE))
cat("\f")

library(mvtnorm) # dmvnorm
library(MASS) # mvrnorm
library(Matrix) # norm
library(BayesTree)
library(zoo) # coredata ->  to replicate
library(parallel)

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

threshold = 0.1 # used for stopping criterion 

set.seed(0)

source("data_sim1.R")
source("analysis_sim1.R")
source("derivatives_sim1.R")
source("propensity_sim1.R")
source("regression_sim1.R")


# +++++++++++++++++++++++++++++++
# Initializations 
# +++++++++++++++++++++++++++++++

# Some Newton-Raphson args 
delta = 0.001 # used in computing the derivatives
mu = 0.01 # used in parameter update (H = H + mu*I_p)
lambda = 1 # step size in update rule: beta = beta - lambda*H^(-1)*J {or lambda/iter}

# bart vs dmvnorm
method_A = "dmvnorm"

# bart vs GLM
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

f_pca = paste0("PCA_n", n, "_p", p, "_case", case, "_", ".csv") 


# +++++++++++++++++++++++++++++++
# Run the analysis 
# +++++++++++++++++++++++++++++++
collection_diff_0 = c() # collect all the frob_diff between the true_param and the estimated_param
collection_diff_1 = c()
collection_diff_2 = c()
collection_diff_pca = c()

# f <- function(i){
for (i in 1:iterations){
  cat("( ****** iteration #", i, "****** ) \n")

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
   alpha[[i]] = cbind(A[i, ], A[i, ]^2)
  }

  # p(a)/p(a|c) -> propensity_sim2.R
  Sigma = cov(dat)
  obj_A = get_fA_fAC(n, p, case, C, A, Y, dat, method_A, samp_size) # bart vs dmvnorm
  fA_fAC = obj_A$fA_fAC
  A_sample = obj_A$A_sample

  # E[Y | A, C]
  obj_Y = get_regression(n, p, C, A, Y, dat, method_Y, samp_size, A_sample) # bart vs GLM
  y_ac = obj_Y$y_ac
  Yrep = obj_Y$Yrep

  # +++++++++++++++++
  # Run the analysis
  # +++++++++++++++++
  # beta = matrix(round(rnorm(p*d, 0, 1), 3), nrow=p)
  beta = matrix(cbind(rep(1, p), rep(-1, p)), nrow=p)

  results = run_analysis(beta, A_sample, threshold, delta, mu, lambda, A, C, Y,
                         y_ac, Yrep, Sigma, alpha, fA_fAC, weight, augment)

  est_beta = results$est_beta
  est_beta = matrix(est_beta, nrow = p)

  diff_0 = norm(true_beta - est_beta, type="f")

  diff_1 =  abs(norm(est_beta, type="f") - norm(true_beta, type="f"))

  est_beta_frob = est_beta%*%solve( t(est_beta)%*%est_beta )%*%t(est_beta)
  true_beta_frob = true_beta%*%solve( t(true_beta)%*%true_beta )%*%t(true_beta)
  diff_2 = norm(est_beta_frob - true_beta_frob, type="f")


  collection_diff_0 = c(collection_diff_0, diff_0)
  collection_diff_1 = c(collection_diff_1, diff_1)
  collection_diff_2 = c(collection_diff_2, diff_2)
  cat("\n Euclid_diff_0: ", diff_0)
  cat("\n Euclid_diff_1: ", diff_1)
  cat("\n Euclid_diff_2: ", diff_2, " \n \n ")
  cat("+++++++++++++++++++++++++++++++++++++++++ \n \n ")

  # # perform pca
  # A_sdv = svd(A)
  # beta_pca = A_sdv$v[, 1:d]
  # beta_pca_frob = beta_pca%*%solve( t(beta_pca)%*%beta_pca )%*%t(beta_pca)
  # true_beta_frob = true_beta%*%solve( t(true_beta)%*%true_beta )%*%t(true_beta)
  # diff_pca = norm(beta_pca_frob - true_beta_frob, type="f")
  # collection_diff_pca = c(collection_diff_pca, diff_pca)

}

collection_diff = data.frame(normDiff = collection_diff_0,
                             diffNorm = collection_diff_1,
                             frobNorm = collection_diff_2)
write.csv(collection_diff, file=file_name, row.names=FALSE)
# write.csv(collection_diff_pca, file=f_pca, row.names=FALSE)



