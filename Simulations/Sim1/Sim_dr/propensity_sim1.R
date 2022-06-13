


get_fA_fAC <- function(n, p, case, C, A, Y, dat, method,  samp_size, miss_A){
 
 # Estimate the covariance matrix 
 Sigma = cov(dat)
 
 # p(A)
 # Razie: keep in mind that p(A) can be anything but maybe it would be more efficient if mean A is said to colMeans(A)
 sigA = Sigma[5:(4+p), 5:(4+p)]
 pA = dmvnorm(A, mean=rep(0, p), sigma=sigA)
 # pA = rep(1/n, n)
 
 # Pick samples from p(A)
 A_sample = mvrnorm(samp_size, rep(0, p), sigA)
 # idx = sample(1:n, samp_size, replace = TRUE)
 # A_sample = A[idx, ]
 
 # p(A | C)
 sig11 = Sigma[1:4, 1:4]
 sig12 = Sigma[5:(4+p), 1:4]
 sig21 = Sigma[1:4, 5:(4+p)]
 sig22 = Sigma[5:(4+p), 5:(4+p)]
 sigAC = sig22 - sig12%*%solve(sig11)%*%sig21 
 if (method == "bart"){
  mean_A = matrix(0, nrow = n, ncol = p)
  for (i in 1:p){
   mod_A = bart(C, A[, i], C, ndpost=300, nskip=100, verbose=FALSE)
   if (i == 5 | i == 6){
    mean_A[, i] = pnorm(colMeans(mod_A$yhat.test))
   }else{
    mean_A[, i] = mod_A$yhat.test.mean
   }
  }
  pA_C = rep(0, n)
  for (i in 1:n) pA_C[i] = dmvnorm(A[i, ], mean = mean_A[i, ], sigma=sigAC)
  
 }else if (method == "dmvnorm" & case == 2){
  muA = C%*%t(sig12%*%solve(sig11))
  pA_C = rep(0, n)
  
  if(miss_A){
    muA[, 1] = 0.5*muA[, 1]
    muA[, 2] = 0.5*muA[, 2]
    muA[, 3] = 0.5*muA[, 3]
  }
  for (i in 1:n) pA_C[i] = dmvnorm(A[i, ], mean = muA[i, ], sigma=sigAC)
  
 }else if ( method == "dmvnorm" & p == 6){ # case = 1
  C1 = C[, 1]
  C2 = C[, 2]
  C3 = C[, 3]
  C4 = C[, 4]
  mu_A1 = C1 + C2 + C3 + C4
  mu_A2 = - C1 + C2 - C3 + C4 
  mu_A = data.frame(A1=mu_A1, A2=mu_A2)
  mu_A = apply(mu_A, 2, as.vector)
  fA1A2 = rep(0, n)
  fA3 = rep(0, n)
  fA4 = rep(0, n)
  fA5 = rep(0, n)
  fA6 = rep(0, n)
  # for (i in 1:n){
  #  fA1A2[i] = dmvnorm(A[i, 1:2], mean = mu_A[i, ], sigma=sigAC[1:2, 1:2])
  #  fA3[i] = dnorm(A[i,3], mean=abs(A[i,1]+A[i,2]), sd=sqrt(abs(A[i,1]))) + dnorm(A[i,3], mean=-abs(A[i,1]+A[i,2]), sd=sqrt(abs(A[i,1])))
  #  fA4[i] = dnorm(A[i,4], mean=sqrt(abs(A[i,1]+A[i,2])), sd=sqrt(abs(A[i,2]))) + dnorm(A[i,4], mean=-sqrt(abs(A[i,1]+A[i,2])), sd=sqrt(abs(A[i,2])))
  # }
  for (i in 1:n){
   fA1A2[i] = dmvnorm(A[i, 1:2], mean = mu_A[i, ], sigma=sigAC[1:2, 1:2])
   fA3[i] = dnorm(A[i,3], mean=abs(A[i,1]+A[i,2]), sd=sqrt(abs(A[i, 1]))) 
   fA4[i] = dnorm(A[i,4], mean=sqrt(abs(A[i,1]+A[i,2])), sd=sqrt(abs(A[i, 2])))
  }
  idx_A5 = which(A[, 5] == 1)
  fA5[idx_A5] = exp(A[idx_A5,2])/(1 + exp(A[idx_A5,2]))
  fA5[-idx_A5] = 1 - exp(A[-idx_A5,2])/(1 + exp(A[-idx_A5,2]))
  idx_A6 = which(A[, 6] == 1)
  fA6[idx_A6] = pnorm(A[idx_A6,2])
  fA6[-idx_A6] = 1 - pnorm(A[-idx_A6,2])
  pA_C = fA1A2*fA3*fA4*fA5*fA5*fA6
  
 }else { # method == "dmvnorm" & p == 12 & case = 1
  C1 = C[, 1]
  C2 = C[, 2]
  C3 = C[, 3]
  C4 = C[, 4]
  mu_A1 = C1 + C2 + C3 + C4
  mu_A2 = - C1 + C2 - C3 + C4 
  mu_A7 = C1
  mu_A8 = C2
  mu_A9 = C3
  mu_A10 = - C1 + C2
  mu_A11 = - C2 + C3 
  mu_A12 = - C3 + C4
  mu_A = data.frame(A1=mu_A1, A2=mu_A2, A7=mu_A7, A8=mu_A8, A9=mu_A9, A10=mu_A10, A11=mu_A11, A12=mu_A12)
  mu_A = apply(mu_A, 2, as.vector)
  fA1A2A7_A12 = rep(0, n)
  fA3 = rep(0, n)
  fA4 = rep(0, n)
  fA5 = rep(0, n)
  fA6 = rep(0, n)
  for (i in 1:n){
   fA1A2A7_A12[i] = dmvnorm(A[i, c(1:2, 7:12)], mean = mu_A[i, ], sigma=sigAC[c(1:2, 7:12), c(1:2, 7:12)])
   fA3[i] = dnorm(A[i,3], mean=abs(A[i,1]+A[i,2]), sd=sqrt(abs(A[i, 1])))
   fA4[i] = dnorm(A[i,4], mean=sqrt(abs(A[i,1]+A[i,2])), sd=sqrt(abs(A[i, 2])))
  }
  idx_A5 = which(A[, 5] == 1)
  fA5[idx_A5] = exp(A[idx_A5,2])/(1 + exp(A[idx_A5,2]))
  fA5[-idx_A5] = 1 - exp(A[-idx_A5,2])/(1 + exp(A[-idx_A5,2]))
  idx_A6 = which(A[, 6] == 1)
  fA6[idx_A6] = pnorm(A[idx_A6,2])
  fA6[-idx_A6] = 1 - pnorm(A[-idx_A6,2])
  pA_C = fA1A2A7_A12*fA3*fA4*fA5*fA5*fA6
 }
 
 
 # p(a)/p(a|c)
 fA_fAC = pA/pA_C
 cat('\n pA/pA_C: ', summary(round(fA_fAC, 3)), '\n')
 object = boxplot(fA_fAC, plot=FALSE)
 up_bound = min(object$out)
 # up_bound = 10*quantile(fA_fAC, 0.75)
 if (up_bound > max(fA_fAC)) up_bound = max(fA_fAC)
 cat("percentage of rows being greater than up_bound: ", mean(fA_fAC > up_bound), "\n")
 fA_fAC[fA_fAC <= 0.0001] = 0.0001
 fA_fAC[fA_fAC >= up_bound] = up_bound
 cat(' pA/pA_C: ', summary(round(fA_fAC, 3)), '\n')
 
 
 return(list(fA_fAC=fA_fAC, A_sample=A_sample)) 
 
}


