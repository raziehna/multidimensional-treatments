
get_fA_fAC <- function(C, A, Y, dat, method, samp_size){
 
 n = dim(C)[1]
 p = dim(A)[2]
 k = dim(C)[2]
 
 # Estimate the covariance matrix 
 Sigma = cov(dat)
 
 # p(A) *************
 sigA = Sigma[(k+1):(k+p), (k+1):(k+p)]
 pA = rep(1/n, n)
 
 # Pick samples from p(A)
 idx = sample(1:n, samp_size, replace = TRUE)
 A_sample = A[idx, ]
 
 
 # p(A | C) *************
 sig11 = Sigma[1:k, 1:k]
 sig12 = Sigma[(k+1):(k+p), 1:k]
 sig21 = Sigma[1:k, (k+1):(k+p)]
 sig22 = Sigma[(k+1):(k+p), (k+1):(k+p)]
 sigAC = sig22 - sig12%*%solve(sig11)%*%sig21 
 if (method == "bart"){
  mean_AC = matrix(0, nrow = n, ncol = p)
  for (i in 1:p){
   xtrain = apply(C, 2, as.double)
   mod = bart(x.train=xtrain, y.train=A[, i], x.test=xtrain, ndpost=100, nskip=100, verbose=FALSE)
   mean_AC[, i] = mod$yhat.test.mean
  }
  pA_C = rep(0, n)
  for (i in 1:n) pA_C[i] = dmvnorm(A[i, ], mean = mean_AC[i, ], sigma=sigAC)
  # truncate 
  pA_C[pA_C < 0.0001] = 0.0001
  
  }else{ # "dmvnorm"
  muA = as.matrix(C)%*%t(sig12%*%solve(sig11))
  pA_C = rep(0, n)
  for (i in 1:n) pA_C[i] = dmvnorm(A[i, ], mean = muA[i, ], sigma=sigAC)
  }
 
 # p(a)/p(a|c) *************
 fA_fAC = pA/pA_C
 cat('\n pA/pA_C: ', summary(round(fA_fAC, 3)), '\n')
 object = boxplot(fA_fAC, plot=FALSE)
 sort_out = sort(object$out)
 max = sort_out[97]
 min = sort_out[58]
 fA_fAC[fA_fAC < min] = min
 fA_fAC[fA_fAC > max] = max

 return(list(fA_fAC=fA_fAC, 
             A_sample=A_sample)) 
 
}


