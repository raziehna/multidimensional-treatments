
# ++++++++++++++++++++++++++++++++++
# Data generation function 
# ++++++++++++++++++++++++++++++++++
generate_data <- function(n, p, case){
  
 # Generate Covariate
  C = mvrnorm(n, rep(0, 4), diag(4))
  C1 = C[, 1]
  C2 = C[, 2]
  C3 = C[, 3]
  C4 = C[, 4]
  colnames(C) = paste0("C", 1:4)
  
  # Generate Treatment
  if (case == 1 & p == 6){
   mu = rep(0, 2)
   sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2)
   A = mvrnorm(n, mu, sigma)
   eps_1 = rnorm(n, 0, 1)
   eps_2 = rnorm(n, 0, 1)
   A1 = A[, 1] + C1 + C2 + C3 + C4 
   A2 = A[, 2] - C1 + C2 - C3 + C4
   A3 = abs(rnorm(n, abs(A1 + A2), sd = sqrt(abs(A1))))
   A4 = abs(rnorm(n, sqrt(abs(A1 + A2)), sd = sqrt(abs(A2))))
   p_A5 = exp(A2)/(1 + exp(A2))
   A5 = rbinom(n, 1, p_A5)
   p_A6 = pnorm(A2)
   A6 = rbinom(n, 1, p_A6)
   A = cbind(A1, A2, A3, A4, A5, A6) 
  } else if (case == 1 & p == 12){
   mu = rep(0, 8)
   sigma = diag(8)
   for (i in 1:(8-1)){
    for (j in (i+1):8){
     s = 0.5^{abs(i - j)}
     sigma[i, j] = s
     sigma[j, i] = s
    }
   }
   A = mvrnorm(n, mu, sigma)
   A1 = A[, 1] + C1 + C2 + C3 + C4 
   A2 = A[, 2] - C1 + C2 - C3 + C4
   A7 = A[, 3] + C1
   A8 = A[, 4] + C2
   A9 = A[, 5] + C3
   A10 = A[, 6] -C1 + C2
   A11 = A[, 7] -C2 + C3 
   A12 = A[, 8] - C3 + C4
   A3 = abs(rnorm(n, abs(A1 + A2), sd = sqrt(abs(A1))))
   A4 = abs(rnorm(n, sqrt(abs(A1 + A2)), sd = sqrt(abs(A2))))
   p_A5 = exp(A2)/(1 + exp(A2))
   A5 = rbinom(n, 1, p_A5)
   p_A6 = pnorm(A2)
   A6 = rbinom(n, 1, p_A6)
   A = cbind(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12)
   
  }else{ # case == 2
   mu = rep(0, p)
   sigma = diag(p)
   for (i in 1:(p-1)){
    for (j in (i+1):p){
     s = 0.5^{abs(i - j)}
     sigma[i, j] = s
     sigma[j, i] = s
    }
   }
   A = mvrnorm(n, mu, sigma)
   A[,1] = A[,1] + C1 + C2 + C3 + C4 
   A[,2] = A[,2] - C1 + C2 - C3 + C4 
   A[,3] = A[,3] + C1 - C2 - C3 + C4 
   A[,4] = A[,4] - C1 + C2 + C3 - C4 
   A[,5] = A[,5] + C1 + C2 - C3 + C4 
   A[,6] = A[,6] - C1 + C2 + C3 + C4 
   if (p == 12){
    A[,7] = A[,7] + C1 
    A[,8] = A[,8] + C2
    A[,9] = A[,9] + C3
    A[,10] = A[,10] - C1
    A[,11] = A[,11] - C2
    A[,12] = A[,12] - C3
   }
  }
  colnames(A) = paste0("A", 1:p)
  
  # Pick True Parameters
  if (p == 6){
   b1 = c(1, 1, 1, 1, 1, 1)/sqrt(6)
   b2 = c(1, -1, 1, -1, 1, -1)/sqrt(6)
   b3 = c(-1, 1, -1, 1, -1, 1)
  }else{
   b1 = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)/sqrt(6)
   b2 = c(1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 0, 0)/sqrt(6)
   b3 = c(-1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1)
  }
  true_beta = matrix(c(b1, b2), ncol = 2)
  
  # Generate Response
  eps = rnorm(n, 0, 1)
  Y = A%*%b1 + A^2%*%b2 + rowSums(C) + rowSums(A)*rowSums(C) + eps 
  # Y = A%*%b1 + A^2%*%b2 - 4*rowSums(C) - 4*rowSums(A)*rowSums(C) + eps
  # Y = A%*%b1 + A^2%*%b2 + 4*rowSums(C) + 2*(A%*%b3)*rowSums(C) + eps
  
  hist(Y)
  
  dat = data.frame(C, A, Y)
  return(list(dat = dat, 
              A = A, 
              Y = Y, 
              C = C, 
              true_beta = true_beta))
}

# Used for standardization
fnMatSqrtInverse <- function(mA){
  ei = eigen(mA)
  d = ei$values
  d = (d+abs(d))/2
  d2 = 1/sqrt(d)
  d2[d == 0] = 0
  return(ei$vectors %*% diag(d2) %*% t(ei$vectors))
}