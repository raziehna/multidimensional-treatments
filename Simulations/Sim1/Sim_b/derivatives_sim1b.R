

pseudoinverse <- function(beta, eps){
 beta%*%solve(t(beta)%*%beta + eps*diag(ncol(beta)))%*%t(beta)
}


estimate_kernel <- function(Y, A_sample, A, beta, h, bandw){
 N = nrow(A)
 S = nrow(A_sample)
 D = ncol(beta)
 
 count_denom = 0
 # if (bandw == TRUE) h = sd(A%*%beta)*((0.75*N)^(-1/(D+4)))
 
 y_est = rep(0, S)
 for (s in 1:S){

  K_u = 1
  for (d in 1:D){
   # diff = as.numeric(A_sample[s, ]^d%*%beta[, d]) - A^d%*%beta[, d]  # diff: d times N
   diff = as.numeric(A_sample[s, ]%*%beta[, d]) - A%*%beta[, d]  # diff: d times N
   
   if (bandw == TRUE) h = sd(A%*%beta[, d])*((0.75*N)^(-1/(d+4)))
   
   # u = diff[, d]/h
   u = diff/h
   k1_u = (1/h)*(3/4)*(1 - u^2)*(abs(u) < 1)
   # k4_u = (15/18)*(1 - (7/3)*u^2)*k1_u
   # K_u = K_u*k4_u
   K_u = K_u*k1_u
  }
  num = sum(Y*K_u)
  denom = sum(K_u)
  
  if (denom == 0){
   y_est[s] = Y[s]
   cat("\n K_u is all zero, s is", s, " \n" )
   count_denom = count_denom + 1
  }else y_est[s] = num/denom
 }
 # cat("count_denom (K_u being all zero) = ", count_denom, "\n")
 return(y_est)
}


get_r_beta <- function(beta, A_sample, p, d, A, C, Y, y_ac, Yrep, 
                       Sigma, alpha, fA_fAC, weight, augment, h, bandw){
  
  n = length(Y) 
  beta = matrix(beta, nrow=p)
  d = ncol(beta)

  # Get U and E[U | A, C]
  y_aTb = estimate_kernel(Y, A, A, beta, h, bandw)
  U = list() # list of p*d matrices
  U_AC = list() 
  for(i in 1:n){
    U[[i]] = (Y[i] - y_aTb[i])*(alpha[[i]])
    U_AC[[i]] = (y_ac[i] - y_aTb[i])*(alpha[[i]])
  }
  
  # Get Eq[ E[U|A,C] | C ] = \int_A E[U | A, C] p(A) = (1/samp_size)*
  if (augment == TRUE){
   samp_size = nrow(A_sample)
   
   
   alpha_samp = list()
   for (i in 1:samp_size){
    alpha_samp[[i]] = 0.01*A_sample[i, ]%*%t(A_sample[i, ])
    # alpha[[i]] = cbind(A[i, ], A[i, ])
   }
   
   
   # # alpha_samp = matrix( rep( A_sample , d ) , nrow =  nrow(A_sample) , byrow = FALSE )
   # alp = rep(1, samp_size)
   # for(j in 1:d){
   #  # alp = cbind(alp, sqrt(abs(A_sample)))
   #  alp = cbind(alp, (0.1)^d*A_sample)
   # }
   # alp = alp[, -1]
   # alpha_samp = matrix( alp , nrow =  nrow(A_sample) , byrow = FALSE )
   
   
   y_aTb_samp = estimate_kernel(Y, A_sample, A, beta, h, bandw)
   # y_aTb_samp = rep(0, samp_size)
   
   q_U_AC = list()
   for (i in 1:n){
    y_ac_sample = Yrep[((i-1)*samp_size+1):(i*samp_size)]
    q = 0
    for (j in 1:samp_size){
     q = q + (y_ac_sample[j] - y_aTb_samp[j])*alpha_samp[[j]]
    }
    q_U_AC[[i]] = q/samp_size
    # q_U_AC[[i]] = t(t(y_ac_sample - y_aTb_samp)%*%alpha_samp)/samp_size 
    # q_U_AC[[i]] = matrix( q_U_AC[[i]], nrow = p)
   }
  }
  
  # Get r(beta)
  # temp = matrix(rep(0, p*d), nrow = p)
  temp = matrix(rep(0, p*p), nrow = p)
  for (i in 1:n){
    if (weight == FALSE){
      s = U[[i]]
    }else{
       if (augment == FALSE){
         s = fA_fAC[i]*U[[i]] 
       }else{
          s = fA_fAC[i]*(U[[i]] - U_AC[[i]]) + q_U_AC[[i]]
       }
    }
    temp = s + temp
  }
  r_beta = temp/n
  r_beta = as.vector(r_beta)
  return(r_beta)
}


get_first_derivatives <- function(beta, A_sample, p, d, A, C, Y, y_ac, Yrep, 
                                  Sigma, alpha, fA_fAC, delta, weight, augment, h, bandw){
  beta = as.vector(beta)
  K = length(beta)
  
  deriv = c()
  for (k in 1:K){
    e = matrix(rep(0, K), ncol = d)
    # e = as.matrix(rep(0, K), ncol = d)
    e[k] = 1
    
    bk_plus = beta + delta*e 
    bk_minus = beta - delta*e 
    
    r_bk_plus = get_r_beta(bk_plus, A_sample, p, d, A, C, Y, y_ac, Yrep, 
                           Sigma, alpha, fA_fAC, weight, augment, h, bandw)
    r_bk_minus = get_r_beta(bk_minus, A_sample, p, d, A, C, Y, y_ac, Yrep, 
                            Sigma, alpha, fA_fAC, weight, augment, h, bandw)

    # first_d = (norm(r_bk_plus, type="f")^2 - norm(r_bk_minus, type="f")^2)/(2*delta)
    first_d = (sum((r_bk_plus)^2) - sum((r_bk_minus)^2))/(2*delta)
    deriv = c(deriv, first_d)
  }
  return(deriv)
}


get_second_derivatives <- function(beta, A_sample, p, d, A, C, Y, y_ac, Yrep, 
                                   Sigma, alpha, fA_fAC, delta, weight, augment, h, bandw){

  beta = as.vector(beta) 
  K = length(beta)
  
  second_deriv = matrix(rep(0, K^2), nrow = K)
  for (k in 1:K){
    e = matrix(rep(0, K), ncol = d)
    # e = as.matrix(rep(0, K), ncol = d)
    e[k] = 1
    
    bk_plus = beta + delta*e
    bk_minus = beta - delta*e 
    
    r_bk_plus = get_first_derivatives(bk_plus, A_sample, p, d, A, C, Y, y_ac, Yrep, 
                                      Sigma, alpha, fA_fAC, delta, weight, augment, h, bandw)
    r_bk_minus = get_first_derivatives(bk_minus, A_sample, p, d, A, C, Y, y_ac, Yrep, 
                                       Sigma, alpha, fA_fAC, delta, weight, augment, h, bandw)
    
    second_deriv[k,] = {r_bk_plus - r_bk_minus}/(2*delta)
  }
  return(second_deriv)
}

update_beta <- function(beta, first_deriv, second_deriv, mu, step_size){
  pd = length(first_deriv)
  second_deriv = second_deriv + mu*diag(pd)
  vec_beta_new= as.vector(beta) - step_size*solve(second_deriv, tol = 1e-20)%*%first_deriv
  beta_new = matrix(vec_beta_new, nrow = nrow(beta))
  return(beta_new)
}


