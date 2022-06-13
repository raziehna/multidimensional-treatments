

pseudoinverse <- function(beta, eps){
 beta%*%solve(t(beta)%*%beta + eps*diag(ncol(beta)))%*%t(beta)
}

get_r_beta <- function(beta, A_sample, p, d, A, C, Y, y_ac, Yrep,
                       Sigma, alpha, fA_fAC, weight, augment){
  n = length(Y) 
  beta = matrix(beta, nrow=p)
  d = ncol(beta)
 
  # Get U and E[U | A, C]
  y_aTb = rep(0, n)
  for (i in 1:d){
    y_aTb = y_aTb + A^i%*%beta[, i]
  }
  U = list() # list of p*d matrices
  U_AC = list() 
  for(i in 1:n){
    U[[i]] = (Y[i] - y_aTb[i])*(alpha[[i]])
    U_AC[[i]] = (y_ac[i] - y_aTb[i])*(alpha[[i]])
  }
  
  # Get Eq[ E[U|A,C] | C ] = \int_A E[U | A, C] p(A) = (1/samp_size)*
  if (augment == TRUE){
   samp_size = nrow(A_sample)
   
   # alpha_samp = matrix( rep( A_sample , d ) , nrow =  nrow(A_sample) , byrow = FALSE )
   alp = rep(1, samp_size) 
   for(j in 1:d) alp = cbind(alp, (A_sample)^j) 
   alp = alp[, -1]
   alpha_samp = matrix( alp , nrow =  nrow(A_sample) , byrow = FALSE )
   y_aTb_samp = rep(0, samp_size)
   for (i in 1:d) y_aTb_samp = y_aTb_samp + A_sample^i%*%beta[, i]
   q_U_AC = list()
   for (i in 1:n){
    y_ac_sample = Yrep[((i-1)*samp_size+1):(i*samp_size)]
    q_U_AC[[i]] = t(t(y_ac_sample - y_aTb_samp)%*%alpha_samp)/samp_size 
    q_U_AC[[i]] = matrix( q_U_AC[[i]], nrow = p)
   }
  }
  
  # Get r(beta)
  temp = matrix(rep(0, p*d), nrow = p)
  for (i in 1:n){
    if (weight == FALSE){
      s = U[[i]]
    }else{
       if (augment == FALSE){
         s = fA_fAC[i]*U[[i]] 
       }else{
         # s = fA_fAC[i]*(2*U[[i]] - U_AC[[i]])
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
                                  Sigma, alpha, fA_fAC, delta, weight, augment){

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
                           Sigma, alpha, fA_fAC, weight, augment)
    r_bk_minus = get_r_beta(bk_minus, A_sample, p, d, A, C, Y, y_ac, Yrep,
                            Sigma, alpha, fA_fAC, weight, augment)
    
    # r_bk_plus = matrix(r_bk_plus, nrow = p)
    # r_bk_minus = matrix(r_bk_minus, nrow = p)
    
    # first_d = (norm(r_bk_plus, type="f")^2 - norm(r_bk_minus, type="f")^2)/(2*delta)
    first_d = (sum((r_bk_plus)^2) - sum((r_bk_minus)^2))/(2*delta)
    # first_d = sum((r_bk_plus - r_bk_minus)^2)/(2*delta)
    deriv = c(deriv, first_d)
  }
  return(deriv)
}


get_second_derivatives <- function(beta, A_sample, p, d, A, C, Y, y_ac, Yrep,
                                   Sigma, alpha, fA_fAC, delta, weight, augment){

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
                                      Sigma, alpha, fA_fAC, delta, weight, augment)
    r_bk_minus = get_first_derivatives(bk_minus, A_sample, p, d, A, C, Y, y_ac, Yrep,
                                       Sigma, alpha, fA_fAC, delta, weight, augment)
    
    second_deriv[k,] = {r_bk_plus - r_bk_minus}/(2*delta)
  }
  return(second_deriv)
}


update_beta <- function(beta, first_deriv, second_deriv, mu, lambda){
  p = length(first_deriv)
  second_deriv = second_deriv + mu*diag(p)
  vec_beta_new= as.vector(beta) - lambda*solve(second_deriv, tol = 1e-24)%*%first_deriv
  beta_new = matrix(vec_beta_new, nrow = nrow(beta))
  # beta_new = as.matrix(vec_beta_new, nrow = nrow(beta))
  return(beta_new)
}


