
pseudoinverse <- function(beta, eps){
 beta%*%solve(t(beta)%*%beta + eps*diag(ncol(beta)))%*%t(beta)
}

estimate_bart <- function(Y, Asamp, A, beta, h, bandw){
 mod_Y = bart(A, as.double(Y), Asamp, ndpost=100, nskip=100, verbose=FALSE)
 y_est = mod_Y$yhat.test.mean
 return(y_est)
}

estimate_kernel <- function(Y, Asamp, A, beta, h, bandw){
 N = nrow(A)
 S = nrow(Asamp)
 D = ncol(beta)
 
 if (bandw == TRUE) h = sd(as.matrix(A)%*%beta)*((0.75*N)^(-1/(D+4)))
 
 y_est = rep(0, S)
 for (s in 1:S){

  K_u = 1
  for (d in 1:D){
   # diff = as.numeric(Asamp[s, ]^d%*%beta[, d]) - A^d%*%beta[, d]  # diff: d times N
   diff = as.numeric(as.matrix(Asamp[s, ])%*%beta[, d]) - as.matrix(A)%*%beta[, d]  # diff: d times N
   
   # if (bandw == TRUE) h = sd(A%*%beta[, d])*((0.75*N)^(-1/(d+4)))
   
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
   y_est[s] = 0
   # cat("\n K_u is all zero \n" )
  }else y_est[s] = num/denom
 }
 return(y_est)
}

get_r_beta <- function(beta, A_sample, p, d, A, C, Y, Y_ac, Yrep, 
                       Sigma, alpha, fA_fAC, weight, augment, h, bandw){
  
  n = length(Y) 
  beta = matrix(beta, nrow=p)
  d = ncol(beta)

  # Get U and E[U | A, C]
  y_aTb = estimate_bart(Y, A, A, beta, h, bandw)
  
  U = list() 
  U_AC = list() 
  for(i in 1:n){
    U[[i]] = (Y[i] - y_aTb[i])*(alpha[[i]])
    U_AC[[i]] = (Y_ac[i] - y_aTb[i])*(alpha[[i]])
  }
  
  # Get Eq[ E[U|A,C] | C ] = \int_A E[U | A, C] p(A) = (1/samp_size)*
  if (augment == TRUE){
   samp_size = nrow(A_sample)
   
   alp = rep(1, samp_size)
   for(j in 1:d) alp = cbind(alp, (A_sample)^j)
   alp = alp[, -1]
   alpha_samp = as.matrix(alp)
   
   y_aTb_samp = estimate_bart(Y, A_sample, A, beta, h, bandw)
   
   for (i in 1:d) y_aTb_samp = y_aTb_samp + as.matrix(A_sample)^i%*%beta[, i]
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
          s = fA_fAC[i]*(U[[i]] - U_AC[[i]]) + q_U_AC[[i]]
       }
    }
    temp = s + temp
  }
  r_beta = temp/n
  r_beta = as.vector(r_beta)
  return(r_beta)
}


get_first_derivatives <- function(beta, A_sample, p, d, A, C, Y, Y_ac, Yrep, 
                                  Sigma, alpha, fA_fAC, delta, weight, augment, h, bandw){
  beta = as.vector(beta)
  K = length(beta)
  
  deriv = c()
  for (k in 1:K){
    e = matrix(rep(0, K), ncol = d)
    e[k] = 1
    
    bk_plus = beta + delta*e 
    bk_minus = beta - delta*e 
    
    r_bk_plus = get_r_beta(bk_plus, A_sample, p, d, A, C, Y, Y_ac, Yrep, 
                           Sigma, alpha, fA_fAC, weight, augment, h, bandw)
    r_bk_minus = get_r_beta(bk_minus, A_sample, p, d, A, C, Y, Y_ac, Yrep, 
                            Sigma, alpha, fA_fAC, weight, augment, h, bandw)s
    first_d = (sum((r_bk_plus)^2) - sum((r_bk_minus)^2))/(2*delta)
    deriv = c(deriv, first_d)
  }
  return(deriv)
}


get_second_derivatives <- function(beta, A_sample, p, d, A, C, Y, Y_ac, Yrep, 
                                   Sigma, alpha, fA_fAC, delta, weight, augment, h, bandw){
  beta = as.vector(beta) 
  K = length(beta)
  
  second_deriv = matrix(rep(0, K^2), nrow = K)
  for (k in 1:K){
    e = matrix(rep(0, K), ncol = d)
    e[k] = 1
    
    bk_plus = beta + delta*e
    bk_minus = beta - delta*e 
    
    r_bk_plus = get_first_derivatives(bk_plus, A_sample, p, d, A, C, Y, Y_ac, Yrep, 
                                      Sigma, alpha, fA_fAC, delta, weight, augment, h, bandw)
    r_bk_minus = get_first_derivatives(bk_minus, A_sample, p, d, A, C, Y, Y_ac, Yrep, 
                                       Sigma, alpha, fA_fAC, delta, weight, augment, h, bandw)
    
    second_deriv[k,] = {r_bk_plus - r_bk_minus}/(2*delta)
  }
  return(second_deriv)
}

update_beta <- function(beta, first_deriv, second_deriv, mu, lambda){
  p = length(first_deriv)
  second_deriv = second_deriv + mu*diag(p)
  vec_beta_new= as.vector(beta) - lambda*solve(second_deriv, tol = 1e-20)%*%first_deriv
  beta_new = matrix(vec_beta_new, nrow = nrow(beta))
  return(beta_new)
}


