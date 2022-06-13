
run_analysis <- function(beta, A_sample, threshold, delta, mu, lambda, A, C, Y, 
                         Y_ac, Yrep, Sigma, alpha, fA_fAC, weight, augment, h, bandw){

 # n = nrow(dat)
 p = nrow(beta)
 d = ncol(beta)
 
 collection_beta = as.vector(beta) # collecting all updates of beta
 
 count = 1
 repeat{
  cat("\n repeat # ", count, "\n") 
  
  first_deriv = get_first_derivatives(beta, A_sample, p, d, A, C, Y, Y_ac, Yrep, 
                                      Sigma, alpha, fA_fAC, delta, weight, augment, h, bandw)
  second_deriv = get_second_derivatives(beta, A_sample, p, d, A, C, Y, Y_ac, Yrep, 
                                        Sigma, alpha, fA_fAC, delta, weight, augment, h, bandw)
  
  step_size =  lambda/count
  beta_new = update_beta(beta, first_deriv, second_deriv, mu, step_size)
  
  diff_parameters = sqrt(sum((beta_new - beta)^2))
  cat("diff_parameters: ", diff_parameters, " \n")
  
  beta = beta_new
  collection_beta = cbind(collection_beta, as.vector(beta))
  if (diff_parameters <= threshold | count > 20) break
  count = count + 1
 }
 return(list(est_beta=beta_new, 
             count=count))
}


