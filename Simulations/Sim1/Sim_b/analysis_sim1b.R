
run_analysis <- function(beta, A_sample, threshold, delta, mu, lambda, A, C, Y, 
                         y_ac, Yrep, Sigma, alpha, fA_fAC, weight, augment, h, bandw){

  p = nrow(beta)
  d = ncol(beta)
  
  count = 1
  # start <- proc.time()
  repeat{
    cat("\n repeat # ", count, "\n") 
   
    # Get derivatives of r(beta) 
    first_deriv = get_first_derivatives(beta, A_sample, p, d, A, C, Y, y_ac, Yrep, 
                                        Sigma, alpha, fA_fAC, delta, weight, augment, h, bandw)
    second_deriv = get_second_derivatives(beta, A_sample, p, d, A, C, Y, y_ac, Yrep, 
                                          Sigma, alpha, fA_fAC, delta, weight, augment, h, bandw)
    
    # Update beta
    # if (count > 5) lambda = 0.1
    
    # step_size =  lambda/(count)^(1/4)
    step_size =  lambda/(count)
    # step_size = 0.01
    (beta_new = update_beta(beta, first_deriv, second_deriv, mu, step_size))
    
    # Stopping criteria
    # (1) beta_new - beta

    # beta_frob = beta%*%solve( t(beta)%*%beta )%*%t(beta)
    # beta_new_frob = beta_new%*%solve( t(beta_new)%*%beta_new )%*%t(beta_new)
    # diff_parameters = norm(beta_frob - beta_new_frob, type="f")
    diff_parameters = sqrt(sum((beta_new - beta)^2))
    cat("diff_parameters: ", diff_parameters, " \n")
    
    # (2) first_deriv < eps
    first_deriv_new = get_first_derivatives(beta_new, A_sample, p, d, A, C, Y, y_ac, Yrep,
                                            Sigma, alpha, fA_fAC, delta, weight, augment, h, bandw)
    first_deriv_norm = round(sqrt(sum(first_deriv_new^2)), 3)
    cat("first_deriv_norm: ", first_deriv_norm, " \n ")

    beta = beta_new
    
    if (diff_parameters <= threshold | count > 20) break
    # if ((first_deriv_norm <= threshold) && (diff_parameters <= threshold)) break
    count = count + 1
  }
  # end <- proc.time()
  # (end - start)
  
  return(list(est_beta=beta_new))
}



