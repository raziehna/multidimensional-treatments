
get_regression <- function(n, p, C, A, Y, dat, method, samp_size, A_sample){
 
 Crep = matrix(rep(C, each=samp_size), nrow=n*samp_size)
 Arep = coredata(A_sample)[rep(seq(nrow(A_sample)), n), ]
 
 # (1) GLM +++++++++++++++++++++++++++
 if (method == "GLM"){
  
  Csum = rowSums(C)
  # AC = rowSums(A)*Csum
  AC = matrix(0, n, p)
  for (i in 1:p) AC[, i] = A[, i]*Csum
  colnames(AC) = paste0("CA", 1:p)
  A2 = A^2
  colnames(A2)= paste0("Asq", 1:p)
  datY = data.frame(A, A2, C, Csum, AC, Y)
  fmla_Y = as.formula(paste0(paste0("Y~-1+"), paste0("A", 1:p, collapse = "+"), "+", 
                             paste0(colnames(A2), collapse = "+"), "+ C1 + C2 + C3 + C4 +",
                             paste0(colnames(AC), collapse = "+"), collapse = "+"))
  glm_y = glm(fmla_Y, datY, family="gaussian")
  beta_y = glm_y$coefficients
  if (sum(is.na(beta_y)) != 0){
   beta_y[5] = beta_y[5]/2
   beta_y[6] = beta_y[6]/2
   if (p == 6){
     beta_y[11] = beta_y[5]
     beta_y[12] = beta_y[6]
   }else{
     beta_y[17] = beta_y[5]
     beta_y[18] = beta_y[6]
    }
   }
  
  # E[Y | A, C]
  y_ac = A%*%beta_y[1:p] + A^2%*%beta_y[(p+1):(2*p)] + C%*%beta_y[(2*p+1):(2*p+4)] + A%*%beta_y[(2*p+5):length(beta_y)]*Csum
  
  # E[Y | Asamp, C]
  Yrep = Arep%*%beta_y[1:p] + Arep^2%*%beta_y[(p+1):(2*p)] + Crep%*%beta_y[(2*p+1):(2*p+4)] + Arep%*%beta_y[(2*p+5):length(beta_y)]*rowSums(Crep)
  
 }else{ # method == bart
  # (2) bart ++++++++++++++++++++++++++
  # Estiamte E[Y | A, C]
  mod_Y = bart(cbind(C, A), as.double(Y), cbind(C, A), ndpost=500, nskip=100, verbose=FALSE)
  y_ac = mod_Y$yhat.test.mean

  # Estimate E[Y | Asamp, C]
  mod_Yrep = bart(cbind(C, A), as.double(Y), cbind(Crep, Arep), ndpost=500, nskip=100, verbose=FALSE)
  Yrep = mod_Yrep$yhat.test.mean
 }
 
 return(list(y_ac = y_ac, 
             Yrep = Yrep))
}
 
 