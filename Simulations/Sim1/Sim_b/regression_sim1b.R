
get_regression <- function(n, p, C, A, Y, dat, method, samp_size, A_sample, true_beta){
 
 Crep = matrix(rep(C, each=samp_size), nrow=n*samp_size)
 Arep = coredata(A_sample)[rep(seq(nrow(A_sample)), n), ]
 
 # (1) GLM +++++++++++++++++++++++++++
 if (method == "Kernels"){
  AC = rowSums(A)*rowSums(C)
  Csum = rowSums(C)
  A2 = A^2
  colnames(A2)= paste0("Asq", 1:p)
  datY = data.frame(A, A2, Csum, AC, Y)
  fmla_Y = as.formula(paste0(paste0("Y~-1+"), paste0("A", 1:p, collapse = "+"), "+", 
                             paste0(colnames(A2), collapse = "+"), "+", "Csum + AC", collapse = "+"))

  # E[Y | A, C] - kernels (ckertype = {epanechnikov vs gaussian})
  model.np <- npreg(fmla_Y, data=datY, bwmethod="cv.aic", regtype="ll", ckertype="epanechnikov")
  model.np$R2
  y_ac = fitted(model.np)
  
  # Estimate E[Y | Asamp, C] - kernels
  AC_rep = rowSums(Arep)*rowSums(Crep)
  Csum_rep = rowSums(Crep)
  Arep_2 = Arep^2
  colnames(Arep_2)= paste0("Asq", 1:p)
  newdatY = data.frame(Arep, Arep_2, Csum=Csum_rep, AC=AC_rep)
  Yrep = predict(model.np, newdata = newdatY)
  
  
 }else if (method == "bart"){ # method == bart
  # (2) bart ++++++++++++++++++++++++++
  # Estiamte E[Y | A, C]
  mod_Y = bart(cbind(C, A), as.double(Y), cbind(C, A), ndpost=100, nskip=100, verbose=FALSE)
  y_ac = mod_Y$yhat.test.mean

  # Estimate E[Y | Asamp, C]
  mod_Yrep = bart(cbind(C, A), as.double(Y), cbind(Crep, Arep), ndpost=100, nskip=100, verbose=FALSE)
  Yrep = mod_Yrep$yhat.test.mean
  
  
 }else{ # method = GLM
  AC = rowSums(A)*rowSums(C)
  Csum = rowSums(C)
  A2 = A^2
  colnames(A2)= paste0("Asq", 1:p)
  datY = data.frame(A, A2, Csum, AC, Y)
  fmla_Y = as.formula(paste0(paste0("Y~-1+"), paste0("A", 1:p, collapse = "+"), "+", 
                             paste0(colnames(A2), collapse = "+"), "+", "Csum + AC", collapse = "+"))
  glm_y = glm(fmla_Y, datY, family="gaussian")
  beta_y = glm_y$coefficients
  # E[Y | A, C]
  y_ac = A%*%beta_y[1:p] + A^2%*%beta_y[(p+1):(2*p)] + beta_y[2*p+1]*rowSums(C) + beta_y[2*p+2]*rowSums(A)*rowSums(C)

  # Estimate E[Y | Asamp, C]
  Yrep = Arep%*%beta_y[1:p] + Arep^2%*%beta_y[(p+1):(2*p)] + beta_y[2*p+1]*rowSums(Crep) + beta_y[2*p+2]*rowSums(Arep)*rowSums(Crep)


  # beta_y = true_beta
  
  # E[Y | A, C] 
  # y_ac = A%*%beta_y[1:p] + (A%*%beta_y[(p+1):(2*p)])^2 + rowSums(C) + rowSums(A)*rowSums(C)
  
  # Estimate E[Y | Asamp, C] 
  # Yrep = Arep%*%beta_y[1:p] + (Arep%*%beta_y[(p+1):(2*p)])^2 + rowSums(Crep) + rowSums(Arep)*rowSums(Crep)
  
 }
 
 return(list(y_ac = y_ac, 
             Yrep = Yrep))
}
 
 