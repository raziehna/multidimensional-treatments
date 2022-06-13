
get_regression <- function(C, A, Y, dat, method, samp_size, A_sample){
 
 n = dim(C)[1]
 p = dim(A)[2]
 
 Crep = C[rep(seq_len(n), each=samp_size), ]
 Arep = coredata(A_sample)[rep(seq(nrow(A_sample)), n), ]
 
 rownames(Crep) = 1:nrow(Crep)
 rownames(Arep) = 1:nrow(Arep)
 
 # (1) Kernels +++++++++++++++++++++++++++
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
  Y_ac = fitted(model.np)
  
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
  Y_ac = mod_Y$yhat.test.mean
 
  cat(head(Y_ac))
  
  # Estimate E[Y | Asamp, C]
  mod_Yrep = bart(cbind(C, A), as.double(Y), cbind(Crep, Arep), ndpost=100, nskip=100, verbose=FALSE)
  Yrep = mod_Yrep$yhat.test.mean
  
 }else{ # method = GLM
  Csum = rowSums(C - colMeans(C))
  AC = rowSums(A)*Csum
  A2 = A^2
  colnames(A2)= paste0("Asq", 1:p)
  datY = data.frame(A, A2, Csum, AC, Y)
  fmla_Y = as.formula(paste0(paste0("Y~-1+"), paste0("A", 1:p, collapse = "+"), "+", 
                             paste0(colnames(A2), collapse = "+"), "+", "Csum + AC", collapse = "+"))
  glm_y = glm(fmla_Y, datY, family="gaussian")
  beta_y = glm_y$coefficients
  
  # E[Y | A, C] 
  Y_ac = as.matrix(A)%*%beta_y[1:p] + as.matrix(A^2)%*%beta_y[(p+1):(2*p)] + beta_y[2*p+1]*Csum + beta_y[2*p+2]*rowSums(A)*Csum
  
  cat(head(Y_ac))
  
  # Estimate E[Y | Asamp, C] 
  Yrep = as.matrix(Arep)%*%beta_y[1:p] + as.matrix(Arep^2)%*%beta_y[(p+1):(2*p)] + beta_y[2*p+1]*Csum + beta_y[2*p+2]*rowSums(Arep)*Csum
 }
 
 return(list(Y_ac = Y_ac, 
             Yrep = Yrep))
}
 
 