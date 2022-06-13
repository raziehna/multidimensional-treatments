
library(lattice)

# p = 6
Yq_sorted = A%*%true_beta[, 1] + A^2%*%true_beta[, 2]


# # AIPW 
# bandw = F
# h = 2
# est_beta_aipw = matrix(c(0.5388236, 0.7856611, 0.4609090, 0.9159269, 0.5152465, 
#                          0.4189074, -0.5388875, -0.7857738, -0.4608907, -0.9155873, 
#                          -0.5151376, -0.4190317), byrow = F, p)
Yq_aipw = estimate_kernel(Y, A, A, est_beta, h, bandw) 


eff_0 = matrix(0, n, n)
eff_aipw = matrix(0, n, n)
for (i in 1:n) eff_0[, i] = Yq_sorted - as.numeric(Yq_sorted[i])  
for (i in 1:n) eff_aipw[, i] = Yq_aipw - as.numeric(Yq_aipw[i]) 


new.palette=colorRampPalette(c("black","white","black"),space="rgb")

levelplot(eff_0, col.regions=new.palette(20), 
          main = paste0("True Effects"), 
          at=seq(-100, 100, 10))

levelplot(eff_aipw, col.regions=new.palette(20), 
          main = paste0("AIPW Effects"), 
          at=seq(-100, 100, 10))









