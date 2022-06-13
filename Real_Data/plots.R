
library(lattice)
library(latticeExtra)

# setwd("?")

d = 1
p = 5

# Get beta
est_beta_reg = as.matrix(read.csv(paste0("run 2/beta_IPW_d", d, "_p", p, ".csv")))

# Get A_low
A_sort = A[order(A[, 1]), ]
A_low = as.matrix(A)%*%est_beta_reg
A_low_sorted = as.matrix(A_sort)%*%est_beta_reg

# Get E[Y | A_low]
mod_Yq_reg = bart(A_low, as.double(Y), A_low_sorted, ndpost=100, nskip=100, verbose=FALSE)
Yq_reg = mod_Yq_reg$yhat.test.mean

# Get the heatmap
eff_reg = matrix(0, n, n)
for (i in 1:n) eff_reg[, i] = Yq_reg - as.numeric(Yq_reg[i])  

# Plot the heatmap
new.palette = colorRampPalette(c("red","white","blue"),space="rgb")

levelplot(eff_reg,
          panel = function(...){
           panel.levelplot(...)
           panel.abline(0, 1, col="white")
          },
          col.regions=new.palette(20),
          main = paste0("p = ?"), xlab = "", ylab = "", 
          at=seq(-10, 10, 2)
)