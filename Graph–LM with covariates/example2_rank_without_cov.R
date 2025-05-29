library(VGAM)
library(MASS)

# Generate incidence matrix for a complete graph
complete_incidence_mat <- function(k) {
  idx <- 1
  M <- matrix(0, nrow = k * (k - 1) / 2, ncol = k)
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      M[idx, i] <- 1
      M[idx, j] <- -1
      idx <- idx + 1
    }
  }
  M
}

##
M = complete_incidence_mat(8)
##
mm = seq (10,1000, by=10)
rplct = 1000

SCL = c(1,10,100)

avg_distance_complete_normal = matrix(NA, nrow = length(SCL), ncol = length(mm))

for (vv in 1:length(SCL)) {
  
  scl = SCL[vv]
  distance_complete = matrix(NA, nrow =rplct, ncol = length(mm))
  
  for (j in 1:rplct) {
    
    for (i in 1:length(mm)) {
      m = mm[i]
      ## defining mean vector
      mu = 0
      for (b in 1:7) {
        mu1=0
        for (a in 1:(8-b)) {
          mu1[((a-1)*m+1):(a*m)] = rep(a*2, m)
        }
        mu = c(mu, mu1)
      }
      ## mu vector
      mu = mu[-1]/scl
      
      ## the covariate matrix
      X = matrix(sample(c(-1,1),28*m*2,replace=TRUE) , nrow = 28*m, ncol = 2)
      ## data vector
      YY = rnorm(28*m, mu, 1) + rowSums(X)
      # M matrix
      MM = M[rep(seq_along(rep(m, 28)), rep(m, 28)), ]
      #### S vector
      S = t(MM)%*%YY
      ## lap_mat
      lap_mat <- t(MM)%*%MM 
      ## estimator
      rank_hat <- rank(ginv(lap_mat)%*%S)
      #print(rank_hat)
      ## distance
      distance_complete[j,i] = as.numeric(identical(as.integer(rank_hat), c(8:1)))
    }
  }
  ## 
  #distance_complete
  avg_distance_complete_normal[vv,] = colMeans(distance_complete)
  
}

# save data
#save(mm, avg_distance_complete_normal, file = "example2_rank_binary_cov.RData")

# plot
pdf(file = "example2_rank_binary_cov.PDF")
# plot against m
plot(mm, avg_distance_complete_normal[1,], type = 'l', pch = 18, col = 'black', 
     ylim = c(0, max(avg_distance_complete_normal)), 
     ylab="Correct ranking probability", xlab="m", 
     cex.main=1.25, cex.lab=1.5, cex.axis=1.1, lty = 2)

# Add a line for path graph
lines(mm, avg_distance_complete_normal[2,], type = 'l', pch = 18, col = "blue", 
      lty=2)

# Add a line for star graph
lines(mm, avg_distance_complete_normal[3,], type = 'l', pch = 18, col = "red", 
      lty=2)

# Add legends
legend(845,0.95, legend=c(expression(paste(gamma, " = ", 0)), expression(paste(gamma, " = ", 1)), expression(paste(gamma, " = ", 2))),
       col=c("black", "blue", "red"), lty=2, cex= 1)
#,title="Graphs", bty = "n") # remove box
dev.off()



