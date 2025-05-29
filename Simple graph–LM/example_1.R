library(VGAM)
library(MASS)
scl =1

# incident matrix of complete graph
complete_incidence_mat <- function(k) {
  num_edges <- k * (k - 1) / 2
  M <- matrix(0, nrow = num_edges, ncol = k)
  edge <- 1
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      M[edge, i] <- 1
      M[edge, j] <- -1
      edge <- edge + 1
    }
  }
  return(M)
}
##
M = complete_incidence_mat(8)
##
mm = seq (10,1000, by=5)
rplct = 1000

# normal distributed errors

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
    ## data vector
    YY = rnorm(28*m, mu, 1)
    # M matrix
    MM = M[rep(seq_along(rep(m, 28)), rep(m, 28)), ]
    #### S vector
    S = t(MM)%*%YY
    ## lap_mat
    lap_mat <- t(MM)%*%MM #path_cauchy_mat(m)
    ## estimator
    mu_hat <- ginv(lap_mat)%*%S
    #print(mu_hat)
    ## distance
    distance_complete[j,i] = sqrt(sum((mu_hat-(c(7,5,3,1,-1,-3,-5,-7)/scl))^2))
  }
}
## 
#distance_complete
avg_distance_complete_normal = colMeans(distance_complete)


# t3 distributed errors

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
    ## data vector
    YY = rt(28*m, df= 3)/sqrt(3) + mu 
    #YY = rcauchy(28*m, mu, 1)
    # M matrix
    MM = M[rep(seq_along(rep(m, 28)), rep(m, 28)), ]
    #### S vector
    S = t(MM)%*%YY
    ## lap_mat
    lap_mat <- t(MM)%*%MM #path_cauchy_mat(m)
    ## estimator
    mu_hat <- ginv(lap_mat)%*%S
    #print(mu_hat)
    ## distance
    distance_complete[j,i] = sqrt(sum((mu_hat-(c(7,5,3,1,-1,-3,-5,-7)/scl))^2))
  }
}
## 
#distance_complete
avg_distance_complete_t3dist = colMeans(distance_complete)



## t with dof 2 errors

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
    ## data vector
    YY = rt(28*m, df= 2) + mu 
    #YY = rcauchy(28*m, mu, 1)
    # M matrix
    MM = M[rep(seq_along(rep(m, 28)), rep(m, 28)), ]
    #### S vector
    S = t(MM)%*%YY
    ## lap_mat
    lap_mat <- t(MM)%*%MM #path_cauchy_mat(m)
    ## estimator
    mu_hat <- ginv(lap_mat)%*%S
    #print(mu_hat)
    ## distance
    distance_complete[j,i] = sqrt(sum((mu_hat-(c(7,5,3,1,-1,-3,-5,-7)/scl))^2))
  }
}
## 
#distance_complete
avg_distance_complete_t2dist = colMeans(distance_complete)

# save data
#save(mm, avg_distance_complete_normal, avg_distance_complete_t3dist, avg_distance_complete_t2dist, file = "example1.RData")

## plot
pdf(file = "mse_complete_without_cov.PDF")
# plot against m
plot(mm, avg_distance_complete_normal, type = 'l', pch = 18, col = 'black', 
     ylim = c(0, max(c(avg_distance_complete_normal, avg_distance_complete_t2dist, avg_distance_complete_t3dist))), 
     ylab="Estimated vectorised RMSE", xlab="m", 
     cex.main=1.25, cex.lab=1.5, cex.axis=1.1, lty = 2)

# Add a line for path graph
lines(mm, avg_distance_complete_t3dist, type = 'l', pch = 18, col = "blue", 
      lty=2)

# Add a line for star graph
lines(mm, avg_distance_complete_t2dist, type = 'l', pch = 18, col = "red", 
      lty=2)

# Add legends
legend("topright", legend=c("Normal", "t(3)", "t(2)"),
       col=c("black", "blue", "red"), lty=2, cex= 1)
#,title="Graphs", bty = "n") # remove box
dev.off()

