library(MASS)

## incident matrix 
ER_incidence_mat <- function(k,p)
{
  M <- matrix(0, nrow = 1, ncol = k)
  for (b in 1:(k-1)) {
    bb = k-b
    M1 <- matrix(0, nrow = bb, ncol = k)
    ## 
    for (i in 1:bb) {
      M1[i,b] <- rbinom(1,1,p)
      M1[i,i+b] <- -M1[i,b]
    }
    M <- rbind(M, M1)
  }
  return(M[-1,])
}

# implementation
kk <- seq(100, 1000, by = 100)
rplct = 1000
distance_ER = matrix(NA, nrow =rplct, ncol = length(kk))

for (j in 1:rplct) {
  
  for (i in 1:length(kk)) {
    k = kk[i]
    ## defining mean vector
    mu = 0
    for (b in 1:(k-1)) {
      mu = c(mu, 2*c(1:(k-b)))
    }
    ## mu vector
    mu = mu[-1]
    ## data vector
    YY = rnorm(k*(k-1)/2, mean = mu, sd = 1)
    # M matrix
    MM = ER_incidence_mat(k,(log(k))^3/k)
    #### S vector
    S = t(MM)%*%YY
    ## lap_mat
    lap_mat <- t(MM)%*%MM # Laplacian matrix of Erdos Reyni graph
    ## estimator
    mu_hat <- ginv(lap_mat)%*%S
    ## distance
    distance_ER[j,i] = max(abs((mu_hat+c(seq(-k+1, k-1, by = 2)))))
  }
}

#distance_complete
distance_ER_avg = colMeans(distance_ER)

# save data
save(distance_ER, distance_ER_avg,kk, file = "inf_norm_ER_logk3_normal_high_dim_200_rep.RData")

