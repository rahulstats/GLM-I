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

## inverse of laplacian
complete_lap_matrix_inv <- function(k)
{
  matt = matrix(rep(-1,k*k), nrow = k, ncol = k)
  diag(matt) = rep(k-1, k)
  return(matt/k^2)
}

# implementation
kk <- seq(100, 1000, by = 100)
rplct =1000
distance_complete = matrix(NA, nrow =rplct, ncol = length(kk))

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
    MM = complete_incidence_mat(k)
    #### S vector
    S = t(MM)%*%YY
    ## lap_mat
    #lap_mat <- t(MM)%*%MM #path_lap_mat(m)
    ## estimator
    mu_hat <- complete_lap_matrix_inv(k)%*%S
    ## distance
    distance_complete[j,i] = max(abs((mu_hat+c(seq(-k+1, k-1, by = 2)))))
  }
}

#distance_complete
distance_complete_avg = colMeans(distance_complete)

# save data
save(distance_complete, distance_complete_avg,kk, file = "inf_norm_complete_normal_high_dim_200_rep.RData")

