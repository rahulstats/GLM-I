library(VGAM)
library(MASS)

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
# Precompute once
M <- complete_incidence_mat(8)
true_mu <- c(7, 5, 3, 1, -1, -3, -5, -7)
n_nodes <- length(true_mu)
scl <- 1
true_mu_scaled <- true_mu / scl

mm <- seq(10, 1000, by = 5)
rplct <- 1000
m_max <- max(mm)
n_edges <- nrow(M)
n_mu_blocks <- 28  

# Generate full-length mu for max sample size
generate_mu_full <- function(m_max) {
  mu <- vector("numeric", length = n_mu_blocks * m_max)
  idx <- 1
  for (b in 1:7) {
    for (a in 1:(8 - b)) {
      mu[idx:(idx + m_max - 1)] <- rep(a * 2, m_max)
      idx <- idx + m_max
    }
  }
  return(mu / scl)
}

# Core simulation function with reused samples
simulate_distances_fixed_sample <- function(error_func) {
  distance_mat <- matrix(NA, nrow = rplct, ncol = length(mm))
  
  for (j in 1:rplct) {
    # Full length sample for max m
    mu_full <- generate_mu_full(m_max)
    Y_full <- error_func(length(mu_full)) + mu_full
    
    for (i in seq_along(mm)) {
      m <- mm[i]
      
      # Subset data
      idx <- rep(1:(n_mu_blocks), each = m)
      sel <- (idx - 1) * m_max + rep(1:m, times = n_mu_blocks)
      YY <- Y_full[sel]
      
      # Corresponding incidence matrix
      MM <- M[rep(1:n_edges, each = m), ]
      S <- crossprod(MM, YY)
      lap_mat <- crossprod(MM)
      mu_hat <- ginv(lap_mat) %*% S
      
      distance_mat[j, i] <- sqrt(sum((mu_hat - true_mu_scaled)^2))
    }
  }
  return(colMeans(distance_mat))
}

# Run simulations
avg_distance_normal <- simulate_distances_fixed_sample(function(n) rnorm(n, mean = 0, sd = 1))
avg_distance_t3 <- simulate_distances_fixed_sample(function(n) rt(n, df = 3) / sqrt(3))
avg_distance_t2 <- simulate_distances_fixed_sample(function(n) rt(n, df = 2))

##########################################
## save data
##########################################
save(mm, avg_distance_normal, avg_distance_t3, avg_distance_t2, file = "example1.RData")
#########################################
## plot
########################################
## normal rvs observations
pdf(file = "mse_complete_without_cov.PDF") ## output pdf declaration
# plot against m
plot(mm, avg_distance_normal, type = 'l', pch = 18, col = 'black', 
     ylim = c(0, max(c(avg_distance_normal, avg_distance_t2, avg_distance_t3))), 
     ylab="Estimated vectorised RMSE", xlab="m", 
     cex.main=1.25, cex.lab=1.5, cex.axis=1.1, lty = 2)

# Add a line for path graph
lines(mm, avg_distance_t3, type = 'l', pch = 18, col = "blue", 
      lty=2)

# Add a line for star graph
lines(mm, avg_distance_t2, type = 'l', pch = 18, col = "red", 
      lty=2)

# Add legends
legend("topright", legend=c("Normal", "t(3)", "t(2)"),
       col=c("black", "blue", "red"), lty=2, cex= 1)
#,title="Graphs", bty = "n") # remove box
dev.off()


