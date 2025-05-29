library(MASS)  
library(VGAM)  

# incidence matrix of complete graph
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

# Generate full-length mu for max sample size
generate_mu_full <- function(m_max, scl) {
  mu <- numeric(28 * m_max)
  idx <- 1
  for (b in 1:7) {
    for (a in 1:(8 - b)) {
      mu[idx:(idx + m_max - 1)] <- rep(a * 2, m_max)
      idx <- idx + m_max
    }
  }
  return(mu / scl)
}

# Core simulation with reused largest sample and rank comparison
simulate_rank_accuracy <- function(M, mm, rplct, scl) {
  m_max <- max(mm)
  n_edges <- nrow(M)
  distance_mat <- matrix(NA, nrow = rplct, ncol = length(mm))
  true_rank <- 8:1
  
  for (j in 1:rplct) {
    mu_full <- generate_mu_full(m_max, scl)
    Y_full <- rnorm(length(mu_full), mean = mu_full, sd = 1)
    
    for (i in seq_along(mm)) {
      m <- mm[i]
      
      # Select the first m values in each of the 28 blocks
      idx <- rep(1:28, each = m)
      sel <- (idx - 1) * m_max + rep(1:m, times = 28)
      YY <- Y_full[sel]
      
      # Replicate incidence matrix
      MM <- M[rep(1:n_edges, each = m), ]
      
      # Estimation
      S <- crossprod(MM, YY)
      lap_mat <- crossprod(MM)
      mu_hat <- ginv(lap_mat) %*% S
      rank_hat <- rank(mu_hat, ties.method = "first")
      
      distance_mat[j, i] <- identical(as.integer(rank_hat), true_rank)
    }
  }
  return(colMeans(distance_mat))
}

# Parameters
mm <- seq(10, 1000, by = 10)
rplct <- 1000
SCL <- c(1, 10, 100)
M <- complete_incidence_mat(8)

# Output matrix
avg_distance_complete_normal <- matrix(NA, nrow = length(SCL), ncol = length(mm))

# Run simulation for each scale
for (vv in seq_along(SCL)) {
  scl <- SCL[vv]
  avg_distance_complete_normal[vv, ] <- simulate_rank_accuracy(M, mm, rplct, scl)
}

##########################################
## save data
##########################################
save(mm, avg_distance_complete_normal, file = "example1_rank.RData")
#########################################
## plot
########################################
## normal rvs observations
pdf(file = "corrrect_ranking_without_cov.PDF")
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


