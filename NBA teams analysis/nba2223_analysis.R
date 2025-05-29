library(MASS)

load("nba22_23.RData")

## bootstrap sample size
B = 200
hat_Mu <- matrix(NA, nrow = B, ncol = 30)
hat_Muwo <- matrix(NA, nrow = B, ncol = 30)
  
for(b in 1:B) {
  allgames1 = allgames[sample(nrow(allgames), size=1320, replace=TRUE), ]
  ## outcome vector
  YY = rep(0, 1320)
  ## 
  for (k in 1:1320) {
    YY[k] = allgames1[k,2]-allgames1[k,4]
  }
  ## the M matrix
  M_mat <- matrix(0, nrow = 1320, ncol = 30)
  ## covariate matrix X
  XX <- matrix(0, nrow = 1320, ncol = 3)
  ##
  for (i in 1:30) {
    for (j in 1:30) {
      for (k in 1:1320) {
        if(allgames1[k,1]==tnames[i,2] && allgames1[k,3]==tnames[j,2]) 
        {
          M_mat[k,i]=1; M_mat[k,j]=-1;
          XX[k,1] = -1
          XX[k,2] = (tnames[i,3]-tnames[j,3])
          XX[k,3] = (tnames[i,4]-tnames[j,4])
        }
      }
    }
  }
  ## lap_mat = t(M_mat)%*%M_mat ## alternate way to obtain laplacian matrix
  ## LSE of merit parameter
  hat_Mu[b,] <- (ginv(t(cbind(M_mat, XX))%*%cbind(M_mat, XX))%*%t(cbind(M_mat, XX))%*%YY)[1:30]
  hat_Muwo[b,] <- (ginv(t(M_mat)%*%M_mat)%*%t(M_mat)%*%YY)[1:30]
}


#rank(-hat_Mu[1,])
hat_Mu200 <- hat_Mu
hat_Muwo200 <- hat_Muwo
save(hat_Mu200, hat_Muwo200, file = "hat_mu200.RData")
