library(MASS)

load("nba22_23.RData")

## outcome vector
YY = rep(0, 1320)
## 
for (k in 1:1320) {
  YY[k] = allgames[k,2]-allgames[k,4]
}
## the M matrix
M_mat <- matrix(0, nrow = 1320, ncol = 30)
## covariate matrix X
XX <- matrix(0, nrow = 1320, ncol = 3)
##
for (i in 1:30) {
  for (j in 1:30) {
    for (k in 1:1320) {
      if(allgames[k,1]==tnames[i,2] && allgames[k,3]==tnames[j,2]) 
      {
        M_mat[k,i]=1; M_mat[k,j]=-1;
        XX[k,1] = -1
        XX[k,2] = (tnames[i,3]-tnames[j,3])
        XX[k,3] = (tnames[i,4]-tnames[j,4])
      }
    }
  }
}
## LSE of merit parameter
# hat_Mu <- (ginv(t(cbind(M_mat, XX))%*%cbind(M_mat, XX))%*%t(cbind(M_mat, XX))%*%YY)[1:30]
# hat_Muwo <- (ginv(t(M_mat)%*%M_mat)%*%t(M_mat)%*%YY)
##
df = as.data.frame(cbind(YY, M_mat, XX))
df1 = as.data.frame(cbind(YY, M_mat))
##
library(tidyverse)
##
model <- lm(YY ~., data = df)
summary(model)
##
model1 <- lm(YY ~., data = df1)
summary(model1)
##
((t(M_mat[1:600,])%*%XX[1:600,])/600)
##########################
## rank of teams using different models
#############################

rank1=rank(-((ginv(t(cbind(M_mat, XX))%*%cbind(M_mat, XX))%*%t(cbind(M_mat, XX))%*%YY)[1:30]))

rank2=rank(-((ginv(t(cbind(M_mat, XX[,1:2]))%*%cbind(M_mat, XX[,1:2]))%*%t(cbind(M_mat, XX[,1:2]))%*%YY)[1:30]))

rank3=rank(-((ginv(t(cbind(M_mat, XX[,1]))%*%cbind(M_mat, XX[,1]))%*%t(cbind(M_mat, XX[,1]))%*%YY)[1:30]))

rank4=rank(-(ginv(t(M_mat)%*%M_mat)%*%t(M_mat)%*%YY))

names_team <- tnames[,2]

#ginv(t(cbind(M_mat, XX))%*%cbind(M_mat, XX))[33,33]

#pnorm(-0.02013299,mean = 0, sd=sqrt(7.29931e-06))
##. making rank table
ddf <- rbind.data.frame(names_team, rank4, rank3, rank2, rank1)

library(xtable)

xtable(ddf)

## distance between rank vectors
library(BayesMallows)
## rank matrix
rank_mat <- rbind(rank4, rank3, rank2, rank1)
## rank distance
rank_distance(rankings = rank_mat, rho = rank4, metric = "cayley")
rank_distance(rankings = rank_mat, rho = rank3, metric = "cayley")
rank_distance(rankings = rank_mat, rho = rank2, metric = "cayley")
rank_distance(rankings = rank_mat, rho = rank1, metric = "cayley")
