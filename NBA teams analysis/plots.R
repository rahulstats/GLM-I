library(tidyverse)
library(lme4)
## with covariates
names_team <- tnames[,2]
df <- t(as.matrix(hat_Mu200))
df_rank <- matrix(NA, nrow = 30, ncol = 200)
for (i in 1:200) {
  df_rank[,i]= rank(-df[,i])
}

pdf(file = "nba22_23withcov.pdf")
par(mar=c(3,10.2,2,0))
boxplot(t(df_rank), horizontal=TRUE, outline=FALSE, yaxt='n')
axis(2, at=1:30, labels=names_team, las=2, cex=0.5)
title(xlab = "Rank", line = 2)            # Add x-axis text
dev.off()

## without covariates
df1 <- t(as.matrix(hat_Muwo200))
df1_rank <- matrix(NA, nrow = 30, ncol = 200)
for (i in 1:200) {
  df1_rank[,i]= rank(-df1[,i])
}

pdf(file = "nba22_23withoutcov.pdf")
par(mar=c(3,10.2,2,0))
boxplot(t(df1_rank), horizontal=TRUE,outline=FALSE, yaxt='n')
axis(2, at=1:30, labels=names_team, las=2, cex=0.5)
title(xlab = "Rank", line = 2)            # Add x-axis text
dev.off()
