pdf(file = "example3_large_dim.PDF", width = 7, height = 6)
# plot margin
par(mar=c(5,7,4,1)+.1)
# 
plot(kk, distance_complete_avg, type = 'l', pch = 18, col = 'black', 
     ylim = c(0.15, 0.4),
       #c(0, max(c(distance_complete_avg,distance_ER1_avg,distance_ER2_avg,distance_ER3_avg))), 
     ylab=expression(paste("Maximum difference between LSE and ",  mu[i])), xlab="K", 
     cex.main=1.25, cex.lab=1.3, cex.axis=1.1, lty = 2)

# Add a line for path graph
lines(kk, distance_ER1_avg, type = 'l', pch = 18, col = "blue", 
      lty=2)

# Add a line for star graph
lines(kk, distance_ER2_avg, type = 'l', pch = 18, col = "red", 
      lty=2)

# Add a line for star graph
lines(kk, distance_ER3_avg, type = 'l', pch = 18, col = "green", 
      lty=2)

# Add legends
legend("topright", legend=c("Complete", expression(paste(p, " = ", 0.5)), expression(paste(p, " = ", log(K)^3/K)), expression(paste(p, " = ", sqrt(log(K)^3/K)))),
       col=c("black", "blue", "red", "green"), lty=2, cex= 1)
#,title="Graphs", bty = "n") # remove box
dev.off()
