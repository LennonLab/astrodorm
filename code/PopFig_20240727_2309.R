##An individual based model to examine dormancy, activity, death, birth

rm(list=ls())
getwd()
setwd('C:/Current/Dine/Research/AstroDorm/AstroDormCode')

AS <- read.csv("AverageShielding.csv")
s
sAS <- AS[c(1,100,200,300,400,500,600,700,800,900,1000),]
x <- c(1,100,200,300,400,500,600,700,800,900,1000)

plot(x, sAS[,4], type = "b", pch = 19, col = rgb(0,1,0, alpha = 0.6), 
     ylim=c(0,400), cex = 2, ylab = "Average Total Population (A + Ad)",
     xlab = "Time Step")
segments(x,sAS[,4]-sAS[,5],x,sAS[,4]+sAS[,5], col = rgb(0,1,0, alpha = 0.6))

points(x, sAS[,2], type = "b", pch = 19, col = rgb(0,0,0, alpha = 0.6), cex = 2)
segments(x,sAS[,2]-sAS[,3],x,sAS[,2]+sAS[,3], col = rgb(0,0,0, alpha = 0.6))

points(x, sAS[,6], type = "b", pch = 19, col = rgb(1,0,1, alpha = 0.6), cex = 2)
segments(x,sAS[,6]-sAS[,7],x,sAS[,6]+sAS[,7], col = rgb(1,0,1, alpha = 0.6))

legend(400, 475, title = "Decay in Dormancy",
       legend=c("1 %", "2 %", "3 %"),
       pch = 19,
       col = c(rgb(0,0,0, alpha = 0.6), 
               rgb(0,1,0, alpha = 0.6), rgb(1,0,1, alpha = 0.6)),
       cex = 1.1,
       box.lty = 0,
       bg= rgb(1,1,1,0),
       y.intersp= 0.5)
