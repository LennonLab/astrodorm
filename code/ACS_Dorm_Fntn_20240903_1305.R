##An analytic model to assess dormancy in autocatalytic sets

rm(list=ls())
getwd()
setwd('C:/Current/Dine/Research/AstroDorm/AstroDorm2_Code')

# some notes about this code - it examines the autocatalytic set in a non
# time varying environment - creates some figures that could be modified \
# for publication

# Code to examine dormancy in an autocatalytic set

# I need some reactants m1, m2, m3, m4, Pa

# m1 + m2 -> B
# m3 + m4 -> C
# m5 + m6 -> A
# C + m5 + m6 -> A + C
# A + m1 + m2 -> B + A
# B + m3 + m4 -> C + B
# A <-> Ai

ACSDoMo <- function(ts, ic, rr, dr, sr) {
  
  #ts = time steps, ic = initial concentration, dr = decay rates, 
  #sr = supply rates

  #initialize vectors
  vm1 <- rep(NA, ts)
  vm2 <- rep(NA, ts)
  vm3 <- rep(NA, ts)
  vm4 <- rep(NA, ts)
  vB <- rep(NA, ts)
  vm5 <- rep(NA, ts)
  vm6 <- rep(NA, ts)
  vA <- rep(NA, ts)
  vC <- rep(NA, ts)
  vAi <- rep(NA, ts)

  #starting concentrations
  vm1[1] <- ic[1] 
  vm2[1] <- ic[2] 
  vm3[1] <- ic[3] 
  vm4[1] <- ic[4] 
  vm5[1] <- ic[5] 
  vm6[1] <- ic[6] 
  vA[1] <- ic[7]
  vAi[1] <- ic[8]
  vB[1] <- ic[9]
  vC[1] <- ic[10] 
  
  #reaction rates
  kB <- rr[1] #some reaction rate for m1 + m2 -> B
  kC <- rr[2] #some reaction rate for m3 + m4 -> C
  kA <- rr[3] #some reaction rate for m5 + m6 -> A
  kcA <- rr[4]  #some reaction rate for C + pa + pa -> A + C
  kcC <- rr[5]  #some reaction rate for B + m3 + m4 -> C + B
  kcB <- rr[6]  #some reaction rate for A + m1 + m2 -> B + A 
  kAi <- rr[7]  #dormancy transition for A
  kAa <- rr[8]  # activity transition for Ai
  
  #decay rates
  dB <- dr[1] #some decay rate for B
  dC <- dr[2] #some decay rate for C
  dA <- dr[3] #some decay rate for A
  dAi <- dr[4] #some decay rate for Ai
  
  #supply rates
  sm1 <- sr[1] #supply rate for m1
  sm2 <- sr[2] #supply rate for m2
  sm3 <- sr[3] #supply rate for m3
  sm4 <- sr[4] #supply rate for m4
  sm5 <- sr[5] #supply rate for m5
  sm6 <- sr[6] #supply rate for m6
  
  for (i in 2:length(vB)){   
    vB[i] <- vB[i-1] + (vm1[i-1]*vm2[i-1])*kB + 
      (vA[i-1]*vm1[i-1]*vm2[i-1])*kcB - dB*vB[i-1]
    
    vm1[i] <- vm1[i-1] - (vm1[i-1]*vm2[i-1])*kB -
      (vA[i-1]*vm1[i-1]*vm2[i-1])*kcB + sm1 #
    
    vm2[i] <- vm2[i-1] - (vm1[i-1]*vm2[i-1])*kB -
      (vA[i-1]*vm1[i-1]*vm2[i-1])*kcB +sm2 #
    
    vC[i] <-  vC[i-1] + (vm3[i-1]*vm4[i-1])*kC + 
      (vB[i-1]*vm3[i-1]*vm4[i-1])*kcC - dC*vC[i-1]
    
    vm3[i] <- vm3[i-1] - (vm3[i-1]*vm4[i-1])*kB - 
      (vB[i-1]*vm3[i-1]*vm4[i-1])*kcC + sm3 #
    
    vm4[i] <- vm4[i-1] - (vm3[i-1]*vm4[i-1])*kB - 
      (vB[i-1]*vm3[i-1]*vm4[i-1])*kcC + sm4 #
    
    vA[i] <- vA[i-1] + vAi[i-1]*kAa + (vm5[i-1]*vm6[i-1])*kA + 
      (vC[i-1]*vm5[i-1]*vm6[i-1])*kcA - vA[i-1]*kAi - dA*vA[i-1]
    
    vm5[i] <- vm5[i-1] - (vm5[i-1]*vm6[i-1])*kA - 
      (vC[i-1]*vm5[i-1]*vm6[i-1])*kcA + sm5 #
    
    vm6[i] <- vm5[i-1] - (vm5[i-1]*vm6[i-1])*kA - 
      (vC[i-1]*vm5[i-1]*vm6[i-1])*kcA + sm6
    
    vAi[i] <- vAi[i-1] + vA[i-1]*kAi - vAi[i-1]*kAa - dAi*vAi[i-1]
  }
  
  return(cbind(vm1, vm2, vm3, vm4, vm5, vm6, vA, vAi, vB, vC))

}

### exploring the model
ic_1 <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0)

rr_1 <- c(0.05, 0.05, 0.05, 0.5, 0.5, 0.5, 0.5, 0.05)

dr_1 <- c(0.2, 0.2, 0.2, 0.01)

sr_1 <- c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01)

ex_1 <- ACSDoMo(1000, ic_1, rr_1, dr_1, sr_1)

plot(ex_1[,7], ylim = c(0, 0.4), col = "black")
points(ex_1[,8], col = "red")
points(ex_1[,9], col = "blue")
points(ex_1[,10], col = "orange")

plot(ex_1[,7], ylim = c(0, 0.6), col = "black")
points(ex_1[,8], col = "red")
points(ex_1[,9], col = "blue")
points(ex_1[,10], col = "orange")
points(ex_1[,1], col = "gray")
points(ex_1[,3], col = "green")
points(ex_1[,5], col = "purple")

#taking supply rates to 0 
ic_2 <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0)

rr_2 <- c(0.05, 0.05, 0.05, 0.5, 0.5, 0.5, 0.5, 0.05)

dr_2 <- c(0.2, 0.2, 0.2, 0.01)

sr_2 <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

ex_2 <- ACSDoMo(1000, ic_2, rr_2, dr_2, sr_2)

plot(ex_2[,7], ylim = c(0, 0.6), col = "black")
points(ex_2[,8], col = "red")
points(ex_2[,9], col = "blue")
points(ex_2[,10], col = "orange")
points(ex_2[,1], col = "gray")
points(ex_2[,3], col = "green")
points(ex_2[,5], col = "purple")

## Model with no dormancy in A + No source, no decay
ic_nd <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0)
# m1, m2, m3, m4, m5, m6, A, Ai, B, C

rr_nd <- c(0.05, 0.05, 0.05, 0.5, 0.5, 0.5, 0, 1)
# kB, kC, kA, kcA, kcC, kcB, kAi kAa

dr_nd <- c(0, 0, 0, 0.0)
# dB, dC, dA, dAi

sr_nd <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
# sm1, sm2, sm3, sm4, sm5, sm6

ex_nd <- ACSDoMo(1000, ic_nd, rr_nd, dr_nd, sr_nd)

plot(ex_nd[,7], ylim = c(0, 0.6), col = "black")
#points(ex_nd[,8], col = "red") # no Ai
points(ex_nd[,9]+0.02, col = "blue")
points(ex_nd[,10]+0.04, col = "orange")
points(ex_nd[,1], col = "gray") # m1
points(ex_nd[,3]+0.02, col = "green") # m3
points(ex_nd[,5]+0.04, col = "purple") # m5 #ok - all good

## Model dormancy in A + No source, no decay
ic_d <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0)
# m1, m2, m3, m4, m5, m6, A, Ai, B, C

rr_d <- c(0.05, 0.05, 0.05, 0.5, 0.5, 0.5, 0.5, 0.05)
# kB, kC, kA, kcA, kcC, kcB, kAi kAa

dr_d <- c(0, 0, 0, 0.0)
# dB, dC, dA, dAi

sr_d <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
# sm1, sm2, sm3, sm4, sm5, sm6

ex_d <- ACSDoMo(1000, ic_d, rr_d, dr_d, sr_d)

plot(ex_d[,7], ylim = c(0, 0.6), col = "black")
points(ex_d[,8], col = "red") # no Ai
points(ex_d[,9], col = "blue")
points(ex_d[,10], col = "orange")
points(ex_d[,1], col = "gray") # m1
points(ex_d[,3], col = "green") # m3
points(ex_d[,5], col = "purple") # m5 #ok - all good

plot(ex_nd[,7]-(ex_d[,7]+(ex_d[,8])), col ="black", ylim = c(0, 0.11))
points(ex_nd[,9]-ex_d[,9], col = "blue") #
points(ex_nd[,10]-ex_d[,10], col ="orange")

#### Create Nice Image #### Create Nice Image #### Create Nice Image

## Fig 2a
# png(file = "C:/Current/Dine/Research/AstroDorm/AstroDorm2_Figs/Fig2a.png",
#     width = 1200, height = 1200, units = "px", res = 96*2) 
# 
# plot(ex_nd[,7], xlim = c(0,1075), ylim = c(0, 0.6), col = "black", 
#      xlab = "Model Time Step", 
#      ylab = "Non-Dimensional Concentration",
#      main = "An autocatalytic set of three molecules")
# #points(ex_nd[,8], col = "red") # no Ai
# points(ex_nd[,9]+0.02, col = "blue")
# points(ex_nd[,10]+0.04, col = "orange")
# mtext(expression(italic(C)), 
#       side = 3, line = -3.6, at = 1050, cex = 1, col = "orange")
# mtext(expression(italic(B)), 
#       side = 3, line = -4.2, at = 1050, cex = 1, col = "blue")
# mtext(expression(italic(A[italic(a)])), 
#       side = 3, line = -5, at = 1050, cex = 1)
# 
# dev.off() # turn off graphics

#Fig 2b
# png(file = "C:/Current/Dine/Research/AstroDorm/AstroDorm2_Figs/Fig2b.png",
#     width = 1200, height = 1200, units = "px", res = 96*2) 
# 
# plot(ex_d[,7], xlim =c(0,1075), ylim = c(0, 0.6), col = "black",
#      xlab = "Model Time Step", 
#      ylab = "Non-Dimensional Concentration",
#      main = "An autocatalytic set of three molecules with dormancy")
# points(ex_d[,8], col = "red") #
# points(ex_d[,9], col = "blue")
# points(ex_d[,10], col = "orange")
# mtext(expression(italic(A)[italic(i)]), 
#       side = 3, line = -6.5, at = 1050, cex = 1, col = "Red")
# mtext(expression(italic(C)), 
#       side = 3, line = -4.5, at = 1050, cex = 1, col = "orange")
# mtext(expression(italic(B)), 
#       side = 3, line = -5.3, at = 1050, cex = 1, col = "blue")
# mtext(expression(italic(A[italic(a)])), 
#       side = 3, line = -20.4, at = 1050, cex = 1)
# 
# dev.off()


## Fig 2C
# png(file = "C:/Current/Dine/Research/AstroDorm/AstroDorm2_Figs/Fig2c.png",
#     width = 1200, height = 1200, units = "px", res = 96*2) 
# 
# plot(ex_1[,7],xlim = c(0,1075),  ylim = c(0, 0.6), col = "black", 
#      xlab = "Model Time Step", ylab = "Non-Dimensional Concentration",
#      main = "Dormancy with decay and a food source")
# points(ex_1[,8], col = "red")
# points(ex_1[,9], col = "blue")
# points(ex_1[,10], col = "orange")
# mtext(expression(italic(A)[italic(i)]), 
#       side = 3, line = -11.9, at = 1050, cex = 1, col = "Red")
# mtext(expression(italic(C)), 
#       side = 3, line = -19.2, at = 1050, cex = 1, col = "orange")
# mtext(expression(italic(B)), 
#       side = 3, line = -20.0, at = 1050, cex = 1, col = "blue")
# mtext(expression(italic(A[italic(a)])), 
#       side = 3, line = -20.9, at = 1050, cex = 1)
# 
# dev.off()

## Fig 2d
# png(file = "C:/Current/Dine/Research/AstroDorm/AstroDorm2_Figs/Fig2d.png",
#     width = 1200, height = 1200, units = "px", res = 96*2) 
# 
# plot(ex_2[,7], ylim = c(0, 0.15), xlim = c(0, 1075), col = "black",
#      xlab = "Model Time Step", ylab = "Non-Dimensional Concentration",
#      main = "Dormancy with decay and without food source")
# points(ex_2[,8]+0.005, col = "red") 
# points(ex_2[,9]+0.01, col = "blue")
# points(ex_2[,10]+0.015, col = "orange")
# mtext(expression(italic(A)[italic(i)]), 
#       side = 3, line = -21.2, at = 1050, cex = 1, col = "Red")
# mtext(expression(italic(C)), 
#       side = 3, line = -19.4, at = 1050, cex = 1, col = "orange")
# mtext(expression(italic(B)), 
#       side = 3, line = -20.3, at = 1050, cex = 1, col = "blue")
# mtext(expression(italic(A[italic(a)])), 
#       side = 3, line = -21.9, at = 1050, cex = 1)
# 
# dev.off()

#plot multipanel

png(file = "C:/Current/Dine/Research/AstroDorm/AstroDorm2_Figs/Fig2.png",
    width = 1400, height = 1400, units = "px", res = 96*2)
par(mfrow=c(2,2))


#2a "An autocatalytic set of three molecules"
plot(ex_nd[,7], xlim = c(0,1075), ylim = c(0, 0.6), col = "black", 
          xlab = "Model Time Step",
          ylab = "Non-Dimensional Concentration")
points(ex_nd[,9]+0.02, col = "blue")
points(ex_nd[,10]+0.04, col = "orange")
mtext(expression(italic(C)),
           side = 3, line = -1.8, at = 1060, cex = 1, col = "orange")
mtext(expression(italic(B)),
           side = 3, line = -2.6, at = 1060, cex = 1, col = "blue")
mtext(expression(italic(A[italic(a)])),
           side = 3, line = -3.5, at = 1060, cex = 1)

#2b "An autocatalytic set of three molecules with dormancy"
plot(ex_d[,7], xlim =c(0,1075), ylim = c(0, 0.6), col = "black",
     xlab = "Model Time Step",
     ylab = "Non-Dimensional Concentration")
points(ex_d[,8], col = "red") #
points(ex_d[,9], col = "blue")
points(ex_d[,10], col = "orange")
mtext(expression(italic(A)[italic(i)]),
      side = 3, line = -4.1, at = 1060, cex = 1, col = "Red")
mtext(expression(italic(C)),
      side = 3, line = -2.8, at = 100, cex = 1, col = "orange")
mtext(expression(italic(B)),
      side = 3, line = -2.7, at = 1060, cex = 1, col = "blue")
mtext(expression(italic(A[italic(a)])),
      side = 3, line = -12, at = 1060, cex = 1)

#2c, "Dormancy with decay and a food source"
plot(ex_1[,7],xlim = c(0,1075),  ylim = c(0, 0.6), col = "black",
     xlab = "Model Time Step", ylab = "Non-Dimensional Concentration")
points(ex_1[,8], col = "red")
points(ex_1[,9], col = "blue")
points(ex_1[,10], col = "orange")
mtext(expression(italic(A)[italic(i)]),
      side = 3, line = -7, at = 1050, cex = 1, col = "Red")
mtext(expression(italic(C)),
      side = 3, line = -11.4, at = 1050, cex = 1, col = "orange")
mtext(expression(italic(B)),
      side = 3, line = -12.1, at = 1050, cex = 1, col = "blue")
mtext(expression(italic(A[italic(a)])),
      side = 3, line = -12.9, at = 1050, cex = 1)

#2d "Dormancy with decay and without food source" 
# could try putting the title on two lines
plot(ex_2[,7], ylim = c(0, 0.6), xlim = c(0, 1075), col = "black",
     xlab = "Model Time Step", ylab = "Non-Dimensional Concentration")
points(ex_2[,8]+0.005, col = "red")
points(ex_2[,9]+0.01, col = "blue")
points(ex_2[,10]+0.015, col = "orange")
mtext(expression(italic(A)[italic(i)]),
      side = 3, line = -10.0, at = 50, cex = 1, col = "Red")
mtext(expression(italic(C)),
      side = 3, line = -12.0, at = 1000, cex = 1, col = "orange")
mtext(expression(italic(B)),
      side = 3, line = -12.0, at = 1060, cex = 1, col = "blue")
mtext(expression(italic(A[italic(a)])),
      side = 3, line = -13, at = 1060, cex = 1)

dev.off()

