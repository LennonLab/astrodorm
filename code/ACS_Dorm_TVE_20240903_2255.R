##An analytic model to assess dormancy in autocatalytic sets

rm(list=ls())
getwd()
setwd('C:/Current/Dine/Research/AstroDorm')

# Code to examine dormancy in an autocatalytic set 
# in a time varying environment, no direct decay of food

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

### Running time variant model

## Trial One, Environment varies by 50 time steps
ic_1 <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0)
# m1, m2, m3, m4, m5, m6, A, Ai, B, C

rr_1 <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.5, 0.05)
# kB, kC, kA, kcA, kcC, kcB, kAi kAa

dr_1 <- c(0.2, 0.2, 0.2, 0.01)
# dB, dC, dA, dAi, dD

sr_1 <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
# sm1, sm2, sm3, sm4, sm5, sm6

ex_1 <- ACSDoMo(50, ic_1, rr_1, dr_1, sr_1)


# notes on this run 
# supply of m5 and m6 is lower than other food sources

TVE_50 <- ex_1

for (i in 1:6) {
  if (i %% 2 == 0) {
    ic_m <- TVE_50[50*i,] 
    # m1, m2, m3, m4, m5, m6, A, Ai, B, C
    
    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.5, 0.05)
    # kB, kC, kA, kcA, kcC, kcB, kAi kAa
    
    dr_m <- c(0.2, 0.2, 0.2, 0.01)
    # dB, dC, dA, dAi,
    
    sr_m <- c(0.05, 0.05, 0.05, 0.05, 0.01, 0.01)
    # sm1, sm2, sm3, sm4, sm5, sm6
    
    ex_m <- ACSDoMo(50, ic_m, rr_m, dr_m, sr_m)
    
    TVE_50 <- rbind(TVE_50, ex_m)
    
  } else {
      ic_m <- TVE_50[50*i,] 
      # m1, m2, m3, m4, m5, m6, A, Ai, B, C
      
      rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.5, 0.05)
      # kB, kC, kA, kcA, kcC, kcB, kAi kAa
      
      dr_m <- c(0.4, 0.4, 0.4, 0.01)
      # dB, dC, dA, dAi,
      
      sr_m <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
      # sm1, sm2, sm3, sm4, sm5, sm6
      
      ex_m <- ACSDoMo(50, ic_m, rr_m, dr_m, sr_m)
      
      TVE_50 <- rbind(TVE_50, ex_m)
  }
}

plot(TVE_50[,7], ylim = c(0, 1.2), col = "black") # A
points(TVE_50[,8], col = "red") # Ai
points(TVE_50[,9], col = "blue") # B
points(TVE_50[,10], col = "orange") # C
points(TVE_50[,7] + TVE_50[,8] , col = "cyan") # total A

plot(TVE_50[,7], ylim = c(0, 1.5), col = "black") # A
points(TVE_50[,8], col = "red") # Ai
points(TVE_50[,9], col = "blue") # B
points(TVE_50[,10], col = "orange") # C
points(TVE_50[,1], col = "purple") # m1
points(TVE_50[,3], col = "gray") # m3
points(TVE_50[,5], col = "darkgreen") # m5

# food does not reach steady state

## Trial Two, Environment varies by 25 time steps
# initial conditions the same 
# supply of food for A is lower than others

ex_2 <- ACSDoMo(25, ic_1, rr_1, dr_1, sr_1)

TVE_25 <- ex_2

for (i in 1:6) {
  if (i %% 2 == 0) {
    ic_m <- TVE_25[25*i,] 
    # m1, m2, m3, m4, m5, m6, A, Ai, B, C
    
    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.5, 0.05)
    # kB, kC, kA, kcA, kcC, kcB, kAi kAa
    
    dr_m <- c(0.2, 0.2, 0.2, 0.01)
    # dB, dC, dA, dAi,
    
    sr_m <- c(0.05, 0.05, 0.05, 0.05, 0.01, 0.01)
    # sm1, sm2, sm3, sm4, sm5, sm6
    
    ex_m <- ACSDoMo(25, ic_m, rr_m, dr_m, sr_m)
    
    TVE_25 <- rbind(TVE_25, ex_m)
    
  } else {
    ic_m <- TVE_25[25*i,] 
    # m1, m2, m3, m4, m5, m6, A, Ai, B, C
    
    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.5, 0.05)
    # kB, kC, kA, kcA, kcC, kcB, kAi kAa
    
    dr_m <- c(0.4, 0.4, 0.4, 0.01)
    # dB, dC, dA, dAi,
    
    sr_m <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # sm1, sm2, sm3, sm4, sm5, sm6
    
    ex_m <- ACSDoMo(25, ic_m, rr_m, dr_m, sr_m)
    
    TVE_25 <- rbind(TVE_25, ex_m)
  }
}

plot(TVE_25[,7], ylim = c(0, 1.2), col = "black") # A
points(TVE_25[,8], col = "red") # Ai
points(TVE_25[,9], col = "blue") # B
points(TVE_25[,10], col = "orange") # c
points(TVE_25[,7] + TVE_25[,8], col = "cyan") # total A

plot(TVE_25[,7], ylim = c(0, 1.6), col = "black") # A
points(TVE_25[,8], col = "red") # Ai
points(TVE_25[,1], col = "purple") # m1
points(TVE_25[,3], col = "gray") # m3
points(TVE_25[,5], col = "darkgreen") # m5

# food does not reach steady state 

## Trial Three, A responsively transitions to dormancy
# initial conditions the same 
# supply of food for A is lower than others

TVE_RT <- ex_2

for (i in 1:6) {
  if (i %% 2 == 0) {
    ic_m <- TVE_RT[25*i,] 
    # m1, m2, m3, m4, m5, m6, A, Ai, B, C
    
    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.5, 0.05)
    # kB, kC, kA, kcA, kcC, kcB, kAi kAa
    
    dr_m <- c(0.2, 0.2, 0.2, 0.01)
    # dB, dC, dA, dAi,
    
    sr_m <- c(0.05, 0.05, 0.05, 0.05, 0.01, 0.01)
    # sm1, sm2, sm3, sm4, sm5, sm6
    
    ex_m <- ACSDoMo(25, ic_m, rr_m, dr_m, sr_m)
    
    TVE_RT <- rbind(TVE_RT, ex_m)
    
  } else {
    ic_m <- TVE_RT[25*i,] 
    # m1, m2, m3, m4, m5, m6, A, Ai, B, C
    
    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.9, 0.01)
    # kB, kC, kA, kcA, kcC, kcB, kAi kAa
    
    dr_m <- c(0.4, 0.4, 0.4, 0.01)
    # dB, dC, dA, dAi,
    
    sr_m <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # sm1, sm2, sm3, sm4, sm5, sm6
    
    ex_m <- ACSDoMo(25, ic_m, rr_m, dr_m, sr_m)
    
    TVE_RT <- rbind(TVE_RT, ex_m)
  }
}

plot(TVE_RT[,7], ylim = c(0, 1.0), col = "black") # A
points(TVE_RT[,8], col = "red") # Ai
points(TVE_RT[,9], col = "blue") # B
points(TVE_RT[,10], col = "orange") # C
points(TVE_RT[,8] + TVE_RT[,7], col = "cyan") # total A

plot(TVE_RT[,7], ylim = c(0, 1.4), col = "black") # A
points(TVE_RT[,8], col = "red") # Ai
points(TVE_RT[,5], col = "darkgreen") # m5 -> A
points(TVE_RT[,3], col = "gray") # m3 -> C
points(TVE_RT[,1], col = "purple") # m1 -> B

# food does not reach steady state

plot(TVE_RT[,9]-TVE_25[,9])
plot(TVE_RT[,10]-TVE_25[,10])

# difference in models - What do I make of this?

########## Trial 4 no dormancy

# change dormancy keep everything else the same
# food sources for A have lower supply functions

rr_nd <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0, 1)
# reaction rate
# kB, kC, kA, kcA, kcC, kcB, kAi kAa

dr_nd <- c(0.2, 0.2, 0.2, 0.01)
# decay rate
# dB, dC, dA, dAi,

ex_nd <- ACSDoMo(25, ic_1, rr_nd, dr_nd, sr_1)

TVE_nd <- ex_nd

for (i in 1:6) {
  if (i %% 2 == 0) {
    ic_m <- TVE_nd[25*i,] 
    # m1, m2, m3, m4, m5, m6, A, Ai, B, C
    
    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0, 1)
    # kB, kC, kA, kcA, kcC, kcB, kAi kAa
    
    dr_m <- c(0.2, 0.2, 0.2, 0.01)
    # dB, dC, dA, dAi,
    
    sr_m <- c(0.05, 0.05, 0.05, 0.05, 0.01, 0.01)
    # sm1, sm2, sm3, sm4, sm5, sm6
    
    ex_m <- ACSDoMo(25, ic_m, rr_m, dr_m, sr_m)
    
    TVE_nd <- rbind(TVE_nd, ex_m)
    
  } else {
    ic_m <- TVE_nd[25*i,] 
    # m1, m2, m3, m4, m5, m6, A, Ai, B, C
    
    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0, 1)
    # kB, kC, kA, kcA, kcC, kcB, kAi kAa
    
    dr_m <- c(0.4, 0.4, 0.4, 0.01)
    # dB, dC, dA, dAi,
    
    sr_m <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # sm1, sm2, sm3, sm4, sm5, sm6
    
    ex_m <- ACSDoMo(25, ic_m, rr_m, dr_m, sr_m)
    
    TVE_nd <- rbind(TVE_nd, ex_m)
  }
}

plot(TVE_nd[,7], ylim = c(0, 0.4), col = "black") # A
points(TVE_nd[,8], col = "red") # Ai - fantastic is 0 in this model
points(TVE_nd[,9], col = "blue") # B
points(TVE_nd[,10], col = "orange") # C

#points(TVE_nd[,8]+ TVE_nd[,7], col = "cyan") # total A
# A is lower due to lower supply rates - cool 

plot(TVE_nd[,7], ylim = c(0, 1.6), col = "black") # A
points(TVE_nd[,8], col = "red") # Ai - fantastic is 0 in this model
points(TVE_nd[,1], col = "purple") # m1
points(TVE_nd[,3], col = "gray") # m3
points(TVE_nd[,5], col = "darkgreen") # m5 

# food does not reach steady state 
# A is lower due to lower supply rates

plot(TVE_nd[,9], col = "blue")
points(TVE_25[,9], col = "green")
points(TVE_RT[,9], col = "darkblue")

plot(TVE_nd[,7] + TVE_nd[,8], col = "blue", ylim = c(0, 0.8))
points(TVE_25[,7] + TVE_25[,8], col = "green")
points(TVE_RT[,7] + TVE_RT[,8], col = "darkblue")

