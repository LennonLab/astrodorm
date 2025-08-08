##An analytic model to assess dormancy in autocatalytic sets

rm(list=ls())
getwd()
setwd('C:/Current/Dine/Research/AstroDorm')

# Code to examine dormancy in an autocatalytic set 
# this is based off a 3 member autocatalytic loop, A can engage in dormancy and
# autocatalysis, B can interact with both A and dormant A, A catalyzes C

# Food molecules: m1, m2, m3, m4, m5, m6

# Define reactions

# m1 + m2 -> B
# m3 + m4 -> C
# m5 + m6 -> A
# C + m1 + m2 -> B + C
# A + m3 + m4 -> C + A
# B + m5 + m6 -> B + A
# B + (A + Ai) + m5 + m6 -> B + (Ai + A)
# A + m5 + m6 -> A + A
# A + m5 + m6 -> A + D
# A <-> Ai

ACSBeMo <- function(ts, ic, rr, dr, sr) {
  
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
  vD <- rep(NA, ts)
  
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
  vD[1] <- ic[11] 
  
  #reaction rates
  
  kB <- rr[1] #some reaction rate for m1 + m2 -> B
  kC <- rr[2] #some reaction rate for m3 + m4 -> C
  kA <- rr[3] #some reaction rate for m5 + m6 -> A
  kcB <- rr[4] #some reaction rate for C + m1 + m2 -> B + C
  kcC <- rr[5] #some reaction rate for A + m3 + m4 -> C + A
  kcA <- rr[6] #some reaction rate for B + m5 + m6 -> B + A
  kcBA <- rr[7] #some reaction rate for B + (A + Ai) + m5 + m6 -> B + (Ai + A) + A
  kcAA <- rr[8] #some reaction rate for A + m5 + m6 -> A + A
  kcD <- rr[9] #some reaction rate for A + m5 + m6 -> A + D
  kAi <- rr[10] #dormancy transition for A
  kAa <- rr[11] #activity transition for Ai
  
  #decay rates
  dB <- dr[1] #some decay rate for B
  dC <- dr[2] #some decay rate for C
  dA <- dr[3] #some decay rate for A
  dAi <- dr[4] #some decay rate for Ai
  dD <- dr[5] #some decay rate for D
  dm1 <- dr[6] #some decay rate for m1
  dm2 <- dr[7] #some decay rate for m2
  dm3 <- dr[8] #some decay rate for m3
  dm4 <- dr[9] #some decay rate for m4
  dm5 <- dr[10] #some decay rate for m5
  dm6 <- dr[11] #some decay rate for m6
  
  #supply rates
  sm1 <- sr[1] #supply rate for m1
  sm2 <- sr[2] #supply rate for m2
  sm3 <- sr[3] #supply rate for m3
  sm4 <- sr[4] #supply rate for m4
  sm5 <- sr[5] #supply rate for m5
  sm6 <- sr[6] #supply rate for m6
  
  for (i in 2:length(vB)){   
    vB[i] <- vB[i-1] + (vm1[i-1]*vm2[i-1])*kB + 
      (vC[i-1]*vm1[i-1]*vm2[i-1])*kcB - dB*vB[i-1]
    
    vm1[i] <- vm1[i-1] - (vm1[i-1]*vm2[i-1])*kB -
      (vC[i-1]*vm1[i-1]*vm2[i-1])*kcB + sm1 - vm1[i-1]*dr[6] #
    
    vm2[i] <- vm2[i-1] - (vm1[i-1]*vm2[i-1])*kB -
      (vC[i-1]*vm1[i-1]*vm2[i-1])*kcB +sm2 - vm2[i-1]*dr[7] #
    
    vC[i] <-  vC[i-1] + (vm3[i-1]*vm4[i-1])*kC + 
      (vA[i-1]*vm3[i-1]*vm4[i-1])*kcC - dC*vC[i-1]
    
    vm3[i] <- vm3[i-1] - (vm3[i-1]*vm4[i-1])*kC - 
      (vA[i-1]*vm3[i-1]*vm4[i-1])*kcC + sm3 - vm3[i-1]*dr[8] #
    
    vm4[i] <- vm4[i-1] - (vm3[i-1]*vm4[i-1])*kC - 
      (vA[i-1]*vm3[i-1]*vm4[i-1])*kcC + sm4 - vm4[i-1]*dr[9] #
    
    vA[i] <- vA[i-1] + vAi[i-1]*kAa + (vm5[i-1]*vm6[i-1])*kA + 
      (vB[i-1]*vm5[i-1]*vm6[i-1])*kcA + (vA[i-1]*vm5[i-1]*vm6[i-1])*kcAA +
      (vB[i-1]*(vA[i-1]+vAi[i-1])*vm5[i-1]*vm6[i-1])*kcBA - vA[i-1]*kAi - 
      dA*vA[i-1]
    
    vD[i] <- vD[i-1] + (vA[i-1]*vm5[i-1]*vm6[i-1])*kcD - dD*vD[i-1]
    
    vm5[i] <- vm5[i-1] - (vm5[i-1]*vm6[i-1])*kA - 
      (vB[i-1]*vm5[i-1]*vm6[i-1])*kcA - (vA[i-1]*vm5[i-1]*vm6[i-1])*kcAA -
      (vB[i-1]*(vA[i-1]+vAi[i-1])*vm5[i-1]*vm6[i-1])*kcBA + sm5 - 
      vm5[i-1]*dr[10] #
    
    vm6[i] <- vm6[i-1] - (vm5[i-1]*vm6[i-1])*kA - 
      (vB[i-1]*vm5[i-1]*vm6[i-1])*kcA - (vA[i-1]*vm5[i-1]*vm6[i-1])*kcAA -
      (vB[i-1]*(vA[i-1]+vAi[i-1])*vm5[i-1]*vm6[i-1])*kcBA + sm6 - 
      vm6[i-1]*dr[11] #
    
    vAi[i] <- vAi[i-1] + vA[i-1]*kAi - vAi[i-1]*kAa - dAi*vAi[i-1]
    
  }
  
  return(cbind(vm1, vm2, vm3, vm4, vm5, vm6, vA, vAi, vB, vC, vD))

}

### Running  model
### set conditions

#### First model run
ic <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0)
# m1, m2, m3, m4, m5, m6, A, Ai, B, C, D

rr <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.2, 0.1, 0.4, 0.2, 0.05) 
# kB, kC, kA, kcB, kcC, kcA, kcBA, kcAA, kcAD, kAi, kAa

dr <- c(0.2, 0.2, 0.2, 0.01, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
# dB, dC, dA, dAi, dD, dm1, dm2, dm3, dm4, dm5, dm6

sr <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
# sm1, sm2, sm3, sm4, sm5, sm6

ex <- ACSBeMo(350, ic, rr, dr, sr)

plot(ex[,7], ylim = c(0, 0.4), col = "darkred") # A
points(ex[,8], col = "red") #AI  
points(ex[,9], col = "blue") # B
points(ex[,10], col = "cyan") # C
points(ex[,11], col = "orange") # D

plot(ex[,7], ylim = c(0, 0.8), col = "darkred") # A
points(ex[,8], col = "red") #AI  
points(ex[,1], col = "purple") # m1
points(ex[,3], col = "gray") # m3
points(ex[,5], col = "darkgreen") # m5

# Run in Time varying environment

# In Env 1 A the environment functions so that A is mostly active (explain more)
# In Env 2 A is mostly dormant

ic_tv <- c(0.1, 0.1, 0.1, 0.1, 0.5, 0.5, 0, 0, 0, 0, 0)
# m1, m2, m3, m4, m5, m6, A, Ai, B, C, D

rr_tv <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.2, 0.06, 0.4, 0.01, 0.99) 
# kB, kC, kA, kcB, kcC, kcA, kcBA, kcAA, kcAD, kAi, kAa

dr_tv <- c(0.2, 0.2, 0.8, 0.01, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
# dB, dC, dA, dAi, dD, dm1, dm2, dm3, dm4, dm5, dm6

sr_tv <- c(0.0, 0.0, 0.0, 0.0, 0.05, 0.05)
# sm1, sm2, sm3, sm4, sm5, sm6

ex_tv <- ACSBeMo(50, ic_tv, rr_tv, dr_tv, sr_tv)

BeMo_tv <- ex_tv

for (i in 1:6) { 
  if (i %% 2 == 0) {
    ic_m <- BeMo_tv[50*i,]
    # m1, m2, m3, m4, m5, m6, A, Ai, B, C, D

    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.2, 0.06, 0.4, 0.01, 0.99) 
    # kB, kC, kA, kcB, kcC, kcA, kcBA, kcAA, kcAD, kAi, kAa

    dr_m <- c(0.2, 0.2, 0.8, 0.01, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
    # dB, dC, dA, dAi, dD, dm1, dm2, dm3, dm4, dm5, dm6

    sr_m <- c(0.0, 0.0, 0.0, 0.0, 0.05, 0.05)
    # sm1, sm2, sm3, sm4, sm5, sm6

    ex_m <- ACSBeMo(50, ic_m, rr_m, dr_m, sr_m)

    BeMo_tv <- rbind(BeMo_tv, ex_m)

  } else {
      ic_m <- BeMo_tv[50*i,]
      # m1, m2, m3, m4, m5, m6, A, Ai, B, C, D

      rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.2, 0.06, 0.4, 0.99, 0.01) 
      # kB, kC, kA, kcB, kcC, kcA, kcBA, kcAA, kcAD, kAi, kAa

      dr_m <- c(0.2, 0.2, 0.2, 0.01, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
      # dB, dC, dA, dAi, dD, dm1, dm2, dm3, dm4, dm5, dm6

      sr_m <- c(0.05, 0.05, 0.05, 0.05, 0.0, 0.0)
      # sm1, sm2, sm3, sm4, sm5, sm6

      ex_m <- ACSBeMo(50, ic_m, rr_m, dr_m, sr_m)

      BeMo_tv <- rbind(BeMo_tv, ex_m)
  }
}

plot(BeMo_tv[,7], ylim = c(0, 0.02), col = "darkred") # A
points(BeMo_tv[,8], col = "red") #AI  
points(BeMo_tv[,9], col = "blue") # B
points(BeMo_tv[,10], col = "cyan") # C
points(BeMo_tv[,11], col = "orange") # D 
points(BeMo_tv[,1], col = "purple") # m1
points(BeMo_tv[,3], col = "gray") # m3
points(BeMo_tv[,5], col = "darkgreen") # m5

plot(BeMo_tv[,7], ylim = c(0, 0.6), col = "darkred") # A
points(BeMo_tv[,8], col = "red") # AI
points(BeMo_tv[,1], col = "purple") # m1
points(BeMo_tv[,3], col = "gray") # m3
points(BeMo_tv[,5], col = "darkgreen") # m5
points(BeMo_tv[,9], col = "blue") # B

# 

