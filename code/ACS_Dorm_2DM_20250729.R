##An analytic model to assess dormancy in autocatalytic sets

rm(list=ls())
getwd()
setwd('C:/Current/Dine/Research/AstroDorm')

# Code to examine dormancy in an autocatalytic set 
# in a time varying environment with degradation of food sources

# I need some reactants m1, m2, m3, m4, Pa

# m1 + m2 -> B
# m3 + m4 -> C
# m5 + m6 -> A
# C + m5 + m6 -> A + C
# A + m1 + m2 -> B + A
# B + m3 + m4 -> C + B
# A <-> Ai
# B <-> Bi

ACS2Do <- function(ts, ic, rr, dr, sr) {
  
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
  vBi <- rep(NA, ts)
  
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
  vBi[1] <- ic[10]
  vC[1] <- ic[11] 
  
  #reaction rates
  kB <- rr[1] #some reaction rate for m1 + m2 -> B
  kC <- rr[2] #some reaction rate for m3 + m4 -> C
  kA <- rr[3] #some reaction rate for m5 + m6 -> A
  kcA <- rr[4] #some reaction rate for C + pa + pa -> A + C
  kcC <- rr[5] #some reaction rate for B + m3 + m4 -> C + B
  kcB <- rr[6] #some reaction rate for A + m1 + m2 -> B + A 
  kAi <- rr[7] #dormancy transition for A
  kAa <- rr[8] #activity transition for Ai
  kBi <- rr[9] #dormancy transition for B
  kBa <- rr[10] #activity transition for Bi
  
  #decay rates
  dB <- dr[1] #some decay rate for B
  dC <- dr[2] #some decay rate for C
  dA <- dr[3] #some decay rate for A
  dAi <- dr[4] #some decay rate for Ai
  dBi <- dr[5] #some decay rate for Bi
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
    vB[i] <- vB[i-1] + vBi[i-1]*kBa  + (vm1[i-1]*vm2[i-1])*kB + 
      (vA[i-1]*vm1[i-1]*vm2[i-1])*kcB - vB[i-1]*kBi - dB*vB[i-1]
    
    vm1[i] <- vm1[i-1] - (vm1[i-1]*vm2[i-1])*kB -
      (vA[i-1]*vm1[i-1]*vm2[i-1])*kcB + sm1 - vm1[i-1]*dr[6] #
    
    vm2[i] <- vm2[i-1] - (vm1[i-1]*vm2[i-1])*kB -
      (vA[i-1]*vm1[i-1]*vm2[i-1])*kcB +sm2 - vm2[i-1]*dr[7] #
    
    vC[i] <-  vC[i-1] + (vm3[i-1]*vm4[i-1])*kC + 
      (vB[i-1]*vm3[i-1]*vm4[i-1])*kcC - dC*vC[i-1]
    
    vm3[i] <- vm3[i-1] - (vm3[i-1]*vm4[i-1])*kB - 
      (vB[i-1]*vm3[i-1]*vm4[i-1])*kcC + sm3 - vm3[i-1]*dr[8] #
    
    vm4[i] <- vm4[i-1] - (vm3[i-1]*vm4[i-1])*kB - 
      (vB[i-1]*vm3[i-1]*vm4[i-1])*kcC + sm4 - vm4[i-1]*dr[9] #
    
    vA[i] <- vA[i-1] + vAi[i-1]*kAa + (vm5[i-1]*vm6[i-1])*kA + 
      (vC[i-1]*vm5[i-1]*vm6[i-1])*kcA - vA[i-1]*kAi - dA*vA[i-1]
    
    vm5[i] <- vm5[i-1] - (vm5[i-1]*vm6[i-1])*kA - 
      (vC[i-1]*vm5[i-1]*vm6[i-1])*kcA + sm5 - vm5[i-1]*dr[10] #
    
    vm6[i] <- vm5[i-1] - (vm5[i-1]*vm6[i-1])*kA - 
      (vC[i-1]*vm5[i-1]*vm6[i-1])*kcA + sm6 - vm6[i-1]*dr[11] #
    
    vAi[i] <- vAi[i-1] + vA[i-1]*kAi - vAi[i-1]*kAa - dAi*vAi[i-1]
    
    vBi[i] <- vBi[i-1] + vB[i-1]*kBi - vBi[i-1]*kBa - dBi*vBi[i-1]
    
  }
  
  return(cbind(vm1, vm2, vm3, vm4, vm5, vm6, vA, vAi, vB, vBi, vC))

}

### Running time variant model
## Environment varies by 50 time steps
## set initial conditions

#### First model run - no dormancy in A, no dormancy in B
ic_nd <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0)
# m1, m2, m3, m4, m5, m6, A, Ai, B, Bi, C

rr_nd <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.0, 1.0, 0.0, 1.0)
# kB, kC, kA, kcA, kcC, kcB, kAi, kAa, kBi, kBa

dr_nd <- c(0.2, 0.2, 0.2, 0.01, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
# dB, dC, dA, dAi, dBi, dm1, dm2, dm3, dm4, dm5, dm6

sr_nd <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
# sm1, sm2, sm3, sm4, sm5, sm6

ex_nd <- ACS2Do(50, ic_nd, rr_nd, dr_nd, sr_nd)

TDM_nd <- ex_nd

for (i in 1:6) {
  if (i %% 2 == 0) {
    ic_m <- TDM_nd[50*i,]
    # m1, m2, m3, m4, m5, m6, A, Ai, B, Bi, C

    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.0, 1.0, 0.0, 1.0)
    # kB, kC, kA, kcA, kcC, kcB, kAi, kAa, kBi, kBa

    dr_m <- c(0.2, 0.2, 0.2, 0.01, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    # dB, dC, dA, dAi, dBi, dm1, dm2, dm3, dm4, dm5, dm6

    sr_m <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    # sm1, sm2, sm3, sm4, sm5, sm6

    ex_m <- ACS2Do(50, ic_m, rr_m, dr_m, sr_m)

    TDM_nd <- rbind(TDM_nd, ex_m)

  } else {
      ic_m <- TDM_nd[50*i,]
      # m1, m2, m3, m4, m5, m6, A, Ai, B, Bi, C

      rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.0, 1.0, 0.0, 1.0) 
      # kB, kC, kA, kcA, kcC, kcB, kAi, kAa, kBi, kBa

      dr_m <- c(0.4, 0.4, 0.4, 0.01, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
      # dB, dC, dA, dAi, dBi, dm1, dm2, dm3, dm4, dm5, dm6

      sr_m <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
      # sm1, sm2, sm3, sm4, sm5, sm6

      ex_m <- ACS2Do(50, ic_m, rr_m, dr_m, sr_m)

      TDM_nd <- rbind(TDM_nd, ex_m)
  }
}

plot(TDM_nd[,7], ylim = c(0, 0.2), col = "darkred") # A
points(TDM_nd[,8], col = "red") #AI  
points(TDM_nd[,9], col = "blue") # B
points(TDM_nd[,10], col = "cyan") # BI
points(TDM_nd[,11], col = "orange") # C

# All dynamics the same - sweet 

plot(TDM_nd[,7], ylim = c(0, 0.8), col = "darkred") # A
points(TDM_nd[,8], col = "red") # AI
points(TDM_nd[,1], col = "purple") # m1
points(TDM_nd[,3], col = "gray") # m3
points(TDM_nd[,5], col = "darkgreen") # m5

# All food the same

### Running time variant model, two dormant molecules

#### First model run - dormancy in A and dormancy in B
ic_2d <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0)
# m1, m2, m3, m4, m5, m6, A, Ai, B, Bi, C

rr_2d <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.75, 0.05, 0.75, 0.05)
# kB, kC, kA, kcA, kcC, kcB, kAi, kAa, kBi, kBa

dr_2d <- c(0.2, 0.2, 0.2, 0.01, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
# dB, dC, dA, dAi, dBi, dm1, dm2, dm3, dm4, dm5, dm6

sr_2d <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
# sm1, sm2, sm3, sm4, sm5, sm6

ex_2d <- ACS2Do(50, ic_2d, rr_2d, dr_2d, sr_2d)

TDM_2d <- ex_2d

for (i in 1:6) {
  if (i %% 2 == 0) {
    ic_m <- TDM_2d[50*i,]
    # m1, m2, m3, m4, m5, m6, A, Ai, B, Bi, C
    
    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.75, 0.05, 0.75, 0.05)
    # kB, kC, kA, kcA, kcC, kcB, kAi, kAa, kBi, kBa
    
    dr_m <- c(0.2, 0.2, 0.2, 0.01, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    # dB, dC, dA, dAi, dBi, dm1, dm2, dm3, dm4, dm5, dm6
    
    sr_m <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    # sm1, sm2, sm3, sm4, sm5, sm6
    
    ex_m <- ACS2Do(50, ic_m, rr_m, dr_m, sr_m)
    
    TDM_2d <- rbind(TDM_2d, ex_m)
    
  } else {
    ic_m <- TDM_2d[50*i,]
    # m1, m2, m3, m4, m5, m6, A, Ai, B, Bi, C
    
    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.75, 0.05, 0.75, 0.05) 
    # kB, kC, kA, kcA, kcC, kcB, kAi, kAa, kBi, kBa
    
    dr_m <- c(0.4, 0.4, 0.4, 0.01, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    # dB, dC, dA, dAi, dBi, dm1, dm2, dm3, dm4, dm5, dm6
    
    sr_m <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # sm1, sm2, sm3, sm4, sm5, sm6
    
    ex_m <- ACS2Do(50, ic_m, rr_m, dr_m, sr_m)
    
    TDM_2d <- rbind(TDM_2d, ex_m)
  }
}

plot(TDM_2d[,7], ylim = c(0, 0.4), col = "darkred") # A
points(TDM_2d[,8], col = "red") #AI  
points(TDM_2d[,9], col = "blue") # B
points(TDM_2d[,10], col = "cyan") # BI
points(TDM_2d[,11], col = "orange") # C

# dormancy in A suppresses B - which suppresses C - feeds back into A

plot(TDM_2d[,7], ylim = c(0, 0.8), col = "darkred") # A
points(TDM_2d[,8], col = "red") # AI
points(TDM_2d[,1], col = "purple") # m1
points(TDM_2d[,3], col = "gray") # m3
points(TDM_2d[,5], col = "darkgreen") # m5


