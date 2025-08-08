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
  kcA <- rr[4] #some reaction rate for C + pa + pa -> A + C
  kcC <- rr[5] #some reaction rate for B + m3 + m4 -> C + B
  kcB <- rr[6] #some reaction rate for A + m1 + m2 -> B + A 
  kAi <- rr[7] #dormancy transition for A
  kAa <- rr[8] #activity transition for Ai
  
  #decay rates
  dB <- dr[1] #some decay rate for B
  dC <- dr[2] #some decay rate for C
  dA <- dr[3] #some decay rate for A
  dAi <- dr[4] #some decay rate for Ai
  dm1 <- dr[5] #some decay rate for m1
  dm2 <- dr[6] #some decay rate for m2
  dm3 <- dr[7] #some decay rate for m3
  dm4 <- dr[8] #some decay rate for m4
  dm5 <- dr[9] #some decay rate for m5
  dm6 <- dr[10] #some decay rate for m6
  
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
      (vA[i-1]*vm1[i-1]*vm2[i-1])*kcB + sm1 - vm1[i-1]*dr[5] #
    
    vm2[i] <- vm2[i-1] - (vm1[i-1]*vm2[i-1])*kB -
      (vA[i-1]*vm1[i-1]*vm2[i-1])*kcB +sm2 - vm2[i-1]*dr[6] #
    
    vC[i] <-  vC[i-1] + (vm3[i-1]*vm4[i-1])*kC + 
      (vB[i-1]*vm3[i-1]*vm4[i-1])*kcC - dC*vC[i-1]
    
    vm3[i] <- vm3[i-1] - (vm3[i-1]*vm4[i-1])*kB - 
      (vB[i-1]*vm3[i-1]*vm4[i-1])*kcC + sm3 - vm3[i-1]*dr[7] #
    
    vm4[i] <- vm4[i-1] - (vm3[i-1]*vm4[i-1])*kB - 
      (vB[i-1]*vm3[i-1]*vm4[i-1])*kcC + sm4 - vm4[i-1]*dr[8] #
    
    vA[i] <- vA[i-1] + vAi[i-1]*kAa + (vm5[i-1]*vm6[i-1])*kA + 
      (vC[i-1]*vm5[i-1]*vm6[i-1])*kcA - vA[i-1]*kAi - dA*vA[i-1]
    
    vm5[i] <- vm5[i-1] - (vm5[i-1]*vm6[i-1])*kA - 
      (vC[i-1]*vm5[i-1]*vm6[i-1])*kcA + sm5 - vm5[i-1]*dr[9] #
    
    vm6[i] <- vm5[i-1] - (vm5[i-1]*vm6[i-1])*kA - 
      (vC[i-1]*vm5[i-1]*vm6[i-1])*kcA + sm6 - vm6[i-1]*dr[10] #
    
    vAi[i] <- vAi[i-1] + vA[i-1]*kAi - vAi[i-1]*kAa - dAi*vAi[i-1]
  }
  
  return(cbind(vm1, vm2, vm3, vm4, vm5, vm6, vA, vAi, vB, vC))

}

### Running time variant model
## Environment varies by 50 time steps
## set initial conditions

#### First model run - no dormancy in A
ic_nd <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0)
# m1, m2, m3, m4, m5, m6, A, Ai, B, C

rr_nd <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.0, 1.0)
# kB, kC, kA, kcA, kcC, kcB, kAi kAa

dr_nd <- c(0.2, 0.2, 0.2, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
# dB, dC, dA, dAi, dm1, dm2, dm3, dm4, dm5, dm6

sr_nd <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
# sm1, sm2, sm3, sm4, sm5, sm6

ex_nd <- ACSDoMo(50, ic_nd, rr_nd, dr_nd, sr_nd)

TVE_nd <- ex_nd

for (i in 1:6) {
  if (i %% 2 == 0) {
    ic_m <- TVE_nd[50*i,]
    # m1, m2, m3, m4, m5, m6, A, Ai, B, C

    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.0, 1.0)
    # kB, kC, kA, kcA, kcC, kcB, kAi kAa

    dr_m <- c(0.2, 0.2, 0.2, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    # dB, dC, dA, dAi, dm1, dm2, dm3, dm4, dm5, dm6

    sr_m <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    # sm1, sm2, sm3, sm4, sm5, sm6

    ex_m <- ACSDoMo(50, ic_m, rr_m, dr_m, sr_m)

    TVE_nd <- rbind(TVE_nd, ex_m)

  } else {
      ic_m <- TVE_nd[50*i,]
      # m1, m2, m3, m4, m5, m6, A, Ai, B, C

      rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.0, 1.0)
      # kB, kC, kA, kcA, kcC, kcB, kAi kAa

      dr_m <- c(0.4, 0.4, 0.4, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
      # dB, dC, dA, dAi, dm1, dm2, dm3, dm4, dm5, dm6

      sr_m <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
      # sm1, sm2, sm3, sm4, sm5, sm6

      ex_m <- ACSDoMo(50, ic_m, rr_m, dr_m, sr_m)

      TVE_nd <- rbind(TVE_nd, ex_m)
  }
}

plot(TVE_nd[,7], ylim = c(0, 0.2), col = "black") # A
#points(TVE_nd[,8], col = "red") #AI - no dormant A in this model run - checks 
points(TVE_nd[,9], col = "blue") # B
points(TVE_nd[,10], col = "orange") # C
points(TVE_nd[,8]+ TVE_nd[,7], col = "cyan", ylim =c(0, 0.2)) #Total A

# All dynamics the same - sweet 

plot(TVE_nd[,7], ylim = c(0, 0.8), col = "black") # A
points(TVE_nd[,8], col = "red") # AI
points(TVE_nd[,1], col = "purple") # m1
points(TVE_nd[,3], col = "gray") # m3
points(TVE_nd[,5], col = "darkgreen") # m5

# All food the same (beautiful)

################## Second model run - Dormancy in A (A is mostly dormant)
ic_md <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0)
# m1, m2, m3, m4, m5, m6, A, Ai, B, C

rr_md <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.75, 0.05)
# kB, kC, kA, kcA, kcC, kcB, kAi kAa

dr_md <- c(0.2, 0.2, 0.2, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
# dB, dC, dA, dAi, dm1, dm2, dm3, dm4, dm5, dm6

sr_md <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
# sm1, sm2, sm3, sm4, sm5, sm6

ex_md <- ACSDoMo(50, ic_md, rr_md, dr_md, sr_md)

TVE_md <- ex_md

for (i in 1:6) {
  if (i %% 2 == 0) {
    ic_m <- TVE_md[50*i,]
    # m1, m2, m3, m4, m5, m6, A, Ai, B, C
    
    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.75, 0.05)
    # kB, kC, kA, kcA, kcC, kcB, kAi kAa
    
    dr_m <- c(0.2, 0.2, 0.2, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    # dB, dC, dA, dAi, dm1, dm2, dm3, dm4, dm5, dm6
    
    sr_m <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    # sm1, sm2, sm3, sm4, sm5, sm6
    
    ex_m <- ACSDoMo(50, ic_m, rr_m, dr_m, sr_m)
    
    TVE_md <- rbind(TVE_md, ex_m)
    
  } else {
    ic_m <- TVE_md[50*i,]
    # m1, m2, m3, m4, m5, m6, A, Ai, B, C
    
    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.75, 0.05)
    # kB, kC, kA, kcA, kcC, kcB, kAi kAa
    
    dr_m <- c(0.2, 0.2, 0.2, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    # dB, dC, dA, dAi, dm1, dm2, dm3, dm4, dm5, dm6
    
    sr_m <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # sm1, sm2, sm3, sm4, sm5, sm6
    
    ex_m <- ACSDoMo(50, ic_m, rr_m, dr_m, sr_m)
    
    TVE_md <- rbind(TVE_md, ex_m)
  }
}

plot(TVE_md[,7], ylim = c(0, 0.6), col = "black") # A
points(TVE_md[,8], col = "red") #AI 
points(TVE_md[,9], col = "blue") # B
points(TVE_md[,10], col = "orange") # C
points(TVE_md[,8]+ TVE_md[,7], col = "cyan") #Total A

plot(TVE_md[,7], ylim = c(0, 0.8), col = "black") # A
points(TVE_md[,8], col = "red") # AI
points(TVE_md[,1], col = "purple") # m1
points(TVE_md[,3], col = "gray") # m3
points(TVE_md[,5], col = "darkgreen") # m5
points(TVE_md[,9], col = "blue") # B

################## Third model run - responsive transition to dormancy
ic_rt <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0)
# m1, m2, m3, m4, m5, m6, A, Ai, B, C

rr_rt <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.05, 0.95)
# kB, kC, kA, kcA, kcC, kcB, kAi kAa

dr_rt <- c(0.2, 0.2, 0.2, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
# dB, dC, dA, dAi, dm1, dm2, dm3, dm4, dm5, dm6

sr_rt <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
# sm1, sm2, sm3, sm4, sm5, sm6

ex_rt <- ACSDoMo(50, ic_rt, rr_rt, dr_rt, sr_rt)

TVE_rt <- ex_rt

for (i in 1:6) {
  if (i %% 2 == 0) {
    ic_m <- TVE_rt[50*i,]
    # m1, m2, m3, m4, m5, m6, A, Ai, B, C
    
    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.05, 0.95)
    # kB, kC, kA, kcA, kcC, kcB, kAi kAa
    
    dr_m <- c(0.2, 0.2, 0.2, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    # dB, dC, dA, dAi, dm1, dm2, dm3, dm4, dm5, dm6
    
    sr_m <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    # sm1, sm2, sm3, sm4, sm5, sm6
    
    ex_m <- ACSDoMo(50, ic_m, rr_m, dr_m, sr_m)
    
    TVE_rt <- rbind(TVE_rt, ex_m)
    
  } else {
    ic_m <- TVE_rt[50*i,]
    # m1, m2, m3, m4, m5, m6, A, Ai, B, C
    
    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.95, 0.05)
    # kB, kC, kA, kcA, kcC, kcB, kAi kAa
    
    dr_m <- c(0.2, 0.2, 0.2, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    # dB, dC, dA, dAi, dm1, dm2, dm3, dm4, dm5, dm6
    
    sr_m <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # sm1, sm2, sm3, sm4, sm5, sm6
    
    ex_m <- ACSDoMo(50, ic_m, rr_m, dr_m, sr_m)
    
    TVE_rt <- rbind(TVE_rt, ex_m)
  }
}

plot(TVE_rt[,7], ylim = c(0, 0.6), col = "black") # A
points(TVE_rt[,8], col = "red") #AI 

plot(TVE_rt[,9], col = "blue", ylim = c(0, 0.2)) # B
points(TVE_rt[,10], col = "orange") # C
points(TVE_rt[,8]+ TVE_rt[,7], col = "cyan") #Total A

plot(TVE_rt[,7], ylim = c(0, 0.8), col = "black") # A
points(TVE_rt[,8], col = "red") # AI
points(TVE_rt[,1], col = "purple") # m1
points(TVE_rt[,3], col = "gray") # m3
points(TVE_rt[,5], col = "darkgreen") # m5
points(TVE_rt[,9], col = "blue") # B

################## Fourth model run - 
# responsive transition to dormancy
# and no loss in dormancy
ic_nl <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0)
# m1, m2, m3, m4, m5, m6, A, Ai, B, C

rr_nl <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.05, 0.95)
# kB, kC, kA, kcA, kcC, kcB, kAi kAa

dr_nl <- c(0.2, 0.2, 0.2, 0.0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
# dB, dC, dA, dAi, dm1, dm2, dm3, dm4, dm5, dm6

sr_nl <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
# sm1, sm2, sm3, sm4, sm5, sm6

ex_nl <- ACSDoMo(50, ic_nl, rr_nl, dr_nl, sr_nl)

TVE_nl <- ex_nl

for (i in 1:6) {
  if (i %% 2 == 0) {
    ic_m <- TVE_nl[50*i,]
    # m1, m2, m3, m4, m5, m6, A, Ai, B, C
    
    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.05, 0.95)
    # kB, kC, kA, kcA, kcC, kcB, kAi kAa
    
    dr_m <- c(0.2, 0.2, 0.2, 0.0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    # dB, dC, dA, dAi, dm1, dm2, dm3, dm4, dm5, dm6
    
    sr_m <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    # sm1, sm2, sm3, sm4, sm5, sm6
    
    ex_m <- ACSDoMo(50, ic_m, rr_m, dr_m, sr_m)
    
    TVE_nl <- rbind(TVE_nl, ex_m)
    
  } else {
    ic_m <- TVE_nl[50*i,]
    # m1, m2, m3, m4, m5, m6, A, Ai, B, C
    
    rr_m <- c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 1.0, 0.0)
    # kB, kC, kA, kcA, kcC, kcB, kAi kAa
    
    dr_m <- c(0.2, 0.2, 0.2, 0.0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
    # dB, dC, dA, dAi, dm1, dm2, dm3, dm4, dm5, dm6
    
    sr_m <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # sm1, sm2, sm3, sm4, sm5, sm6
    
    ex_m <- ACSDoMo(50, ic_m, rr_m, dr_m, sr_m)
    
    TVE_nl <- rbind(TVE_nl, ex_m)
  }
}

plot(TVE_nl[,7], ylim = c(0, 0.6), col = "black") # A
points(TVE_nl[,8], col = "red") #AI 
points(TVE_nl[,9], col = "blue") # B
points(TVE_nl[,10], col = "orange") # C
points(TVE_nl[,8]+ TVE_nl[,7], col = "cyan") #Total A

plot(TVE_nl[,7], ylim = c(0, 0.8), col = "black") # A
points(TVE_nl[,8], col = "red") # AI
points(TVE_nl[,1], col = "purple") # m1
points(TVE_nl[,3], col = "gray") # m3
points(TVE_nl[,5], col = "darkgreen") # m5
points(TVE_nl[,9], col = "blue") # B
points(TVE_nl[,9]+ TVE_nl[,7] +TVE_nl[,10], col = "magenta") # T aACS

plot(TVE_nl[,7]+TVE_nl[,8], ylim = c(0, 0.6), col = "black") # A responsive nl
points(TVE_rt[,7]+TVE_rt[,8], col = "blue") # A responsive l
points(TVE_md[,7]+TVE_md[,8], col = "darkgreen") # A mostly dormant
points(TVE_nd[,7]+TVE_nd[,8], col = "orange") # A not dormant

plot(TVE_nl[,9], ylim = c(0, 0.2), col = "black") # A responsive nl
points(TVE_rt[,9], col = "blue") # A responsive l
points(TVE_md[,9], col = "darkgreen") # A mostly dormant
points(TVE_nd[,9], col = "orange") # A not dormant

plot(TVE_nl[,10], ylim = c(0, 0.2), col = "black") # A responsive nl
points(TVE_rt[,10], col = "blue") # A responsive l
points(TVE_md[,10], col = "darkgreen") # A mostly dormant
points(TVE_nd[,10], col = "orange") # A not dormant

plot(TVE_nl[,7], ylim = c(0, 0.2), col = "black") # A responsive nl
points(TVE_rt[,7], col = "blue") # A responsive l
points(TVE_md[,7], col = "darkgreen") # A mostly dormant
points(TVE_nd[,7], col = "orange") # A not dormant

par(mfrow=c(2,2)) #set up the plotting space
plot(TVE_nd[,7], ylim = c(0, 0.6), col = "darkgreen",lty = 1, type = "l",
     lwd = 3,  ylab = "Concentration", xlab = "Time Step", 
     main = "No Dormancy")
points(TVE_nd[,9], col = "blue", lty = 1, type = "l", lwd = 3) # B
points(TVE_nd[,10], col = "darkmagenta", lty = 1, type = "l", lwd = 3) # C
points(TVE_nd[,8], col = "darkgreen", lty = 3, type = "l", lwd = 3)
points(TVE_nd[,7], col = "darkgreen", lty = 1, type = "l", lwd = 3)

plot(TVE_md[,7], ylim = c(0, 0.6), col = "darkgreen",lty = 1, type = "l",
     lwd = 3,  ylab = "Concentration", xlab = "Time Step", 
     main = "Sochastic Dormancy")
points(TVE_md[,8], col = "darkgreen", lty = 3, type = "l", lwd = 3)
points(TVE_md[,9], col = "blue", lty = 1, type = "l", lwd = 3) # B
points(TVE_md[,10], col = "darkmagenta", lty = 1, type = "l", lwd = 3) # C

plot(TVE_rt[,7], ylim = c(0, 0.6), col = "darkgreen",lty = 1, type = "l",
     lwd = 3,  ylab = "Concentration", xlab = "Time Step", 
     main = "Responsive transition, some Loss in Dormancy")
points(TVE_rt[,8], col = "darkgreen", lty = 3, type = "l", lwd = 3)
points(TVE_rt[,9], col = "blue", lty = 1, type = "l", lwd = 3) # B
points(TVE_rt[,10], col = "darkmagenta", lty = 1, type = "l", lwd = 3) # C

plot(TVE_nl[,7], ylim = c(0, 0.6), col = "darkgreen",lty = 1, type = "l",
    lwd = 3,  ylab = "Concentration", xlab = "Time Step", 
    main = "Responsive transition, No Loss in Dormancy")
points(TVE_nl[,8], col = "darkgreen", lty = 3, type = "l", lwd = 3)
points(TVE_nl[,9], col = "blue", lty = 1, type = "l", lwd = 3) # B
points(TVE_nl[,10], col = "darkmagenta", lty = 1, type = "l", lwd = 3) # C

