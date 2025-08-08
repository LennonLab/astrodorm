##An analytic model to assess dormancy in autocatalytic sets

rm(list=ls())
getwd()
setwd('C:/Current/Dine/Research/AstroDorm')

# Code to examine dormancy in an autocatalytic set

# I need some reactants m1, m2, m3, m4, Pa

# m1 + m2 -> B
# m3 + m4 -> C
# m5 + m6 -> A
# C + m5 + m6 -> A + C
# A + m1 + m2 -> B + A
# B + m3 + m4 -> C + B
# A <-> Ai

#initialize vectors
vm1 <- rep(NA, 10000)
vm2 <- rep(NA, 10000)
vm3 <- rep(NA, 10000)
vm4 <- rep(NA, 10000)
vB <- rep(NA, 10000)
vm5 <- rep(NA, 10000)
vm6 <- rep(NA, 10000)
vA <- rep(NA, 10000)
vC <- rep(NA, 10000)
vAi <- rep(NA, 10000)

#starting concentrations
vm1[1] <- 0.5  # some concentration
vm2[1] <- 0.5  # some concentration
vm3[1] <- 0.5  # some concentration
vm4[1] <- 0.5  # some concentration
vB[1] <- 0  # some concentration
vm5[1] <- 0.5  # some concentration
vm6[1] <- 0.5  # some concentration
vA[1] <- 0  # some concentration
vC[1] <- 0  # some concentration
vAi[1] <- 0  # some concentration

#reaction rates
kB <- 0.05 #some reaction rate for m1 + m2 -> B
kC <- 0.05 #some reaction rate for m3 + m4 -> C
kA <- 0.05 #some reaction rate for m5 + m6 -> A
kcA <- 0.5 #some reaction rate for C + pa + pa -> A + C
kcC <- 0.5 #some reaction rate for B + m3 + m4 -> C + B
kcB <- 0.5 #some reaction rate for A + m1 + m2 -> B + A 
kAi <- 0.1 #dormancy transition for A
kAa <- 0.05 # activity transition for Ai

#can add dormancy for other complex molecules

#decay rates
dB <- 0.01 #some decay rate for B
dC <- 0.01 #some decay rate for C
dA <- 0.01 #some decay rate for A
dAi <- 0.0001 #some decay rate for Ai

#source rates
sm1 <- 0.01 #some source rate for m1
sm2 <- 0.01 #some source rate for m2
sm3 <- 0.01 #some source rate for m3
sm4 <- 0.01 #some source rate for m4
sm5 <- 0.01 #some source rate for m5
sm6 <- 0.01 #some source rate for m6

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
  
} #

plot(vm1, ylim = c(0,1.1), col = "blue")
plot(vm2, col= "red")
plot(vB, col = "gray")
plot(vA)
plot(vm5)
plot(vC)
plot(vm3)
plot(vm4)
plot(vAi)

# reaction approaches equilibrium faster in presence of catalyst

# Add in dormancy

# for (i in 2:length(vB)){   
#   vB[i] <-  vB[i-1] + (vm1[i-1]*vm2[i-1])*k
#   vm1[i] <- 0.0005 + vm1[i-1] - (vm1[i-1]*vm2[i-1])*k #add in supply for fun
#   vm2[i] <- 0.0005 + vm2[i-1] - (vm1[i-1]*vm2[i-1])*k #add in supply for fun
# }
# 
# #plot(vm1)
# plot(vm1)
# plot(vm2)
# plot(vB)

# Code for InDiMo is below for reference

# 
# InDiMo <- function(ts, ip, cm, ad, da, de, ds, br, ae) {
#   # ts = time steps #can continue add parameters
#   # ip = initial population
#   # cm = catalytic molecules, cm cannot be 0
#   # ad = stochastic active to dormant transition (%)
#   # da = stochastic dormant to active transition (%)
#   # de = death rate (%)
#   # ds = dormancy stabilizer (divides death rate by some factor)
#   # br = birth rate (%)
#   # ae = autocatalytic enhancement (multiplies birth rate)
# 
#   #Create a grid
#   H <- matrix(0, 20, 20)
# 
#   #create molecules (50) 
#   H[sample(1:400, ip, replace = FALSE)] <- 1
#   length(which(H[,] == 1)) #100
# 
#   #create Catalytic molecules
#   B <- cbind(sample(1:20, cm, replace = TRUE), #
#              sample(1:20, cm, replace = TRUE)) #catalytic molecules
# 
#   # Random movement in B - Create random movement across 100 time steps
# 
#   #create array of time steps
#   Bt <- array(dim = c(ts,cm,2)) #1000 ts, 200 catlaytic molecules, 2 dimensions
# 
#   #initialize array
#   Bt[1,,] <- B
# 
#   #move "B" across the environment
#   for (i in 2:ts){   # 
#    Bt[i,,] <- Bt[i-1,,] + cbind(sample(-1:1, cm, replace = TRUE),
#                                 sample(-1:1, cm, replace = TRUE))
#    if (min(Bt[i,,]) < 1 & max(Bt[i,,]) <= 20){
#      Bt[i,,][which(Bt[i,,] < 1 )] <- Bt[i,,][which(Bt[i,,] < 1 )] + 1
#     } else if (min(Bt[i,,]) < 1 & max(Bt[i,,]) >= 20){
#       Bt[i,,][which(Bt[i,,] < 1 )] <- Bt[i,,][which(Bt[i,,] < 1 )] + 1
#       Bt[i,,][which(Bt[i,,] > 20 )] <- Bt[i,,][which(Bt[i,,] > 20 )] - 1
#     } else if (max(Bt[i,,]) > 20) {
#       Bt[i,,][which(Bt[i,,] > 20 )] <- Bt[i,,][which(Bt[i,,] > 20 )] - 1
#     }
#     #print(c(i,max(Bt[i,,])))
#   }
# 
#   #Make movement more easily compatible with molecule positions
#   BGt <- array(dim = c(ts,20,20))
# 
#   #create array
#   for (i in 1:ts){
#    #print(i)
#    BGt[i,,] <- matrix(0, 20, 20)
#    BGt[i,,][Bt[i,,]] <- 1
# }
# 
#   #make array of autocatalytic molecules
#   G <- array(dim = c(ts,20,20))
# 
#   G[1,,] <- H
# 
#   for (i in 2:ts){
#   
#     D <- G[i,,]
#   
#     #death loop
#     for (j in 1:20){
#       for (k in 1:20){
#         dr <- sample(1:100, 1)
#         if (G[i-1,j,k] == 1 & dr <= de){
#           D[j,k] <- 0
#         } else if (G[i-1,j,k] == 2 & dr <= de/ds){ #dormancy stabilizer
#           D[j,k] <- 0 
#         } else {
#           D[j,k] <- G[i-1,j,k]
#         }
#       }
#     }
#   
#     E <- D
#   
#     #dormancy loop
#     for (j in 1:20){
#        for (k in 1:20){
#         w <- sample(1:100,1) # sample dormant to active transition
#         v <- sample(1:100,1) # sample active to dormant transition
#       
#         if (BGt[i,j,k] == 1 & D[j,k] > 0) {
#           E[j,k] = 1
#         } else if (D[j,k] == 1 & v > ad) {  #act to dormant trans
#           E[j,k] = 1
#         } else if (D[j,k] == 1 & v <= ad) { #act to dormant trans
#           E[j,k] = 2
#         } else if (D[j,k] == 2 & w > da) { #dorm to act
#           E[j,k] = 2
#         } else if (D[j,k] == 2 & w <= da) { #dorm to act
#           E[j,k] = 1 #dormancy transition loop
#         
#         } 
#       }
#     }
# 
#     #autocatalysis loop
#     l <- NULL
#       for (j in 1:20){
#         for (k in 1:20){
#         s <- sample(1:100, 1)
#         if (E[j,k] == 2) {
#           G[i,j,k] <- 2 #no autocatalysis in dormancy
#         } else if (BGt[i,j,k] == 1 & E[j,k] == 1 & s*ae <= br) { 
#           #catalytic factor
#           l <- rbind(l,c(j,k))
#           G[i,j,k] = 1 
#         } else if (E[j,k] == 1 & s <= br){
#           l <- rbind(l,c(j,k))
#           G[i,j,k] <- 1
#         } else {
#           G[i,j,k] <- E[j,k]
#         }
#       }
#     }
#   
#     if (length(l[,1]) == 100) {
#       G[i,,] <- G[i,,]
#       } else if (length(l[,1]) >= length(which(G[i,,] == 0))) {
#       G[i,,][which(G[i,,] == 0)] <- 1  
#       } else{
#       G[i,,][sample(which(G[i,,] == 0), length(l[,1]))] <- 1
#       }
#   }
#   return(G) #returns array of dormant and active molecules
#   # death has to be a little greater than 0.5*b for metastability
# }  
