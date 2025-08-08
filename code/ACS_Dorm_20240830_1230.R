##An individual based model to examine dormancy, activity, death, birth

rm(list=ls())
getwd()
setwd('C:/Current/Dine/Research/AstroDorm')

# Code to examine dormancy in an autocatalytic set

# I need some reactants m1, m2, m3, m4, Pa

# m1 + m2 -> B
# m3 + m4 -> C
# C + Pa -> A + C
# A + m1 + m2 -> B + A
# B + m3 + m4 -> C + B
# A <-> Ai

#trial to make B
vm1 <- rep(NA, 1000)
vm2 <- rep(NA, 1000)
vB <- rep(NA, 1000)

vm1[1] <- 0.5  # some concentration
vm2[1] <- 0.5  # some concentration
vB[1] <- 0 # some concentration

k <- 0.1 #some reaction rate for m1 + m2 -> B

for (i in 2:length(vB)){   
  vB[i] <-  vB[i-1] + (vm1[i-1]*vm2[i-1])*k
  vm1[i] <- vm1[i-1] - (vm1[i-1]*vm2[i-1])*k #
  vm2[i] <- vm2[i-1] - (vm1[i-1]*vm2[i-1])*k #
}

#plot(vm1)
plot(vm1, ylim = c(0,1.1), col = "blue")
points(vm2, col= "red")
points(vB, col = "gray")

for (i in 2:length(vB)){   
  vB[i] <-  vB[i-1] + (vm1[i-1]*vm2[i-1])*k
  vm1[i] <- 0.0005 + vm1[i-1] - (vm1[i-1]*vm2[i-1])*k #add in supply for fun
  vm2[i] <- 0.0005 + vm2[i-1] - (vm1[i-1]*vm2[i-1])*k #add in supply for fun
}

#plot(vm1)
plot(vm1)
plot(vm2)
plot(vB)

### Need to build that out further


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
