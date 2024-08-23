##An individual based model to examine dormancy, activity, death, birth

rm(list=ls())
getwd()
setwd('C:/Current/Dine/Research/AstroDorm')

InDiMo <- function(ts, ip, cm, ad, da, de, ds, br, ae) {
  # ts = time steps #can continue add parameters
  # ip = initial population
  # cm = catalytic molecules, cm cannot be 0
  # ad = stochastic active to dormant transition (%)
  # da = stochastic dormant to active transition (%)
  # de = death rate (%)
  # ds = dormancy stabilizer (divides death rate by some factor)
  # br = birth rate (%)
  # ae = autocatalytic enhancement (multiplies birth rate)

  #Create a grid
  H <- matrix(0, 20, 20)

  #create molecules (50) 
  H[sample(1:400, ip, replace = FALSE)] <- 1
  length(which(H[,] == 1)) #100

  #create Catalytic molecules
  B <- cbind(sample(1:20, cm, replace = TRUE), #
             sample(1:20, cm, replace = TRUE)) #catalytic molecules

  # Random movement in B - Create random movement across 100 time steps

  #create array of time steps
  Bt <- array(dim = c(ts,cm,2)) #1000 ts, 200 catlaytic molecules, 2 dimensions

  #initialize array
  Bt[1,,] <- B

  #move "B" across the environment
  for (i in 2:ts){   # 
   Bt[i,,] <- Bt[i-1,,] + cbind(sample(-1:1, cm, replace = TRUE),
                                sample(-1:1, cm, replace = TRUE))
   if (min(Bt[i,,]) < 1 & max(Bt[i,,]) <= 20){
     Bt[i,,][which(Bt[i,,] < 1 )] <- Bt[i,,][which(Bt[i,,] < 1 )] + 1
    } else if (min(Bt[i,,]) < 1 & max(Bt[i,,]) >= 20){
      Bt[i,,][which(Bt[i,,] < 1 )] <- Bt[i,,][which(Bt[i,,] < 1 )] + 1
      Bt[i,,][which(Bt[i,,] > 20 )] <- Bt[i,,][which(Bt[i,,] > 20 )] - 1
    } else if (max(Bt[i,,]) > 20) {
      Bt[i,,][which(Bt[i,,] > 20 )] <- Bt[i,,][which(Bt[i,,] > 20 )] - 1
    }
    #print(c(i,max(Bt[i,,])))
  }

  #Make movement more easily compatible with molecule positions
  BGt <- array(dim = c(ts,20,20))

  #create array
  for (i in 1:ts){
   #print(i)
   BGt[i,,] <- matrix(0, 20, 20)
   BGt[i,,][Bt[i,,]] <- 1
}

  #make array of autocatalytic molecules
  G <- array(dim = c(ts,20,20))

  G[1,,] <- H

  for (i in 2:ts){
  
    D <- G[i,,]
  
    #death loop
    for (j in 1:20){
      for (k in 1:20){
        dr <- sample(1:100, 1)
        if (G[i-1,j,k] == 1 & dr <= de){
          D[j,k] <- 0
        } else if (G[i-1,j,k] == 2 & dr <= de/ds){ #dormancy stabilizer
          D[j,k] <- 0 
        } else {
          D[j,k] <- G[i-1,j,k]
        }
      }
    }
  
    E <- D
  
    #dormancy loop
    for (j in 1:20){
       for (k in 1:20){
        w <- sample(1:100,1) # sample dormant to active transition
        v <- sample(1:100,1) # sample active to dormant transition
      
        if (BGt[i,j,k] == 1 & D[j,k] > 0) {
          E[j,k] = 1
        } else if (D[j,k] == 1 & v > ad) {  #act to dormant trans
          E[j,k] = 1
        } else if (D[j,k] == 1 & v <= ad) { #act to dormant trans
          E[j,k] = 2
        } else if (D[j,k] == 2 & w > da) { #dorm to act
          E[j,k] = 2
        } else if (D[j,k] == 2 & w <= da) { #dorm to act
          E[j,k] = 1 #dormancy transition loop
        
        } 
      }
    }

    #autocatalysis loop
    l <- NULL
      for (j in 1:20){
        for (k in 1:20){
        s <- sample(1:100, 1)
        if (E[j,k] == 2) {
          G[i,j,k] <- 2 #no autocatalysis in dormancy
        } else if (BGt[i,j,k] == 1 & E[j,k] == 1 & s*ae <= br) { 
          #catalytic factor
          l <- rbind(l,c(j,k))
          G[i,j,k] = 1 
        } else if (E[j,k] == 1 & s <= br){
          l <- rbind(l,c(j,k))
          G[i,j,k] <- 1
        } else {
          G[i,j,k] <- E[j,k]
        }
      }
    }
  
    if (length(l[,1]) == 100) {
      G[i,,] <- G[i,,]
      } else if (length(l[,1]) >= length(which(G[i,,] == 0))) {
      G[i,,][which(G[i,,] == 0)] <- 1  
      } else{
      G[i,,][sample(which(G[i,,] == 0), length(l[,1]))] <- 1
      }
  }
  return(G) #returns array of dormant and active molecules
  # death has to be a little greater than 0.5*b for metastability
}  
# 
# #90ad, 10de
# IDM_90ad_10de <- InDiMo(1000, 200, 50, 90, 5, 10, 5, 21, 6) 
# #ts, ip, cm, ad, da, de, ds, br, ae
# 
# #model checks
# a_90ad_10de <- rep(NA, length(IDM_90ad_10de[,1,1]))
# 
# for (i in 1:length(IDM_90ad_10de[,1,1])) {
#   a_90ad_10de[i] <- length(which(IDM_90ad_10de[i,,] == 1))
# }
# 
# d_90ad_10de <- rep(NA, length(IDM_90ad_10de[,1,1]))
# 
# for (i in 1:length(IDM_90ad_10de[,1,1])) {
#   d_90ad_10de[i] <- length(which(IDM_90ad_10de[i,,] == 2))
# }
# 
# plot(a_90ad_10de) # 
# plot(d_90ad_10de) # 
# plot(a_90ad_10de/(a_90ad_10de+d_90ad_10de)) #prop active
# plot(a_90ad_10de+d_90ad_10de) # total pop 
# 
# #75ad, 10de
# IDM_75ad_10de <- InDiMo(1000, 200, 50, 75, 5, 10, 5, 21, 6) 
# #ts, ip, cm, ad, da, de, ds, br, ae
# 
# #model checks
# a_75ad_10de <- rep(NA, length(IDM_75ad_10de[,1,1]))
# 
# for (i in 1:length(IDM_75ad_10de[,1,1])) {
#   a_75ad_10de[i] <- length(which(IDM_75ad_10de[i,,] == 1))
# }
# 
# d_75ad_10de <- rep(NA, length(IDM_75ad_10de[,1,1]))
# 
# for (i in 1:length(IDM_75ad_10de[,1,1])) {
#   d_75ad_10de[i] <- length(which(IDM_75ad_10de[i,,] == 2))
# }
# 
# plot(a_75ad_10de) # 
# plot(d_75ad_10de) # 
# plot(a_75ad_10de/(a_75ad_10de+d_75ad_10de)) #prop active
# plot(a_75ad_10de+d_75ad_10de) # total pop 
# 
# #50ad, 10de
# IDM_50ad_10de <- InDiMo(1000, 200, 50, 50, 5, 10, 5, 21, 6) 
# #ts, ip, cm, ad, da, de, ds, br, ae
# 
# #model checks
# a_50ad_10de <- rep(NA, length(IDM_50ad_10de[,1,1]))
# 
# for (i in 1:length(IDM_50ad_10de[,1,1])) {
#   a_50ad_10de[i] <- length(which(IDM_50ad_10de[i,,] == 1))
# }
# 
# d_50ad_10de <- rep(NA, length(IDM_50ad_10de[,1,1]))
# 
# for (i in 1:length(IDM_50ad_10de[,1,1])) {
#   d_50ad_10de[i] <- length(which(IDM_50ad_10de[i,,] == 2))
# }
# 
# plot(a_50ad_10de) # 
# plot(d_50ad_10de) # 
# plot(a_50ad_10de/(a_50ad_10de+d_50ad_10de)) #prop active
# plot(a_50ad_10de+d_50ad_10de) # total pop 
# 
# #90ad, 15de
# IDM_90ad_15de <- InDiMo(1000, 200, 50, 90, 5, 15, 5, 21, 6) 
# #ts, ip, cm, ad, da, de, ds, br, ae
# 
# #model checks
# a_90ad_15de <- rep(NA, length(IDM_90ad_15de[,1,1]))
# 
# for (i in 1:length(IDM_90ad_15de[,1,1])) {
#   a_90ad_15de[i] <- length(which(IDM_90ad_15de[i,,] == 1))
# }
# 
# d_90ad_15de <- rep(NA, length(IDM_90ad_15de[,1,1]))
# 
# for (i in 1:length(IDM_90ad_15de[,1,1])) {
#   d_90ad_15de[i] <- length(which(IDM_90ad_15de[i,,] == 2))
# }
# 
# plot(a_90ad_15de) # 
# plot(d_90ad_15de) # 
# plot(a_90ad_15de/(a_90ad_15de+d_90ad_15de)) #prop active
# plot(a_90ad_15de+d_90ad_15de) # total pop 
# 
# #75ad, 15de
# IDM_75ad_15de <- InDiMo(1000, 200, 50, 75, 5, 15, 5, 21, 6) 
# #ts, ip, cm, ad, da, de, ds, br, ae
# 
# #model checks
# a_75ad_15de <- rep(NA, length(IDM_75ad_15de[,1,1]))
# 
# for (i in 1:length(IDM_75ad_15de[,1,1])) {
#   a_75ad_15de[i] <- length(which(IDM_75ad_15de[i,,] == 1))
# }
# 
# d_75ad_15de <- rep(NA, length(IDM_75ad_15de[,1,1]))
# 
# for (i in 1:length(IDM_75ad_15de[,1,1])) {
#   d_75ad_15de[i] <- length(which(IDM_75ad_15de[i,,] == 2))
# }
# 
# plot(a_75ad_15de) # 
# plot(d_75ad_15de) # 
# plot(a_75ad_15de/(a_75ad_15de+d_75ad_15de)) #prop active
# plot(a_75ad_15de+d_75ad_15de) # total pop 
# 
# #50ad, 15de
# IDM_50ad_15de <- InDiMo(1000, 200, 50, 50, 5, 15, 5, 21, 6) 
# #ts, ip, cm, ad, da, de, ds, br, ae
# 
# #model checks
# a_50ad_15de <- rep(NA, length(IDM_50ad_15de[,1,1]))
# 
# for (i in 1:length(IDM_50ad_15de[,1,1])) {
#   a_50ad_15de[i] <- length(which(IDM_50ad_15de[i,,] == 1))
# }
# 
# d_50ad_15de <- rep(NA, length(IDM_50ad_15de[,1,1]))
# 
# for (i in 1:length(IDM_50ad_15de[,1,1])) {
#   d_50ad_15de[i] <- length(which(IDM_50ad_15de[i,,] == 2))
# }
# 
# plot(a_50ad_15de) # 
# plot(d_50ad_15de) # 
# plot(a_50ad_15de/(a_50ad_15de+d_50ad_15de)) #prop active
# plot(a_50ad_15de+d_50ad_15de, col = "black") # total pop 
# #points(a_75ad_15de+d_75ad_15de, col = "green")
# #points(a_90ad_15de+d_90ad_15de, col = "blue")
# 
# #25ad, 15de
# IDM_25ad_15de <- InDiMo(1000, 200, 50, 25, 5, 15, 5, 21, 6) 
# #ts, ip, cm, ad, da, de, ds, br, ae
# 
# #model checks
# a_25ad_15de <- rep(NA, length(IDM_25ad_15de[,1,1]))
# 
# for (i in 1:length(IDM_25ad_15de[,1,1])) {
#   a_25ad_15de[i] <- length(which(IDM_25ad_15de[i,,] == 1))
# }
# 
# d_25ad_15de <- rep(NA, length(IDM_25ad_15de[,1,1]))
# 
# for (i in 1:length(IDM_25ad_15de[,1,1])) {
#   d_25ad_15de[i] <- length(which(IDM_25ad_15de[i,,] == 2))
# }
# 
# plot(a_25ad_15de) # 
# plot(d_25ad_15de) # 
# plot(a_25ad_15de/(a_25ad_15de+d_25ad_15de)) #prop active
# plot(a_25ad_15de+d_25ad_15de) # total pop 
# points(a_90ad_15de+d_90ad_15de, col = "blue")
# points(a_50ad_15de+d_50ad_15de, col = "green")
# 
# #25ad, 15de, 15 ds
# IDM_25ad_15de_15ds <- InDiMo(1000, 200, 50, 25, 5, 15, 15, 21, 6) 
# #ts, ip, cm, ad, da, de, ds, br, ae
# 
# #model checks
# a_25ad_15de_15ds <- rep(NA, length(IDM_25ad_15de_15ds[,1,1]))
# 
# for (i in 1:length(IDM_25ad_15de_15ds[,1,1])) {
#   a_25ad_15de_15ds[i] <- length(which(IDM_25ad_15de_15ds[i,,] == 1))
# }
# 
# d_25ad_15de_15ds <- rep(NA, length(IDM_25ad_15de_15ds[,1,1]))
# 
# for (i in 1:length(IDM_25ad_15de_15ds[,1,1])) {
#   d_25ad_15de_15ds[i] <- length(which(IDM_25ad_15de_15ds[i,,] == 2))
# }
# 
# plot(a_25ad_15de_15ds) # 
# plot(d_25ad_15de_15ds) # 
# plot(a_25ad_15de_15ds/(a_25ad_15de_15ds+d_25ad_15de_15ds)) #prop active
# plot(a_25ad_15de_15ds+d_25ad_15de_15ds, ylim= c(0,400)) # total pop 
# #points(a_90ad_15de+d_90ad_15de, col = "blue")
# #points(a_25ad_15de+d_25ad_15de, col = "green") #
# 
# #90ad, 15de, 15 ds
# IDM_90ad_15de_15ds <- InDiMo(1000, 200, 50, 90, 5, 15, 15, 21, 6) 
# #ts, ip, cm, ad, da, de, ds, br, ae
# 
# #model checks
# a_90ad_15de_15ds <- rep(NA, length(IDM_90ad_15de_15ds[,1,1]))
# 
# for (i in 1:length(IDM_90ad_15de_15ds[,1,1])) {
#   a_90ad_15de_15ds[i] <- length(which(IDM_90ad_15de_15ds[i,,] == 1))
# }
# 
# d_90ad_15de_15ds <- rep(NA, length(IDM_90ad_15de_15ds[,1,1]))
# 
# for (i in 1:length(IDM_90ad_15de_15ds[,1,1])) {
#   d_90ad_15de_15ds[i] <- length(which(IDM_90ad_15de_15ds[i,,] == 2))
# }
# 
# plot(a_90ad_15de_15ds) # 
# plot(d_90ad_15de_15ds) # 
# plot(a_90ad_15de_15ds/(a_90ad_15de_15ds+d_90ad_15de_15ds)) #prop active
# plot(a_90ad_15de_15ds+d_90ad_15de_15ds, ylim= c(0,400)) # total pop 
# points(a_25ad_15de_15ds+d_25ad_15de_15ds, col = "blue")
# points(a_25ad_15de+d_25ad_15de, col = "green") #
# points(a_90ad_15de+d_90ad_15de, col = "purple")
# 
# #25ad, 15de, 6 ds
# IDM_25ad_15de_6ds <- InDiMo(1000, 200, 50, 25, 5, 15, 6, 21, 6) 
# #ts, ip, cm, ad, da, de, ds, br, ae
# 
# #model checks
# a_25ad_15de_6ds<- rep(NA, length(IDM_25ad_15de_6ds[,1,1]))
# 
# for (i in 1:length(IDM_25ad_15de_6ds[,1,1])) {
#   a_25ad_15de_6ds[i] <- length(which(IDM_25ad_15de_6ds[i,,] == 1))
# }
# 
# d_25ad_15de_6ds <- rep(NA, length(IDM_25ad_15de_6ds[,1,1]))
# 
# for (i in 1:length(IDM_25ad_15de_6ds[,1,1])) {
#   d_25ad_15de_6ds[i] <- length(which(IDM_25ad_15de_6ds[i,,] == 2))
# }
# 
# plot(a_25ad_15de_6ds) # 
# plot(d_25ad_15de_6ds) # 
# plot(a_25ad_15de_6ds/(a_25ad_15de_6ds+d_25ad_15de_6ds)) #prop active
# plot(a_25ad_15de_6ds+d_25ad_15de_6ds, ylim= c(0,400)) # total pop 
# points(a_25ad_15de_15ds+d_25ad_15de_15ds, col = "blue")
# points(a_25ad_15de+d_25ad_15de, col = "green")
# points(a_90ad_15de+d_90ad_15de, col = "purple")
# 
# #90ad, 15de, 15 ds
# IDM_90ad_15de_16ds <- InDiMo(1000, 200, 50, 90, 5, 15, 15, 21, 6) 
# #ts, ip, cm, ad, da, de, ds, br, ae
# 
# #model checks
# a_90ad_15de_16ds <- rep(NA, length(IDM_90ad_15de_16ds[,1,1]))
# 
# for (i in 1:length(IDM_90ad_15de_16ds[,1,1])) {
#   a_90ad_15de_16ds[i] <- length(which(IDM_90ad_15de_16ds[i,,] == 1))
# }
# 
# d_90ad_15de_16ds <- rep(NA, length(IDM_90ad_15de_16ds[,1,1]))
# 
# for (i in 1:length(IDM_90ad_15de_16ds[,1,1])) {
#   d_90ad_15de_16ds[i] <- length(which(IDM_90ad_15de_16ds[i,,] == 2))
# }
# 
# plot(a_90ad_15de_16ds) # 
# plot(d_90ad_15de_16ds) # 
# plot(a_90ad_15de_16ds/(a_90ad_15de_16ds+d_90ad_15de_16ds)) #prop active
# plot(a_90ad_15de_16ds+d_90ad_15de_16ds, ylim= c(0,400)) # total pop 
# 
# #50ad, 15de, 15ds
# IDM_50ad_15de_15ds <- InDiMo(1000, 200, 50, 50, 5, 15, 15, 21, 6) 
# #ts, ip, cm, ad, da, de, ds, br, ae
# 
# #model checks
# a_50ad_15de_15ds<- rep(NA, length(IDM_50ad_15de_15ds[,1,1]))
# 
# for (i in 1:length(IDM_50ad_15de_15ds[,1,1])) {
#   a_50ad_15de_15ds[i] <- length(which(IDM_50ad_15de_15ds[i,,] == 1))
# }
# 
# d_50ad_15de_15ds <- rep(NA, length(IDM_50ad_15de_15ds[,1,1]))
# 
# for (i in 1:length(IDM_50ad_15de_15ds[,1,1])) {
#   d_50ad_15de_15ds[i] <- length(which(IDM_50ad_15de_15ds[i,,] == 2))
# }
# 
# plot(a_50ad_15de_15ds) # 
# plot(d_50ad_15de_15ds) # 
# plot(a_50ad_15de_15ds/(a_50ad_15de_15ds+d_50ad_15de_15ds)) #prop active
# plot(a_50ad_15de_15ds+d_50ad_15de_15ds, ylim= c(0,400)) # total pop 
# points(a_50ad_15de+d_50ad_15de, col = "green")
# 
# #75ad, 15de, 15ds
# IDM_75ad_15de_15ds <- InDiMo(1000, 200, 50, 75, 5, 15, 15, 21, 6) 
# #ts, ip, cm, ad, da, de, ds, br, ae
# 
# #model checks
# a_75ad_15de_15ds<- rep(NA, length(IDM_75ad_15de_15ds[,1,1]))
# 
# for (i in 1:length(IDM_75ad_15de_15ds[,1,1])) {
#   a_75ad_15de_15ds[i] <- length(which(IDM_75ad_15de_15ds[i,,] == 1))
# }
# 
# d_75ad_15de_15ds <- rep(NA, length(IDM_75ad_15de_15ds[,1,1]))
# 
# for (i in 1:length(IDM_75ad_15de_15ds[,1,1])) {
#   d_75ad_15de_15ds[i] <- length(which(IDM_75ad_15de_15ds[i,,] == 2))
# }
# 
# plot(a_75ad_15de_15ds) # 
# plot(d_75ad_15de_15ds) # 
# plot(a_75ad_15de_15ds/(a_75ad_15de_15ds+d_75ad_15de_15ds)) #prop active
# plot(a_75ad_15de_15ds+d_75ad_15de_15ds, ylim= c(0,400)) # total pop 
# points(a_75ad_15de+d_75ad_15de, col = "green")
# 
# #50ad 15de 8ds
# IDM_50ad_15de_8ds <- InDiMo(1000, 200, 50, 50, 5, 15, 8, 21, 6) 
# #ts, ip, cm, ad, da, de, ds, br, ae
# 
# #model checks 
# a_50ad_15de_8ds<- rep(NA, length(IDM_50ad_15de_8ds[,1,1]))
# 
# for (i in 1:length(IDM_50ad_15de_8ds[,1,1])) {
#   a_50ad_15de_8ds[i] <- length(which(IDM_50ad_15de_8ds[i,,] == 1))
# }
# 
# d_50ad_15de_8ds <- rep(NA, length(IDM_50ad_15de_8ds[,1,1]))
# 
# for (i in 1:length(IDM_50ad_15de_8ds[,1,1])) {
#   d_50ad_15de_8ds[i] <- length(which(IDM_50ad_15de_8ds[i,,] == 2))
# }
# 
# plot(a_50ad_15de_8ds) # 
# plot(d_50ad_15de_8ds) # 
# plot(a_50ad_15de_8ds/(a_50ad_15de_8ds+d_50ad_15de_8ds)) #prop active
# plot(a_50ad_15de_8ds+d_50ad_15de_8ds, ylim= c(0,400)) # total pop 
# points(a_50ad_15de+d_50ad_15de, col = "green")

# #50ad 15de 8ds -  99 times
# IDM_7 <- array(NA, dim = c(99, 1000, 20, 20))
# 
# for (i in 1:99) {
#   IDM_7[i,,,] <- InDiMo(1000, 200, 50, 50, 5, 15, 7, 21, 6)
#   print(i)
# }

#ts, ip, cm, ad, da, de, ds, br, ae

a_7 <- array(NA, dim = c(99, 1000))
d_7 <- array(NA, dim = c(99, 1000))

for (i in 1:length(IDM_7[,1,1,1]))  {
  for (j in 1:length(IDM_7[1,,1,1])) {
    a_7[i,j] <- length(which(IDM_7[i,j,,] == 1))
    d_7[i,j] <- length(which(IDM_7[i,j,,] == 2))
  }
  print(i)
}

apd_7_m <- rep(NA, 1000)
apd_7_sd <- rep(NA, 1000)

for (i in 1:1000){
  apd_7_m[i] <- mean(a_7[,i]+d_7[,i])
  apd_7_sd[i] <- sd(a_7[,i]+d_7[,i])
}


# #50ad 15de 15ds -  99 times 
# IDM_15 <- array(NA, dim = c(99, 1000, 20, 20))
# 
# for (i in 1:99) {
#   IDM_15[i,,,] <- InDiMo(1000, 200, 50, 50, 5, 15, 15, 21, 6)
#   print(i)
# }

#ts, ip, cm, ad, da, de, ds, br, ae

a_15 <- array(NA, dim = c(99, 1000))
d_15 <- array(NA, dim = c(99, 1000))

for (i in 1:length(IDM_15[,1,1,1]))  {
  for (j in 1:length(IDM_15[1,,1,1])) {
    a_15[i,j] <- length(which(IDM_15[i,j,,] == 1))
    d_15[i,j] <- length(which(IDM_15[i,j,,] == 2))
  }
  print(i)
}

apd_15_m <- rep(NA, 1000)
apd_15_sd <- rep(NA, 1000)

for (i in 1:1000){
  apd_15_m[i] <- mean(a_15[,i]+d_15[,i])
  apd_15_sd[i] <- sd(a_15[,i]+d_15[,i])
}

# #50ad 15de 5ds -  99 times 
# IDM_5 <- array(NA, dim = c(99, 1000, 20, 20))
# 
# for (i in 1:99) {
#   IDM_5[i,,,] <- InDiMo(1000, 200, 50, 50, 5, 15, 5, 21, 6)
#   print(i)
# }

#ts, ip, cm, ad, da, de, ds, br, ae

a_5 <- array(NA, dim = c(99, 1000))
d_5 <- array(NA, dim = c(99, 1000))

for (i in 1:length(IDM_5[,1,1,1]))  {
  for (j in 1:length(IDM_5[1,,1,1])) {
    a_5[i,j] <- length(which(IDM_5[i,j,,] == 1))
    d_5[i,j] <- length(which(IDM_5[i,j,,] == 2))
  }
  print(i)
}

apd_5_m <- rep(NA, 1000)
apd_5_sd <- rep(NA, 1000)

for (i in 1:1000){
  apd_5_m[i] <- mean(a_5[,i]+d_5[,i])
  apd_5_sd[i] <- sd(a_5[,i]+d_5[,i])
}

PM <- cbind(apd_15_m, apd_15_sd, apd_7_m, apd_7_sd, apd_5_m, apd_5_sd)

write.csv(PM, "AverageShielding.csv")

#Plot for the Fig 50ad_15de_15ds + 50ad_15de_8ds + 50ad_15de_5ds
# plot(a_50ad_15de_8ds+d_50ad_15de_15ds, pch = 19,
#      ylim = c(0,400), col = rgb(0,0,0, alpha = 0.01),
#      xlab = "Model Time Step",
#      ylab = "Total Population (A + Ad)") # total pop
# 
# for (i in 1:99){
#   points(a_p[i,]+d_p[i,], pch = 19, 
#          col = rgb(0,0,0, alpha = 0.01))
# }
# 
# 
# points(a_50ad_15de_8ds+d_50ad_15de_8ds, pch = 19, 
#        col = rgb(1,0,1, alpha = 0.01)) # total pop 
# points(a_50ad_15de+d_50ad_15de, pch = 19,
#        col = rgb(0,1,0, alpha = 0.01))
# 
# legend(600, 425, title = "Decay in Dormancy",
#        legend=c("1 %", "2 %", "3 %"),  
#        pch = 19,
#        col = c("black","purple", "green"),
#        cex = 1,
#        box.lty = 0,
#        bg= rgb(1,1,1,0),
#        y.intersp= 0.5)
