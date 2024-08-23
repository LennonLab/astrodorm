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

#10cm, 90ad
IDM_10cm_90ad <- InDiMo(1000, 200, 10, 90, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_10cm_90ad <- rep(NA, length(IDM_10cm_90ad[,1,1]))

for (i in 1:length(IDM_10cm_90ad[,1,1])) {
  a_10cm_90ad[i] <- length(which(IDM_10cm_90ad[i,,] == 1))
}

d_10cm_90ad <- rep(NA, length(IDM_10cm_90ad[,1,1]))

for (i in 1:length(IDM_10cm_90ad[,1,1])) {
  d_10cm_90ad[i] <- length(which(IDM_10cm_90ad[i,,] == 2))
}

plot(a_10cm_90ad) # low but growing slowly
plot(d_10cm_90ad) # pop growth
plot(a_10cm_90ad/(a_10cm_90ad+d_10cm_90ad)) #prop active
plot(a_10cm_90ad+d_10cm_90ad) # total pop 

# 50cm; 90ad
IDM_50cm_90ad <- InDiMo(1000, 200, 50, 90, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_50cm_90ad <- rep(NA, length(IDM_50cm_90ad[,1,1]))

for (i in 1:length(IDM_50cm_90ad[,1,1])) {
  a_50cm_90ad[i] <- length(which(IDM_50cm_90ad[i,,] == 1))
}

d_50cm_90ad <- rep(NA, length(IDM_50cm_90ad[,1,1]))

for (i in 1:length(IDM_50cm_90ad[,1,1])) {
  d_50cm_90ad[i] <- length(which(IDM_50cm_90ad[i,,] == 2))
}

plot(a_50cm_90ad) # 
plot(d_50cm_90ad) #
plot(a_50cm_90ad/(a_50cm_90ad+d_50cm_90ad)) #prop active
plot(a_50cm_90ad+d_50cm_90ad) # total pop 

# 150cm; 90ad
IDM_150cm_90ad <- InDiMo(1000, 200, 150, 90, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_150cm_90ad <- rep(NA, length(IDM_150cm_90ad[,1,1]))

for (i in 1:length(IDM_150cm_90ad[,1,1])) {
  a_150cm_90ad[i] <- length(which(IDM_150cm_90ad[i,,] == 1))
}

d_150cm_90ad <- rep(NA, length(IDM_150cm_90ad[,1,1]))

for (i in 1:length(IDM_150cm_90ad[,1,1])) {
  d_150cm_90ad[i] <- length(which(IDM_150cm_90ad[i,,] == 2))
}

plot(a_150cm_90ad) # low but growing slowly
plot(d_150cm_90ad) # pop growth
plot(a_150cm_90ad/(a_150cm_90ad+d_150cm_90ad)) #prop active
plot(a_150cm_90ad+d_150cm_90ad) # total pop 

# 250cm; 90ad
IDM_250cm_90ad <- InDiMo(1000, 200, 250, 90, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_250cm_90ad <- rep(NA, length(IDM_250cm_90ad[,1,1]))

for (i in 1:length(IDM_250cm_90ad[,1,1])) {
  a_250cm_90ad[i] <- length(which(IDM_250cm_90ad[i,,] == 1))
}

d_250cm_90ad <- rep(NA, length(IDM_250cm_90ad[,1,1]))

for (i in 1:length(IDM_250cm_90ad[,1,1])) {
  d_250cm_90ad[i] <- length(which(IDM_250cm_90ad[i,,] == 2))
}

plot(a_250cm_90ad) # 
plot(d_250cm_90ad) #
plot(a_250cm_90ad/(a_250cm_90ad+d_250cm_90ad)) #prop active
plot(a_250cm_90ad+d_250cm_90ad) # total pop 

# 350cm; 90ad
IDM_350cm_90ad <- InDiMo(1000, 200, 350, 90, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_350cm_90ad <- rep(NA, length(IDM_350cm_90ad[,1,1]))

for (i in 1:length(IDM_350cm_90ad[,1,1])) {
  a_350cm_90ad[i] <- length(which(IDM_350cm_90ad[i,,] == 1))
}

d_350cm_90ad <- rep(NA, length(IDM_350cm_90ad[,1,1]))

for (i in 1:length(IDM_350cm_90ad[,1,1])) {
  d_350cm_90ad[i] <- length(which(IDM_350cm_90ad[i,,] == 2))
}

plot(a_350cm_90ad) # low but growing slowly
plot(d_350cm_90ad) # pop growth
plot(a_350cm_90ad/(a_350cm_90ad+d_350cm_90ad)) #prop active
plot(a_350cm_90ad+d_350cm_90ad) # total pop

# 10cm; 75ad
IDM_10cm_75ad <- InDiMo(1000, 200, 10, 75, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_10cm_75ad <- rep(NA, length(IDM_10cm_75ad[,1,1]))

for (i in 1:length(IDM_10cm_75ad[,1,1])) {
  a_10cm_75ad[i] <- length(which(IDM_10cm_75ad[i,,] == 1))
}

d_10cm_75ad <- rep(NA, length(IDM_10cm_75ad[,1,1]))

for (i in 1:length(IDM_10cm_75ad[,1,1])) {
  d_10cm_75ad[i] <- length(which(IDM_10cm_75ad[i,,] == 2))
}

plot(a_10cm_75ad) # low but growing slowly
plot(d_10cm_75ad) # pop growth
plot(a_10cm_75ad/(a_10cm_75ad+d_10cm_75ad)) #prop active
plot(a_10cm_75ad+d_10cm_75ad) # total pop

# 50cm; 75ad
IDM_50cm_75ad <- InDiMo(1000, 200, 50, 75, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_50cm_75ad <- rep(NA, length(IDM_50cm_75ad[,1,1]))

for (i in 1:length(IDM_50cm_75ad[,1,1])) {
  a_50cm_75ad[i] <- length(which(IDM_50cm_75ad[i,,] == 1))
}

d_50cm_75ad <- rep(NA, length(IDM_50cm_75ad[,1,1]))

for (i in 1:length(IDM_50cm_75ad[,1,1])) {
  d_50cm_75ad[i] <- length(which(IDM_50cm_75ad[i,,] == 2))
}

plot(a_50cm_75ad) # low but growing slowly
plot(d_50cm_75ad) # pop growth
plot(a_50cm_75ad/(a_50cm_75ad+d_50cm_75ad)) #prop active
plot(a_50cm_75ad+d_50cm_75ad) # total pop

# 150cm; 75ad
IDM_150cm_75ad <- InDiMo(1000, 200, 150, 75, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_150cm_75ad <- rep(NA, length(IDM_150cm_75ad[,1,1]))

for (i in 1:length(IDM_150cm_75ad[,1,1])) {
  a_150cm_75ad[i] <- length(which(IDM_150cm_75ad[i,,] == 1))
}

d_150cm_75ad <- rep(NA, length(IDM_150cm_75ad[,1,1]))

for (i in 1:length(IDM_150cm_75ad[,1,1])) {
  d_150cm_75ad[i] <- length(which(IDM_150cm_75ad[i,,] == 2))
}

plot(a_150cm_75ad) # low but growing slowly
plot(d_150cm_75ad) # pop growth
plot(a_150cm_75ad/(a_150cm_75ad+d_150cm_75ad)) #prop active
plot(a_150cm_75ad+d_150cm_75ad) # total pop

# 250cm; 75ad
IDM_250cm_75ad <- InDiMo(1000, 200, 250, 75, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_250cm_75ad <- rep(NA, length(IDM_250cm_75ad[,1,1]))

for (i in 1:length(IDM_250cm_75ad[,1,1])) {
  a_250cm_75ad[i] <- length(which(IDM_250cm_75ad[i,,] == 1))
}

d_250cm_75ad <- rep(NA, length(IDM_250cm_75ad[,1,1]))

for (i in 1:length(IDM_250cm_75ad[,1,1])) {
  d_250cm_75ad[i] <- length(which(IDM_250cm_75ad[i,,] == 2))
}

plot(a_250cm_75ad) # low but growing slowly
plot(d_250cm_75ad) # pop growth
plot(a_250cm_75ad/(a_250cm_75ad+d_250cm_75ad)) #prop active
plot(a_250cm_75ad+d_250cm_75ad) # total pop

# 350cm; 75ad
IDM_350cm_75ad <- InDiMo(1000, 200, 350, 75, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_350cm_75ad <- rep(NA, length(IDM_350cm_75ad[,1,1]))

for (i in 1:length(IDM_350cm_75ad[,1,1])) {
  a_350cm_75ad[i] <- length(which(IDM_350cm_75ad[i,,] == 1))
}

d_350cm_75ad <- rep(NA, length(IDM_350cm_75ad[,1,1]))

for (i in 1:length(IDM_350cm_75ad[,1,1])) {
  d_350cm_75ad[i] <- length(which(IDM_350cm_75ad[i,,] == 2))
}

plot(a_350cm_75ad) # low but growing slowly
plot(d_350cm_75ad) # pop growth
plot(a_350cm_75ad/(a_350cm_75ad+d_350cm_75ad)) #prop active
plot(a_350cm_75ad+d_350cm_75ad) # total pop

# 10cm; 50ad
IDM_10cm_50ad <- InDiMo(1000, 200, 10, 50, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_10cm_50ad <- rep(NA, length(IDM_10cm_50ad[,1,1]))

for (i in 1:length(IDM_10cm_50ad[,1,1])) {
  a_10cm_50ad[i] <- length(which(IDM_10cm_50ad[i,,] == 1))
}

d_10cm_50ad <- rep(NA, length(IDM_10cm_50ad[,1,1]))

for (i in 1:length(IDM_10cm_50ad[,1,1])) {
  d_10cm_50ad[i] <- length(which(IDM_10cm_50ad[i,,] == 2))
}

plot(a_10cm_50ad) # low but growing slowly
plot(d_10cm_50ad) # pop growth
plot(a_10cm_50ad/(a_10cm_50ad+d_10cm_50ad)) #prop active
plot(a_10cm_50ad+d_10cm_50ad) # total pop

# 50cm; 50ad
IDM_50cm_50ad <- InDiMo(1000, 200, 50, 50, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_50cm_50ad <- rep(NA, length(IDM_50cm_50ad[,1,1]))

for (i in 1:length(IDM_50cm_50ad[,1,1])) {
  a_50cm_50ad[i] <- length(which(IDM_50cm_50ad[i,,] == 1))
}

d_50cm_50ad <- rep(NA, length(IDM_50cm_50ad[,1,1]))

for (i in 1:length(IDM_50cm_50ad[,1,1])) {
  d_50cm_50ad[i] <- length(which(IDM_50cm_50ad[i,,] == 2))
}

plot(a_50cm_50ad) # low but growing slowly
plot(d_50cm_50ad) # pop growth
plot(a_50cm_50ad/(a_50cm_50ad+d_50cm_50ad)) #prop active
plot(a_50cm_50ad+d_50cm_50ad) # total pop

# 150cm; 50ad
IDM_150cm_50ad <- InDiMo(1000, 200, 150, 50, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_150cm_50ad <- rep(NA, length(IDM_150cm_50ad[,1,1]))

for (i in 1:length(IDM_150cm_50ad[,1,1])) {
  a_150cm_50ad[i] <- length(which(IDM_150cm_50ad[i,,] == 1))
}

d_150cm_50ad <- rep(NA, length(IDM_150cm_50ad[,1,1]))

for (i in 1:length(IDM_150cm_50ad[,1,1])) {
  d_150cm_50ad[i] <- length(which(IDM_150cm_50ad[i,,] == 2))
}

plot(a_150cm_50ad) # low but growing slowly
plot(d_150cm_50ad) # pop growth
plot(a_150cm_50ad/(a_150cm_50ad+d_150cm_50ad)) #prop active
plot(a_150cm_50ad+d_150cm_50ad) # total pop

# 250cm; 50ad
IDM_250cm_50ad <- InDiMo(1000, 200, 250, 50, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_250cm_50ad <- rep(NA, length(IDM_250cm_50ad[,1,1]))

for (i in 1:length(IDM_250cm_50ad[,1,1])) {
  a_250cm_50ad[i] <- length(which(IDM_250cm_50ad[i,,] == 1))
}

d_250cm_50ad <- rep(NA, length(IDM_250cm_50ad[,1,1]))

for (i in 1:length(IDM_250cm_50ad[,1,1])) {
  d_250cm_50ad[i] <- length(which(IDM_250cm_50ad[i,,] == 2))
}

plot(a_250cm_50ad) # low but growing slowly
plot(d_250cm_50ad) # pop growth
plot(a_250cm_50ad/(a_250cm_50ad+d_250cm_50ad)) #prop active
plot(a_250cm_50ad+d_250cm_50ad) # total pop

# 300cm; 50ad
IDM_350cm_50ad <- InDiMo(1000, 200, 350, 50, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_350cm_50ad <- rep(NA, length(IDM_350cm_50ad[,1,1]))

for (i in 1:length(IDM_350cm_50ad[,1,1])) {
  a_350cm_50ad[i] <- length(which(IDM_350cm_50ad[i,,] == 1))
}

d_350cm_50ad <- rep(NA, length(IDM_350cm_50ad[,1,1]))

for (i in 1:length(IDM_350cm_50ad[,1,1])) {
  d_350cm_50ad[i] <- length(which(IDM_350cm_50ad[i,,] == 2))
}

plot(a_350cm_50ad) # low but growing slowly
plot(d_350cm_50ad) # pop growth
plot(a_350cm_50ad/(a_350cm_50ad+d_350cm_50ad)) #prop active
plot(a_350cm_50ad+d_350cm_50ad) # total pop

# 10cm; 25ad
IDM_10cm_25ad <- InDiMo(1000, 200, 10, 25, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_10cm_25ad <- rep(NA, length(IDM_10cm_25ad[,1,1]))

for (i in 1:length(IDM_10cm_25ad[,1,1])) {
  a_10cm_25ad[i] <- length(which(IDM_10cm_25ad[i,,] == 1))
}

d_10cm_25ad <- rep(NA, length(IDM_10cm_25ad[,1,1]))

for (i in 1:length(IDM_10cm_25ad[,1,1])) {
  d_10cm_25ad[i] <- length(which(IDM_10cm_25ad[i,,] == 2))
}

plot(a_10cm_25ad) # low but growing slowly
plot(d_10cm_25ad) # pop growth
plot(a_10cm_25ad/(a_10cm_25ad+d_10cm_25ad)) #prop active
plot(a_10cm_25ad+d_10cm_25ad) # total pop

# 50cm; 25ad
IDM_50cm_25ad <- InDiMo(1000, 200, 50, 25, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_50cm_25ad <- rep(NA, length(IDM_50cm_25ad[,1,1]))

for (i in 1:length(IDM_50cm_25ad[,1,1])) {
  a_50cm_25ad[i] <- length(which(IDM_50cm_25ad[i,,] == 1))
}

d_50cm_25ad <- rep(NA, length(IDM_50cm_25ad[,1,1]))

for (i in 1:length(IDM_50cm_25ad[,1,1])) {
  d_50cm_25ad[i] <- length(which(IDM_50cm_25ad[i,,] == 2))
}

plot(a_50cm_25ad) # low but growing slowly
plot(d_50cm_25ad) # pop growth
plot(a_50cm_25ad/(a_50cm_25ad+d_50cm_25ad)) #prop active
plot(a_50cm_25ad+d_50cm_25ad) # total pop

# 150cm; 25ad
IDM_150cm_25ad <- InDiMo(1000, 200, 150, 25, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_150cm_25ad <- rep(NA, length(IDM_150cm_25ad[,1,1]))

for (i in 1:length(IDM_150cm_25ad[,1,1])) {
  a_150cm_25ad[i] <- length(which(IDM_150cm_25ad[i,,] == 1))
}

d_150cm_25ad <- rep(NA, length(IDM_150cm_25ad[,1,1]))

for (i in 1:length(IDM_150cm_25ad[,1,1])) {
  d_150cm_25ad[i] <- length(which(IDM_150cm_25ad[i,,] == 2))
}

plot(a_150cm_25ad) # 
plot(d_150cm_25ad) # 
plot(a_150cm_25ad/(a_150cm_25ad+d_150cm_25ad)) #prop active
plot(a_150cm_25ad+d_150cm_25ad) # total pop

# 250cm; 25ad
IDM_250cm_25ad <- InDiMo(1000, 200, 250, 25, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_250cm_25ad <- rep(NA, length(IDM_250cm_25ad[,1,1]))

for (i in 1:length(IDM_250cm_25ad[,1,1])) {
  a_250cm_25ad[i] <- length(which(IDM_250cm_25ad[i,,] == 1))
}

d_250cm_25ad <- rep(NA, length(IDM_250cm_25ad[,1,1]))

for (i in 1:length(IDM_250cm_25ad[,1,1])) {
  d_250cm_25ad[i] <- length(which(IDM_250cm_25ad[i,,] == 2))
}

plot(a_250cm_25ad) # 
plot(d_250cm_25ad) # 
plot(a_250cm_25ad/(a_250cm_25ad+d_250cm_25ad)) #prop active
plot(a_250cm_25ad+d_250cm_25ad) # total pop

# 350cm; 25ad
IDM_350cm_25ad <- InDiMo(1000, 200, 350, 25, 5, 10, 5, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_350cm_25ad <- rep(NA, length(IDM_350cm_25ad[,1,1]))

for (i in 1:length(IDM_350cm_25ad[,1,1])) {
  a_350cm_25ad[i] <- length(which(IDM_350cm_25ad[i,,] == 1))
}

d_350cm_25ad <- rep(NA, length(IDM_350cm_25ad[,1,1]))

for (i in 1:length(IDM_350cm_25ad[,1,1])) {
  d_350cm_25ad[i] <- length(which(IDM_350cm_25ad[i,,] == 2))
}

plot(a_350cm_25ad) # 
plot(d_350cm_25ad) # 
plot(a_350cm_25ad/(a_350cm_25ad+d_350cm_25ad)) #prop active
plot(a_350cm_25ad+d_350cm_25ad) # total pop


#create semi fancy plots for kind of nice visualization
m10cm_90ad <- mean(a_10cm_90ad[3:150]/
                     (a_10cm_90ad[3:150]+d_10cm_90ad[3:150]))
m50cm_90ad <- mean(a_50cm_90ad[3:150]/
                     (a_50cm_90ad[3:150]+d_50cm_90ad[3:150]))
m150cm_90ad <- mean(a_150cm_90ad[3:150]/
                      (a_150cm_90ad[3:150]+d_150cm_90ad[3:150]))
m250cm_90ad <- mean(a_250cm_90ad[3:150]/
                      (a_250cm_90ad[3:150]+d_250cm_90ad[3:150]))
m350cm_90ad <- mean(a_350cm_90ad[3:150]/
                      (a_350cm_90ad[3:150]+d_350cm_90ad[3:150]))

m10cm_75ad <- mean(a_10cm_75ad[3:150]/
                     (a_10cm_75ad[3:150]+d_10cm_75ad[3:150]))
m50cm_75ad <- mean(a_50cm_75ad[3:150]/
                     (a_50cm_75ad[3:150]+d_50cm_75ad[3:150]))
m150cm_75ad <- mean(a_150cm_75ad[3:150]/
                      (a_150cm_75ad[3:150]+d_150cm_75ad[3:150]))
m250cm_75ad <- mean(a_250cm_75ad[3:150]/
                      (a_250cm_75ad[3:150]+d_250cm_75ad[3:150]))
m350cm_75ad <- mean(a_350cm_75ad[3:150]/
                      (a_350cm_75ad[3:150]+d_350cm_75ad[3:150]))

m10cm_50ad <- mean(a_10cm_50ad[3:150]/
                     (a_10cm_50ad[3:150]+d_10cm_50ad[3:150]))
m50cm_50ad <- mean(a_50cm_50ad[3:150]/
                     (a_50cm_50ad[3:150]+d_50cm_50ad[3:150]))
m150cm_50ad <- mean(a_150cm_50ad[3:150]/
                      (a_150cm_50ad[3:150]+d_150cm_50ad[3:150]))
m250cm_50ad <- mean(a_250cm_50ad[3:150]/
                      (a_250cm_50ad[3:150]+d_250cm_50ad[3:150]))
m350cm_50ad <- mean(a_350cm_50ad[3:150]/
                      (a_350cm_50ad[3:150]+d_350cm_50ad[3:150]))

m10cm_25ad <- mean(a_10cm_25ad[3:150]/
                     (a_10cm_25ad[3:150]+d_10cm_25ad[3:150]))
m50cm_25ad <- mean(a_50cm_25ad[3:150]/
                      (a_50cm_25ad[3:150]+d_50cm_25ad[3:150]))
m150cm_25ad <- mean(a_150cm_25ad[3:150]/
                      (a_150cm_25ad[3:150]+d_150cm_25ad[3:150]))
m250cm_25ad <- mean(a_250cm_25ad[3:150]/
                      (a_250cm_25ad[3:150]+d_250cm_25ad[3:150]))
m350cm_25ad <- mean(a_350cm_25ad[3:150]/
                      (a_350cm_25ad[3:150]+d_350cm_25ad[3:150]))

# semi fancy plot 
v1 <- cbind(m10cm_90ad, m50cm_90ad, m150cm_90ad, m250cm_90ad, m350cm_90ad)
v2 <- cbind(m10cm_75ad, m50cm_75ad, m150cm_75ad, m250cm_75ad, m350cm_75ad)
v3 <- cbind(m10cm_50ad, m50cm_50ad, m150cm_50ad, m250cm_50ad, m350cm_50ad)
v4 <- cbind(m10cm_25ad, m50cm_25ad, m150cm_25ad, m250cm_25ad, m350cm_25ad)
a <- cbind(10, 50, 150, 250, 350)

plot(a, v1, type = "b", xlab = "Number of catalysts", pch = 19,
     ylab = "Ratio of active/dormant molecules", ylim = c(0, 1.0))
lines(a, v2, pch = 19, type = "b", col = "red")
lines(a, v3, pch = 19, type = "b", col = "orange")
lines(a, v4, pch = 19, type = "b", col = "blue")

legend(150, 0.7, title = "Stochastic Transition \nto Dormancy",
       legend=c("90 %", "75 %", "50 %", "25 %"),  
       cex = 1.1,
       col = c("black", "red", "orange", "blue"),
       pch = 19,
       box.lty = 0,
       bg= rgb(1,1,1,0),
       y.intersp= 0.4)

# #Final pop plot
# fp10cm_90ad <- mean(a_10cm_90ad[900:1000]+d_10cm_90ad[900:1000])
# fp50cm_90ad <- mean(a_50cm_90ad[900:1000]+d_50cm_90ad[900:1000])
# fp150cm_90ad <- mean(a_150cm_90ad[900:1000]+d_150cm_90ad[900:1000])
# fp250cm_90ad <- mean(a_250cm_90ad[900:1000]+d_250cm_90ad[900:1000])
# fp350cm_90ad <- mean(a_350cm_90ad[900:1000]+d_350cm_90ad[900:1000])
# 
# fp10cm_75ad <- mean(a_10cm_75ad[900:1000]+d_10cm_75ad[900:1000])
# fp50cm_75ad <- mean(a_50cm_75ad[900:1000]+d_50cm_75ad[900:1000])
# fp150cm_75ad <- mean(a_150cm_75ad[900:1000]+d_150cm_75ad[900:1000])
# fp250cm_75ad <- mean(a_250cm_75ad[900:1000]+d_250cm_75ad[900:1000])
# fp350cm_75ad <- mean(a_350cm_75ad[900:1000]+d_350cm_75ad[900:1000])
# 
# fp10cm_50ad <- mean(a_10cm_50ad[900:1000]+d_10cm_50ad[900:1000])
# fp50cm_50ad <- mean(a_50cm_50ad[900:1000]+d_50cm_50ad[900:1000])
# fp150cm_50ad <- mean(a_150cm_50ad[900:1000]+d_150cm_50ad[900:1000])
# fp250cm_50ad <- mean(a_250cm_50ad[900:1000]+d_250cm_50ad[900:1000])
# fp350cm_50ad <- mean(a_350cm_50ad[900:1000]+d_350cm_50ad[900:1000])
# 
# fp10cm_25ad <- mean(a_10cm_25ad[900:1000]+d_10cm_25ad[900:1000])
# fp50cm_25ad <- mean(a_50cm_25ad[900:1000]+d_50cm_25ad[900:1000])
# fp150cm_25ad <- mean(a_150cm_25ad[900:1000]+d_150cm_25ad[900:1000])
# fp250cm_25ad <- mean(a_250cm_25ad[900:1000]+d_250cm_25ad[900:1000])
# fp350cm_25ad <- mean(a_350cm_25ad[900:1000]+d_350cm_25ad[900:1000])
# 
# p1 <- cbind(fp10cm_90ad, fp50cm_90ad, fp150cm_90ad, fp250cm_90ad, fp350cm_90ad)
# p2 <- cbind(fp10cm_75ad, fp50cm_75ad, fp150cm_75ad, fp250cm_75ad, fp350cm_75ad)
# p3 <- cbind(fp10cm_50ad, fp50cm_50ad, fp150cm_50ad, fp250cm_50ad, fp350cm_50ad)
# p4 <- cbind(fp10cm_25ad, fp50cm_25ad, fp150cm_25ad, fp250cm_25ad, fp350cm_25ad)
# 
# plot(a, p1, type = "b", xlab = "Number of catalysts", 
#      ylab = "Final Population (A)", ylim = c(0, 400))
# lines(a, p2, type = "b", col = "red")
# lines(a, p3, type = "b", col = "orange")
# lines(a, p4, type = "b", col = "blue")
# 
# legend(215, 300, title = "Stochastic Transition to Dormancy",
#        legend=c("90 %", "75 %", "50 %", "25 %"),  
#        fill = c("black","red", "orange", "blue"),
#        cex = 0.45) 

### Nice figure with ggplot
# require(ggplot2)
# 
# #stack vectors
# RA <- as.data.frame(t(rbind(a,v1,v2,v3,v4)))
# 
# p.mA <- ggplot()+ 
#   geom_point(data = RA, aes(x=V1, y=V2), color = "black", 
#              alpha = 0.5, pch = 16, size = 3)+
#   geom_line(data = RA, aes(x=V1, y=V2), color = "black", 
#              alpha = 0.5, pch = 16, size = 3)+
#   geom_point(data = RA, aes(x=V1, y=V3), color = "black", 
#              alpha = 0.5, pch = 16, size = 3)
#   
#   
# plot(p.mA)
  
  # scale_y_continuous(breaks=c(-2,-1,0,1,2,3),
  #                    labels = c(expression(0.01,0.1,1,10,100,1000)),
  #                    limits = c(-2, 3))+
  # scale_x_continuous(breaks= c(1000,1500,2000,2500),
  #                    labels = c(1000,1500,2000,2500),
  #                    limits = c(900, 2600))+
  # annotate(geom="text", y=-1.65, x=1400, 
  #          label=expression(paste('[U] = 10^(-0.0013*'*italic(z)~'+ 2.59)')),
  #          color="black", size=4)+
  # annotate(geom="text", y=-1.9, x=1400, 
  #          label=expression(paste(R^{2},'= 0.11')),
  #          color="black", size=4)+
  # labs(x="Elevation (m)", 
  #      y=expression(paste('[U] ('*italic("\u00B5")*'g/L)')))+
  # theme(panel.background = element_blank())


# 350cm; 25ad
IDM_350cm_25ad <- InDiMo(1000, 200, 10, 90, 5, 50, 50, 17, 6) 
#ts, ip, cm, ad, da, de, ds, br, ae

#model checks
a_350cm_25ad <- rep(NA, length(IDM_350cm_25ad[,1,1]))

for (i in 1:length(IDM_350cm_25ad[,1,1])) {
  a_350cm_25ad[i] <- length(which(IDM_350cm_25ad[i,,] == 1))
}

d_350cm_25ad <- rep(NA, length(IDM_350cm_25ad[,1,1]))

for (i in 1:length(IDM_350cm_25ad[,1,1])) {
  d_350cm_25ad[i] <- length(which(IDM_350cm_25ad[i,,] == 2))
}

plot(a_350cm_25ad) # 
plot(d_350cm_25ad) # 
plot(a_350cm_25ad/(a_350cm_25ad+d_350cm_25ad)) #prop active
plot(a_350cm_25ad+d_350cm_25ad) # total pop
