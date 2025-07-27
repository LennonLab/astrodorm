set.seed(42)  # For reproducibility

# Simulation parameters
steps <- 50
dt <- 1  # time step
k_background <- 0.001  # slow spontaneous rate
k_cat <- 0.05          # catalyzed rate

# Initialize data frame
state <- data.frame(
  time = 0:steps,
  A = rep(1.0, steps + 1),  # food, constant
  B = rep(1.0, steps + 1),  # food, constant
  C = rep(0.0, steps + 1),
  D = rep(0.0, steps + 1)
)

# Simulation loop
for (t in 1:steps) {
  A <- state$A[t]
  B <- state$B[t]
  C <- state$C[t]
  D <- state$D[t]
  
  # Reaction R1: A + B → C, catalyzed by C
  rate_R1 <- (k_background + k_cat * C) * A * B
  dC <- rate_R1 * dt
  
  # Reaction R2: A + C → D, catalyzed by D
  rate_R2 <- (k_background + k_cat * D) * A * C
  dD <- rate_R2 * dt
  
  # Update abundances
  state$A[t + 1] <- A  # remains constant
  state$B[t + 1] <- B  # remains constant
  state$C[t + 1] <- C + dC
  state$D[t + 1] <- D + dD
}

# View results
print(round(state, 3))

# Plot
matplot(state$time, state[, -1], type = "l", lwd = 2, lty = 1,
        col = c("darkgreen", "orange", "blue", "red"),
        ylab = "Abundance", xlab = "Time",
        main = "Minimal RAF Model with Continuous Abundances")
legend("topleft", legend = colnames(state)[-1],
       col = c("darkgreen", "orange", "blue", "red"),
       lty = 1, lwd = 2)

