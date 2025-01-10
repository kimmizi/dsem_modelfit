# This file computes 10 simulation configurations that have roughly the same runtime
# each one is computed based on a set of nt x np conditions.
# in the final config, each one of the 10 will have the same number of models runnign
# that is 1 normal, 4 misfits. both in lavaan and blavaan (no stan)

# Define the vectors
Timepoints <- c(1, 2, 3, 4, 5, 10, 15)
Person_size <- c(31, 61, 91, 121, 151, 181, 211, 501, 1001, 1501, 2001)

# Create the full matrix
full_matrix <- outer(Timepoints, Person_size, "*")

# Create triangular matrix (upper triangle)
tri_matrix <- full_matrix
tri_matrix[lower.tri(tri_matrix)] <- NA

# Create a data frame of all workloads with their indices
workloads <- data.frame(
  timepoint_idx = row(tri_matrix)[!is.na(tri_matrix)],
  person_idx = col(tri_matrix)[!is.na(tri_matrix)],
  workload = tri_matrix[!is.na(tri_matrix)]
)

# Add actual timepoint and person size values for reference
workloads$timepoint <- Timepoints[workloads$timepoint_idx]
workloads$person_size <- Person_size[workloads$person_idx]

# Sort workloads in descending order
workloads <- workloads[order(-workloads$workload), ]

# Initialize 10 simulation groups
n_sims <- 10
sim_groups <- vector("list", n_sims)
sim_totals <- numeric(n_sims)

# Distribute workloads using a greedy algorithm
for(i in 1:nrow(workloads)) {
  # Find simulation group with minimum total workload
  min_group <- which.min(sim_totals)
  
  # Add workload to that group
  sim_groups[[min_group]] <- rbind(sim_groups[[min_group]], 
                                   workloads[i, ])
  sim_totals[min_group] <- sim_totals[min_group] + workloads$workload[i]
}

# Print summary of distribution
cat("Summary of workload distribution:\n")
for(i in 1:n_sims) {
  cat(sprintf("\nSimulation %d (Total workload: %d):\n", i, sim_totals[i]))
  group_data <- sim_groups[[i]]
  for(j in 1:nrow(group_data)) {
    cat(sprintf("  T=%d, P=%d (workload: %d)\n", 
                group_data$timepoint[j], 
                group_data$person_size[j], 
                group_data$workload[j]))
  }
}

# Print statistics about the distribution
cat("\nDistribution statistics:\n")
cat(sprintf("Mean workload per simulation: %.2f\n", mean(sim_totals)))
cat(sprintf("Standard deviation: %.2f\n", sd(sim_totals)))
cat(sprintf("Coefficient of variation: %.2f%%\n", 100 * sd(sim_totals) / mean(sim_totals)))
cat(sprintf("Max/Min ratio: %.2f\n", max(sim_totals) / min(sim_totals)))