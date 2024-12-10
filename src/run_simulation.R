# This file does the following:
# - Imports requirements and sets hyperparams for simulation run
# - Defines helper functions for simulation
# - Runs model simulations for desired conditions and saves results exp

# 0. Importing relevant libraries and functions ################################

library(mvtnorm)
library(matrixcalc)
library(lavaan)
# any other?

# Get the directory of the current script
current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # RStudio
# Alternatively, for non-interactive scripts:
# current_dir <- dirname(normalizePath(sys.frame(1)$ofile))

# Dynamically source files from the local directory
source(file.path(current_dir, "lavaan_dsem_models_randomintercept.R"))
source(file.path(current_dir,  "gen_data_version03.R"))

# source(file.path(current_dir, "lavaan_dsem_models_comp1.R"))
# source(file.path(current_dir,  "lavaan_dsem_models_comp2.R"))
source(file.path(current_dir,  "gen_data_version04.R"))

# working dir?
#setwd()



# 1. Setting hyperparameters ###################################################

# Simulation params
# Save as filetype (one of local, RDS, csv)
save_as = 'local'

# Number of cores for parallelisation
cores <- 1

# Number of people / time point?
Person_size <- c(91)
#Person_size <- c(91, 121, 151, 181, 211, 501, 1001, 1501, 2001, 2501, 61, 31)

# Number of measurement time points
Timepoints <- 5 # N_timep
#Timepoints <- c(1:5, 10, 15) # N_timep

# Size of crossloading (indicating misfit)
Size_crossloading <- c(0,.3,.6)

# What models to specify and run: within time points/ between time points
Type_crossloading <- c("none","tt","tt1")

# Number of data samples to generate and run model fit on them
N_sim_samples <- 3

# DSEM model params
phi0 <- diag(2)*.7+.3  # cov(eta)
mu0  <- c(0, 0)        # mean(eta)
ar0  <- c(.3, .3)      # ar(1) structure
ly00 <- .6
ly0  <- matrix(c(ly00, ly00, ly00, 0, 0, 0,
                 0, 0, 0, ly00, ly00, ly00), 6, 2, byrow = F) # factor loadings
td   <- diag(6)*(1 - ly00^2) # cond. var(y) -> res var

# Fit indices
fitnom <- c("npar","fmin","chisq","df","pvalue","baseline.chisq","baseline.df",
            "baseline.pvalue","cfi","tli","nnfi","rfi","nfi","pnfi","ifi","rni",
            "logl","unrestricted.logl","aic","bic","ntotal","bic2","rmsea","rmsea.ci.lower",
            "rmsea.ci.upper","rmsea.ci.level","rmsea.pvalue","rmsea.close.h0","rmsea.notclose.pvalue",
            "rmsea.notclose.h0","rmr","rmr_nomean","srmr","srmr_bentler","srmr_bentler_nomean",
            "crmr","crmr_nomean","srmr_mplus","srmr_mplus_nomean","cn_05","cn_01","gfi",
            "agfi","pgfi","mfi","ecvi") 


# Parallelise cores 
# TODO: for core in cores
# temporary: using core=1 for 1 run
core <- 1




# 2. Helper functions ##########################################################

# Helper function to fit a single model, capture fit specific errors for later
# Input: data y0
# Output: model res1
fit_model <- function(y0){ 
  
  res1 <- try(sem(dsem[[N_timep]], data = y0, std.lv = TRUE), silent = TRUE)
  # if there is no try-error
  if(!inherits(res1, "try-error")){ 
    # if model converged: save fit indices
    if(res1@optim$converged == T){ 
      fitmeasures(res1)  
    }
    else { # if model did not converge: save NAs so that simulation doesnt get broken
      rep(NA, length(fitnom))
    }
  } 
  else { # if there is a try-error: save NAs so that simulation doesnt get broken
    rep(NA, length(fitnom))
  }
}


# Helper function to initialize factor loadings
initialize_ly1 <- function(Type_misfit, Size_misfit) {
  if (Type_misfit == "tt") {
    return(matrix(c(0, 0, 0, Size_misfit, 0, 0,
                    0, 0, 0, 0, 0, 0), 6, 2, byrow = F))
  } else if (Type_misfit == "tt1") {
    return(matrix(c(Size_misfit, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0), 6, 2, byrow = F))
  } else {
    return(matrix(c(0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0), 6, 2, byrow = F)) # For "none"
  }
}

# Helper function to initialize simulation specific data set

generate_temp_data <- function(Type_crossloading, N_pers, N_timep, phi0, mu0, ar0, ly0, ly1, td) {
  if (Type_crossloading == "tt") {
    return(gendata01(N_pers, N_timep, phi0, mu0, ar0, ly0, ly1, td))
  } else if (Type_crossloading == "tt1") {
    return(gendata02(N_pers, N_timep, phi0, mu0, ar0, ly0, ly1, td))
  } else {
    return(gendata01(N_pers, N_timep, phi0, mu0, ar0, ly0, ly1, td))
  }
}


# 3. Simulation ################################################################


# Record loop start time
start_measurement_time <- proc.time()

# Main simulation loop
for (N_pers in Person_size) {
  for (N_timep in Timepoints) {
    # Check specification requirement
    if (N_timep * 6 >= N_pers) next
    for (Type_misfit in Type_crossloading) {

      # Set random seed for reproducibility
      set.seed(131212023 + core + N_pers + N_timep * 1000 + as.numeric(as.factor(Type_misfit)) * 10000)

      for (Size_misfit in Size_crossloading) {
        # if type tt or tt1, dont run size=0, else run
        if(Size_misfit == 0 && Type_misfit != "none"){
          cat('skipped over misfit type:', Type_misfit, " and size:", Size_misfit, "\n")
        }else if(Size_misfit != 0 && Type_misfit == "none"){
          cat('skipped over misfit type:', Type_misfit, " and size:", Size_misfit, "\n")
        }else{
          # Initialize factor loadings for the current condition
          ly1 <- initialize_ly1(Type_misfit, Size_misfit)
          #print(Type_misfit)
          #print(Size_misfit)

          # Prepare file name
          Exp_name_info <- paste(N_pers, N_timep, Type_misfit,
                                            Size_misfit, core, "_version03_rand", sep = "_")

          # Initialize empty matrix to store fit indices
          fitm1 <- matrix(NA, N_sim_samples, length(fitnom))

          # Loop over samples
          for (Index_sample in 1:N_sim_samples) {
            temp_dat <- generate_temp_data(Type_misfit, N_pers, N_timep, phi0, mu0, ar0, ly0, ly1, td)

            # Check if data is positive definite
            if (temp_dat[["is_positive_def"]]) {
              y0 <- temp_dat[["y0"]]
              fitm1[Index_sample, ] <- fit_model(y0) # Fit model and store results
            }
          }

          # Print model being saved after simulation block
          cat('Ran misfit_type:', Type_misfit, ' misfit_size:', Size_misfit, ' for', N_sim_samples, 'indep samples', '\n')
          # Save results
          colnames(fitm1) <- fitnom

          # Move up one directory level to dsem_modelfit
          parent_dir <- dirname(current_dir)
          # Define the experiment name
          Exp_day <- format(Sys.time(), "%Y-%m-%d") # Format: YYYY-MM-DD_HH-MM-SS
          Exp_time <- format(Sys.time(), "%H-%M-%S") # Format: YYYY-MM-DD_HH-MM-SS
          # Define the save directory in dsem_modelfit/exp
          save_dir <- file.path(parent_dir, "exp", Exp_day)
          # Create the day directory if it does not exist
          if (!dir.exists(save_dir)) {
            dir.create(save_dir, recursive = TRUE)
          }
          # Save model along its time of saving to have tractability
          
          if(save_as == 'local'){
            # Initialize an empty dataframe to store results across runs
            if (!exists("simulation_results_df")) {
              simulation_results_df <- data.frame()
            }
            
            # Add the current `fitm1` as a new row (ensure it can be coerced into a dataframe)
            if (is.data.frame(fitm1)) {
              fitm1_row <- fitm1
            } else {
              fitm1_row <- as.data.frame(fitm1)
            }
            
            # Optionally add simulation metadata (e.g., Exp_name_info) to the row
            fitm1_row$Exp_name_info <- Exp_name_info
            
            # Append the row to the dataframe
            simulation_results_df <- rbind(simulation_results_df, fitm1_row)
          }else if(save_as == 'RDS'){
            # Save your RDS file in the desired directory with EXP_NAME in the filename
            saveRDS(fitm1, file = file.path(save_dir, paste0(Exp_name_info,"_time_",Exp_time, ".RDS")))
          }else if(save_as == 'csv'){
            csv_path <- file.path(save_dir, paste0(Exp_name_info, "_time_", Exp_time, ".csv"))
            write.csv(fitm1, file = csv_path, row.names = FALSE) 
          }else{
            cat("Format of type", save_as, " cannot be saved.")
          }
        }
      }
    }
  }
}

# Record end time
end_measurement_time <- proc.time()

# Calculate elapsed time
time_taken <- end_measurement_time - start_measurement_time
print(time_taken)


