################################################################################
#### DSEM MODELFIT ####
# run_simulation_all_models
# Authors: Prof. Holger Brandt, Mihai Falcusan, Kim Zierahn
# Date: 12/2024
################################################################################



# This file does the following:
# - Imports requirements and sets hyperparams for simulation run
# - Defines helper functions for simulation
# - Runs model simulations for desired conditions and saves results exp



################################################################################
#### 0. Importing relevant libraries and functions ####
################################################################################

library(mvtnorm)
library(matrixcalc)
library(lavaan)
library(blavaan)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Get the directory of the current script
current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # RStudio
# Alternatively, for non-interactive scripts:
# current_dir <- dirname(normalizePath(sys.frame(1)$ofile))

# Dynamically source files from the local directory
source(file.path(current_dir,  "lavaan_dsem_models_randomintercept.R"))
source(file.path(current_dir,  "gen_data_version03.R"))
source(file.path(current_dir,  "gen_data_version04.R"))
source(file.path(current_dir,  "lavaan_dsem_nullmodels.R"))
source(file.path(current_dir,  "fitfunctions_stan.R"))
source(file.path(current_dir,  "lavaan_dsem_models_randomintercept_tt.R"))
source(file.path(current_dir,  "lavaan_dsem_models_randomintercept_tt1.R"))

# working dir
setwd(current_dir)



################################################################################
#### 1. Setting hyperparameters ####
################################################################################

# Simulation params
# Save as filetype (one of local, RDS, csv)
save_as = 'csv'

# Number of cores for parallelisation
cores <- 1

# Number of people / time point?
Person_size <- c(91)
#Person_size <- c(31, 61, 91, 121, 151, 181, 211, 501, 1001, 1501, 2001, 2501)

# Number of measurement time points
Timepoints <- c(2) # N_timep
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

# frequentist fit indices
fitnom_lavaan <- c("npar","fmin","chisq","df","pvalue","baseline.chisq","baseline.df",
                   "baseline.pvalue","cfi","tli","nnfi","rfi","nfi","pnfi","ifi","rni",
                   "logl","unrestricted.logl","aic","bic","ntotal","bic2","rmsea","rmsea.ci.lower",
                   "rmsea.ci.upper","rmsea.ci.level","rmsea.pvalue","rmsea.close.h0","rmsea.notclose.pvalue",
                   "rmsea.notclose.h0","rmr","rmr_nomean","srmr","srmr_bentler","srmr_bentler_nomean",
                   "crmr","crmr_nomean","srmr_mplus","srmr_mplus_nomean","cn_05","cn_01","gfi",
                   "agfi","pgfi","mfi","ecvi") 

# bayesian fit indices
fitnom_blavaan <- c("npar", "PPP", "MargLogLik", 
                    "BRMSEA", "BGammaHat", "adjBGammaHat", 
                    "BMc", "BCFI", "BTLI", "BNFI", "chisq", "pd")

fitnom_stan <- c("BRMSEA", "BGammahat", "adjBGammahat", 
                 "BMc", "pd_1", "p_1", "loglik", "logliksat_1", "chisq")


# Parallelise cores 
# TODO: for core in cores
# temporary: using core=1 for 1 run
core <- 1



################################################################################
#### 2. Helper functions #####
################################################################################

#### Frequentist model ####

### LAVAAN ###
# Helper function to fit a single model, capture fit specific errors for later
# Input: data y0
# Output: model res1
fit_model_lav <- function(y0, N_timep){ 
  
  res1 <- try(sem(dsem[[N_timep]], data = y0, std.lv = TRUE), silent = TRUE)
  
  # if there is no try-error
  if(!inherits(res1, "try-error")){ 
    # if model converged: save fit indices
    if(res1@optim$converged == T){ 
      fitmeasures(res1)  
    }
    else { # if model did not converge: save NAs so that simulation doesnt get broken
      print("noconvergence")
      return(rep(NA, length(fitnom_lavaan)))
    }
  } 
  else { # if there is a try-error: save NAs so that simulation doesnt get broken
    print("tryerror")
    return(rep(NA, length(fitnom_lavaan)))
  }
}


fit_model_lav_true <- function(y0, N_timep, type_MISS){ 
  
  if(type_MISS=="tt"){
    res1 <- try(sem(dsem.tt[[N_timep]], data = y0, std.lv = TRUE, se = "none"), silent = TRUE)
  }else{
    res1 <- try(sem(dsem.tt1[[N_timep]], data = y0, std.lv = TRUE, se = "none"), silent = TRUE)
  }
  
  # if there is no try-error
  if(!inherits(res1, "try-error")){ 
    # if model converged: save fit indices
    if(res1@optim$converged == T){ 
      fitmeasures(res1)  
    }
    else { # if model did not converge: save NAs so that simulation doesnt get broken
      print("noconvergence")
      return(rep(NA, length(fitnom_lavaan)))
    }
  } 
  else { # if there is a try-error: save NAs so that simulation doesnt get broken
    print("tryerror")
    return(rep(NA, length(fitnom_lavaan)))
  }
}


#### Bayesian models ####

### BLAVAAN ###
fit_model_blav <- function(y0, N_timep){ 
  
  res1 <- try(bsem(dsem[[N_timep]], data = y0, std.lv = TRUE,
                   n.chains = 4, burnin = 1000, sample = 1000), silent = TRUE)
  
  nullmodel <- try(bsem(null_model(N_timep), data = y0, 
                        n.chains = 4, burnin = 1000, sample = 1000), silent = TRUE)
  
  # if there is no try-error
  if(!inherits(res1, "try-error")){ 
    # if model converged: save fit indices
    if(blavInspect(res1, "converged") == TRUE){ 
      fit_indices <- blavFitIndices(res1, baseline.model = nullmodel) 
      BRMSEA <- mean(fit_indices@indices$BRMSEA)
      BGammaHat <- mean(fit_indices@indices$BGammaHat)
      adjBGammaHat <- mean(fit_indices@indices$adjBGammaHat)
      BMc <- mean(fit_indices@indices$BMc)
      BCFI <- mean(fit_indices@indices$BCFI)
      BTLI <- mean(fit_indices@indices$BTLI)
      BNFI <- mean(fit_indices@indices$BNFI)
      npar <- blavInspect(res1, "npar")
      marg_loglik <- blavInspect(res1, "test")[[1]]$stat
      ppp <- blavInspect(res1, "test")[[2]]$stat
      chisq <- mean(fit_indices@details$chisq)
      pd <- fit_indices@details$pD
      
      # store
      modelfit <- c(npar, ppp, marg_loglik, BRMSEA, BGammaHat,
                    adjBGammaHat, BMc, BCFI, BTLI, BNFI, chisq, pd)
      fitresults <- matrix(modelfit, nrow = 1)
      colnames(fitresults) <- fitnom_blavaan
      
      return(fitresults)
    }
    else { # if model did not converge: save NAs so that simulation doesnt get broken
      print("noconvergence")
      return(rep(NA, length(fitnom_blavaan)))
    }
  } 
  else { # if there is a try-error: save NAs so that simulation doesnt get broken
    print("tryerror")
    return(rep(NA, length(fitnom_blavaan)))
  }
} 


### STAN ###
fit_model_stan <- function(ydat, ydat2, timepoints, person_size){ 
  
  # load and compile model
  rt1 <- try(stanc("model3_dsem.stan"))
  sm1 <- try(stan_model(stanc_ret = rt1, verbose=FALSE)) 
  
  # parameter names
  params2 <- c("ly", "sigmaeps", "sigmaeta1", "sigmaeta2", "ty", "ka", "beta1")
  
  data.stan1 <- list(N = person_size, Nt = timepoints, y = ydat)
  fit1 <- try(sampling(sm1, data = data.stan1))
  
  
  # if there is no try-error
  if(!inherits(fit1, "try-error")){ 
    # save fit indices
    modelfit <- getfit(fit1, ydat, ydat2)
    fitresults1 <- matrix(modelfit, nrow = 1)
    colnames(fitresults1) <- fitnom_stan
    return(fitresults1)
  } 
  else { # if there is a try-error: save NAs so that simulation doesnt get broken
    print("tryerror")
    rep(NA, length(fitnom_stan))
  }
} 


#### Other functions ####

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



################################################################################
#### 3. Simulation #####
################################################################################

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
          fitm_lavaan <- matrix(NA, N_sim_samples, length(fitnom_lavaan))
          fitm_lavaan_true <- matrix(NA, N_sim_samples, length(fitnom_lavaan))
          fitm_blavaan <- matrix(NA, N_sim_samples, length(fitnom_blavaan))
          fitm_stan <- matrix(NA, N_sim_samples, length(fitnom_stan))
          
          
          # Loop over samples
          for (Index_sample in 1:N_sim_samples) {
            temp_dat <- generate_temp_data(Type_misfit, N_pers, N_timep, phi0, mu0, ar0, ly0, ly1, td)
            ydat <- temp_dat$y # use array version of data for stan
            ydat2 <- temp_dat$y0 # use the transformed version in the getfit function
            
            
            # Check if data is positive definite
            if (temp_dat[["is_positive_def"]]) {
              y0 <- temp_dat[["y0"]]
              print("Fit model and store results: lavaan")
              fitm_lavaan[Index_sample, ] <- fit_model_lav(ydat2, N_timep) # Fit model and store results
              fitm_lavaan_true[Index_sample, ] <- fit_model_lav_true(ydat2, N_timep, Type_misfit) # Fit true model and store results (for plotting)
              print("Fit model and store results: blavaan")
              fitm_blavaan[Index_sample, ] <- fit_model_blav(ydat2, N_timep) # Fit model and store results
              print("Fit model and store results: stan")
              fitm_stan[Index_sample, ] <- as.numeric(fit_model_stan(ydat, ydat2, N_timep, N_pers)) # Fit model and store results
              
            }
          }
          
          # Print model being saved after simulation block
          cat('Ran misfit_type:', Type_misfit, ' misfit_size:', Size_misfit, ' for', N_sim_samples, 'indep samples', '\n')
          # Save results
          colnames(fitm_lavaan) <- fitnom_lavaan
          colnames(fitm_lavaan_true) <- fitnom_lavaan
          colnames(fitm_blavaan) <- fitnom_blavaan
          colnames(fitm_stan) <- fitnom_stan
          
          
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
            if (!exists("lavaan_simulation_results_df")) {
              lavaan_simulation_results_df <- data.frame()
            }
            
            if (!exists("lavaan_true_simulation_results_df")) {
              lavaan_true_simulation_results_df <- data.frame()
            }
            
            if (!exists("blavaan_simulation_results_df")) {
              blavaan_simulation_results_df <- data.frame()
            }
            
            if (!exists("stan_simulation_results_df")) {
              stan_simulation_results_df <- data.frame()
            }
            
            
            # Add the current `fitm_lavaan` as a new row (ensure it can be coerced into a dataframe)
            if (is.data.frame(fitm_lavaan)) {
              fitm_lavaan_row <- fitm_lavaan
            } else {
              fitm_lavaan_row <- as.data.frame(fitm_lavaan)
            }
            
            if (is.data.frame(fitm_lavaan_true)) {
              fitm_lavaan_true_row <- fitm_lavaan_true
            } else {
              fitm_lavaan_true_row <- as.data.frame(fitm_lavaan_true)
            }
            
            if (is.data.frame(fitm_blavaan)) {
              fitm_blavaan_row <- fitm_blavaan
            } else {
              fitm_blavaan_row <- as.data.frame(fitm_blavaan)
            }
            
            if (is.data.frame(fitm_stan)) {
              fitm_stan_row <- (fitm_stan)
            } else {
              fitm_stan_row <- as.data.frame((fitm_stan))
            }
            
            # Optionally add simulation metadata (e.g., Exp_name_info) to the row
            fitm_lavaan_row$Exp_name_info <- Exp_name_info
            fitm_lavaan_true_row$Exp_name_info <- Exp_name_info
            fitm_blavaan_row$Exp_name_info <- Exp_name_info
            fitm_stan_row$Exp_name_info <- Exp_name_info
            
            # Append the row to the dataframe
            lavaan_simulation_results_df <- rbind(lavaan_simulation_results_df, fitm_lavaan_row)
            lavaan_true_simulation_results_df <- rbind(lavaan_true_simulation_results_df, fitm_lavaan_true_row)
            blavaan_simulation_results_df <- rbind(blavaan_simulation_results_df, fitm_blavaan_row)
            stan_simulation_results_df <- rbind(stan_simulation_results_df, fitm_stan_row)
            
          }else if(save_as == 'RDS'){
            # Save your RDS file in the desired directory with EXP_NAME in the filename
            saveRDS(fitm_lavaan, file = file.path(save_dir, paste0("lavaan_", Exp_name_info,"_time_",Exp_time, ".RDS")))
            saveRDS(fitm_lavaan_true, file = file.path(save_dir, paste0("lavaan_true_", Exp_name_info,"_time_",Exp_time, ".RDS")))
            saveRDS(fitm_blavaan, file = file.path(save_dir, paste0("blavaan_", Exp_name_info,"_time_",Exp_time, ".RDS")))
            saveRDS(fitm_stan, file = file.path(save_dir, paste0("stan_", Exp_name_info,"_time_",Exp_time, ".RDS")))
          }else if(save_as == 'csv'){
            csv_path_lavaan <- file.path(save_dir, paste0("lavaan_", Exp_name_info, ".csv"))
            write.csv(fitm_lavaan, file = csv_path_lavaan, row.names = FALSE) 
            csv_path_lavaan_true <- file.path(save_dir, paste0("lavaan_true_", Exp_name_info, ".csv"))
            write.csv(fitm_lavaan_true, file = csv_path_lavaan_true, row.names = FALSE) 
            csv_path_blavaan <- file.path(save_dir, paste0("blavaan_", Exp_name_info, ".csv"))
            write.csv(fitm_blavaan, file = csv_path_blavaan, row.names = FALSE)
            csv_path_stan <- file.path(save_dir, paste0("stan_", Exp_name_info, ".csv"))
            write.csv(fitm_stan, file = csv_path_stan, row.names = FALSE)
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


