################################################################################
#### DSEM MODELFIT ####
# run_simulation_all_models
# Authors: Prof. Holger Brandt, Mihai Falcusan, Kim Zierahn
# Date: 01/2025
################################################################################



# This file does the following:
# - Imports requirements and sets hyperparams for simulation run based on simulation nr
# - Manually defines 10 simulation conditions based on simulation_workload_calculator.R
# - Runs model simulations for desired conditions and saves results exp

################################################################################
#### 1. Read system configurations from sbatch file
################################################################################

# find cluster working directory
current_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  script_path <- NULL
  
  for (arg in args) {
    if (startsWith(arg, file_flag)) {
      script_path <- sub(file_flag, "", arg)
      break
    }
  }
  
  if (!is.null(script_path)) {
    dirname(normalizePath(script_path))
  } else {
    stop("Could not determine script path")
  }
}, error = function(e) {
  getwd()  # Fallback to working directory
})

print(current_dir)
#current_dir = "C:/Users/mihai/Documents/Faculta/Research Project/dsem_modelfit/src"

# Automatically load all libraries only once ###################################
source(file.path(current_dir,  "load_libraries.R"))

# load all other required files ################################################
source(file.path(current_dir,  "lavaan_dsem_models_randomintercept.R"))
source(file.path(current_dir,  "gen_data_version03.R"))
source(file.path(current_dir,  "gen_data_version04.R"))
source(file.path(current_dir,  "lavaan_dsem_nullmodels.R"))
#source(file.path(current_dir,  "fitfunctions_stan.R"))
source(file.path(current_dir,  "lavaan_dsem_models_randomintercept_tt.R"))
source(file.path(current_dir,  "lavaan_dsem_models_randomintercept_tt1.R"))
source(file.path(current_dir,  "helper_functions_sim.R"))

# Only thing that doesn't need to be parallelised? Hyperparams?

# detect cores and split tasks


cores <- detectCores()
cat('detected cores:', cores)
registerDoParallel(cores) 
# from here on, each core does 1 computation of this
sim_exp <- foreach (core=1:cores) %dopar% {
  ################################################################################
  #### 0. Importing relevant libraries and functions ####
  ################################################################################

  options(mc.cores = parallel::detectCores())
  # working dir
  setwd(current_dir)
  
  ################################################################################
  #### 1. Setting hyperparameters ####
  ################################################################################
  
  # Simulation params
  # Save as filetype (one of local, RDS, csv)
  save_as = 'csv'
  
  # Number of people / time point?
  
  Person_size <- c(31, 61, 91, 121, 151, 181, 211, 501, 1001, 1501, 2001)
  
  # Number of measurement time points
  
  #Timepoints <- c(1:5, 10, 15) # N_timep
  #Timepoints <- c(1,2,3,4,5,10,15)
  
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
                      "BRMSEA", "BGammaHat", "adjBGammaHat", "BMc", 
                      "BCFI_0C", "BTLI_0C", "BNFI_0C",
                      "BCFI_0A", "BTLI_0A", "BNFI_0A",
                      "chisq", "pd")
    
    fitnom_stan <- c("BRMSEA", "BGammahat", "adjBGammahat", 
                    "BMc", "pd_1", "p_1", "loglik", "logliksat_1", "chisq")
    
    
    # temporary: using core=1 for 1 run
    # core <- 1
    
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
            fitm_blavaan <- matrix(NA, N_sim_samples, length(fitnom_blavaan))
            #fitm_stan <- matrix(NA, N_sim_samples, length(fitnom_stan))
            
            
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
                print("Fit model and store results: blavaan")
                fitm_blavaan[Index_sample, ] <- fit_model_blav(ydat2, N_timep) # Fit model and store results
                #print("Fit model and store results: stan")
                #fitm_stan[Index_sample, ] <- as.numeric(fit_model_stan(ydat, ydat2, N_timep, N_pers)) # Fit model and store results
                
              }
            }
            
            # Print model being saved after simulation block
            cat('Ran misfit_type:', Type_misfit, ' misfit_size:', Size_misfit, ' for', N_sim_samples, 'indep samples', '\n')
            # Save results
            colnames(fitm_lavaan) <- fitnom_lavaan
            colnames(fitm_blavaan) <- fitnom_blavaan
            #colnames(fitm_stan) <- fitnom_stan
            
            
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

              if (!exists("blavaan_simulation_results_df")) {
                blavaan_simulation_results_df <- data.frame()
              }
              
              #if (!exists("stan_simulation_results_df")) {
              #  stan_simulation_results_df <- data.frame()
              #}
              
              
              # Add the current `fitm_lavaan` as a new row (ensure it can be coerced into a dataframe)
              if (is.data.frame(fitm_lavaan)) {
                fitm_lavaan_row <- fitm_lavaan
              } else {
                fitm_lavaan_row <- as.data.frame(fitm_lavaan)
              }
              
              if (is.data.frame(fitm_blavaan)) {
                fitm_blavaan_row <- fitm_blavaan
              } else {
                fitm_blavaan_row <- as.data.frame(fitm_blavaan)
              }
              
              #if (is.data.frame(fitm_stan)) {
              #  fitm_stan_row <- (fitm_stan)
              #} else {
              #  fitm_stan_row <- as.data.frame((fitm_stan))
              #}
              
              # Optionally add simulation metadata (e.g., Exp_name_info) to the row
              fitm_lavaan_row$Exp_name_info <- Exp_name_info
              fitm_blavaan_row$Exp_name_info <- Exp_name_info
              #fitm_stan_row$Exp_name_info <- Exp_name_info
              
              # Append the row to the dataframe
              lavaan_simulation_results_df <- rbind(lavaan_simulation_results_df, fitm_lavaan_row)
              blavaan_simulation_results_df <- rbind(blavaan_simulation_results_df, fitm_blavaan_row)
              #stan_simulation_results_df <- rbind(stan_simulation_results_df, fitm_stan_row)
              
            }else if(save_as == 'RDS'){
              # Save your RDS file in the desired directory with EXP_NAME in the filename
              saveRDS(fitm_lavaan, file = file.path(save_dir, paste0("lavaan_", Exp_name_info,"_time_",Exp_time, ".RDS")))
              saveRDS(fitm_blavaan, file = file.path(save_dir, paste0("blavaan_", Exp_name_info,"_time_",Exp_time, ".RDS")))
              #saveRDS(fitm_stan, file = file.path(save_dir, paste0("stan_", Exp_name_info,"_time_",Exp_time, ".RDS")))
            }else if(save_as == 'csv'){
              csv_path_lavaan <- file.path(save_dir, paste0("lavaan_", Exp_name_info, ".csv"))
              write.csv(fitm_lavaan, file = csv_path_lavaan, row.names = FALSE) 
              csv_path_blavaan <- file.path(save_dir, paste0("blavaan_", Exp_name_info, ".csv"))
              write.csv(fitm_blavaan, file = csv_path_blavaan, row.names = FALSE)
              #csv_path_stan <- file.path(save_dir, paste0("stan_", Exp_name_info, ".csv"))
              #write.csv(fitm_stan, file = csv_path_stan, row.names = FALSE)
            }else{
              cat("Format of type", save_as, " cannot be saved.")
            }
          }
        }
      }
    }
  } # end whole simulation

  # Record end time
  end_measurement_time <- proc.time()

  # Calculate elapsed time
  time_taken <- end_measurement_time - start_measurement_time
  cat("Time taken on core=", core, " was ",time_taken)
  
} # end loop over cores