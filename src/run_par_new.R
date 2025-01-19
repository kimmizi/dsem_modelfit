# run_par_new.R

# Plan: 
# 1. Initialise everything inside the file so that there are no problems with the parallel workers (this should work indep of part 2)
# 2. Each file gets a seed from the slurm bash file
# 3. For each seed, we have a new gendata, new sim
# 4. Somehow saves each run in a separate file

# 1.1. Packages and file dependencies #######################################################################

# all required packages should be here and initialised only once at the start
library(mvtnorm)
library(matrixcalc)
library(blavaan)
library(rstan)
# I think this is so that rstan and blavaan can operate efficiently. 
# How might this interact with the other part of the simulation?
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# 1.2. DSEM models ###########################################################################################

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
# We try to read in all required model files
source(file.path(current_dir,  "lavaan_dsem_models_randomintercept.R"))
source(file.path(current_dir,  "lavaan_dsem_models_randomintercept_tt.R"))
source(file.path(current_dir,  "lavaan_dsem_models_randomintercept_tt1.R"))
source(file.path(current_dir,  "lavaan_dsem_nullmodels.R"))

# 1.3. Workload specifications ###############################################################################

# Define all workloads directly in R. these were computed via simulation_workload_calculator.R
# to have an approx. equal workload. when init one of 10 workloads via bash argument,
# the R file reads one of these and simulates all the N_t and N_p it contains

# so a workload is basically a list of dsem model conditions to run. We suspect running N_t=2001 x  N_p =15 is equal to running a few smaller conditions
workloads <- list(
  sim_1 = list(
    total_workload = 30015,
    combinations = data.frame(
      timepoint = c(15),
      person_size = c(2001)
    )
  ),
  sim_2 = list(
    total_workload = 22515,
    combinations = data.frame(
      timepoint = c(15),
      person_size = c(1501)
    )
  ),
  sim_3 = list(
    total_workload = 20616,
    combinations = data.frame(
      timepoint = c(10, 3, 2, 1),
      person_size = c(2001, 121, 91, 61)
    )
  ),
  sim_4 = list(
    total_workload = 20593,
    combinations = data.frame(
      timepoint = c(15, 10, 1, 4, 1, 2),
      person_size = c(1001, 211, 2001, 181, 501, 121)
    )
  ),
  sim_5 = list(
    total_workload = 20576,
    combinations = data.frame(
      timepoint = c(10, 2, 5, 5, 3, 1),
      person_size = c(1501, 1501, 211, 151, 181, 211)
    )
  ),
  sim_6 = list(
    total_workload = 20623,
    combinations = data.frame(
      timepoint = c(10, 3, 3, 1, 1, 2, 2, 1),
      person_size = c(1001, 1501, 1001, 1501, 1001, 181, 61, 121)
    )
  ),
  sim_7 = list(
    total_workload = 20588,
    combinations = data.frame(
      timepoint = c(5, 5, 5, 3, 5, 3, 1, 1),
      person_size = c(2001, 1001, 501, 501, 181, 151, 181, 31)
    )
  ),
  sim_8 = list(
    total_workload = 20621,
    combinations = data.frame(
      timepoint = c(4, 10, 4, 2, 4, 4, 3),
      person_size = c(2001, 501, 1001, 1001, 211, 121, 91)
    )
  ),
  sim_9 = list(
    total_workload = 20595,
    combinations = data.frame(
      timepoint = c(15, 3, 15, 4, 2, 4, 2),
      person_size = c(501, 2001, 211, 501, 501, 151, 151)
    )
  ),
  sim_10 = list(
    total_workload = 20618,
    combinations = data.frame(
      timepoint = c(5, 4, 2, 10, 3, 2, 1, 1),
      person_size = c(1501, 1501, 2001, 181, 211, 211, 151, 91)
    )
  )
)

# 1.4. Helper functions ################################################################################

# Function to initialize factor loadings
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

# Function to initialize simulation specific data set
generate_temp_data <- function(Type_crossloading, N_pers, N_timep, phi0, mu0, ar0, ly0, ly1, td) {
  if (Type_crossloading == "tt") {
    return(gendata01(N_pers, N_timep, phi0, mu0, ar0, ly0, ly1, td))
  } else if (Type_crossloading == "tt1") {
    return(gendata02(N_pers, N_timep, phi0, mu0, ar0, ly0, ly1, td))
  } else if (Type_crossloading == "none"){
    return(gendata01(N_pers, N_timep, phi0, mu0, ar0, ly0, ly1, td))
  } else {
    stop("Invalid misfit type. Must be one of tt, tt1, none")
  }
}

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
      return(rep(NA, length(fitnom_lavaan)))
    }
  } 
  else { # if there is a try-error: save NAs so that simulation doesnt get broken
    return(rep(NA, length(fitnom_lavaan)))
  }
}


fit_model_blav <- function(y0, N_timep){ 
  
  # At the start of the function
  cat("Starting fit_model_blav with N_timep =", N_timep, "\n")
  print(str(y0))  # Check input data structure
  # First check if null models fit successfully
  print("nullmodel_0A")
  nullmodel_0A <- try(bsem(null_model_0A(N_timep), data = y0, 
                           n.chains = 4, burnin = 1000, sample = 1000), silent = TRUE)
  
  if(inherits(nullmodel_0A, "try-error")) {
    cat("Null model 0A failed to fit\n")
    return(rep(NA, length(fitnom_blavaan)))
  }
  
  print("nullmodel_0C")
  nullmodel_0C <- try(bsem(null_model_0C(N_timep), data = y0, 
                           n.chains = 4, burnin = 1000, sample = 1000), silent = TRUE)
  
  if(inherits(nullmodel_0C, "try-error")) {
    cat("Null model 0C failed to fit\n")
    return(rep(NA, length(fitnom_blavaan)))
  }
  
  print("res1")
  res1 <- try(bsem(dsem[[N_timep]], data = y0, std.lv = TRUE,
                   n.chains = 4, burnin = 1000, sample = 1000), silent = TRUE)
  
  if(!inherits(res1, "try-error")){ 
    if(blavInspect(res1, "converged") == TRUE){ 
      # Add error checking for fit indices
      tryCatch({
        fit_indices_C <- blavFitIndices(res1, baseline.model = nullmodel_0C)
        fit_indices_A <- blavFitIndices(res1, baseline.model = nullmodel_0A)
        
        BCFI_0C <- mean(fit_indices_C@indices$BCFI)
        BTLI_0C <- mean(fit_indices_C@indices$BTLI)
        BNFI_0C <- mean(fit_indices_C@indices$BNFI)
        
        BCFI_0A <- mean(fit_indices_A@indices$BCFI)
        BTLI_0A <- mean(fit_indices_A@indices$BTLI)
        BNFI_0A <- mean(fit_indices_A@indices$BNFI)
        
        # remaining fit indices (using either one, they should be the same)
        BRMSEA <- mean(fit_indices_A@indices$BRMSEA)
        BGammaHat <- mean(fit_indices_A@indices$BGammaHat)
        adjBGammaHat <- mean(fit_indices_A@indices$adjBGammaHat)
        BMc <- mean(fit_indices_A@indices$BMc)
        
        npar <- blavInspect(res1, "npar")
        marg_loglik <- blavInspect(res1, "test")[[1]]$stat
        ppp <- blavInspect(res1, "test")[[2]]$stat
        chisq <- mean(fit_indices_A@details$chisq)
        pd <- fit_indices_A@details$pD
        
        modelfit <- c(npar, ppp, marg_loglik, 
                      BRMSEA, BGammaHat, adjBGammaHat, BMc, 
                      BCFI_0C, BTLI_0C, BNFI_0C, 
                      BCFI_0A, BTLI_0A, BNFI_0A, 
                      chisq, pd)
        fitresults <- matrix(modelfit, nrow = 1)
        colnames(fitresults) <- fitnom_blavaan
        
        return(fitresults)
      }, error = function(e) {
        cat("Error in fit indices calculation:", conditionMessage(e), "\n")
        return(rep(NA, length(fitnom_blavaan)))
      })
    } else {
      cat("Model did not converge\n")
      return(rep(NA, length(fitnom_blavaan)))
    }
  } else {
    cat("Model fitting failed\n")
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
    rep(NA, length(fitnom_stan))
  }
} 

# 2.1. Read from slurm bash file: workload, specific seed from array, core? #######################################

# Get simulation number from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript run_simulation.R <simulation_number>")
}
sim_num <- as.numeric(args[1]) # = 5

# Get current simulation workload
current_sim <- workloads[[paste0("sim_", sim_num)]]
if (is.null(current_sim)) {
  stop(paste("No configuration found for simulation", sim_num))
}
#current_sim <- 4

# 2.2. Test and make sure slurm unique_id and worloads are working ################################################

# print rslurm_id
rslurm_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) 
#rslurm_id<-1
cat('rslurm_id=', rslurm_id, '\n')

# Workload id
cat("Workload_id =", sim_num, '\n')

# current_sim 
cat('current_sim=', current_sim, '\n')

# Workload specific conditions
N_t <- current_sim$combinations$timepoint
N_p <- current_sim$combinations$person_size
cat("Workload N_p:", N_p, "\n", "Workload N_t:", N_t)


# 3.0 Global vars ##############################################################################################

Conditions_P <- c(91,121, 151, 181,211,501,1001,1501,2001)
Conditions_T <- c(1:5, 10, 15)

# Setting all hyperparameters excluding N_t and N_p 
# Save as filetype (one of local, RDS, csv)
save_as = 'csv'
# Size of crossloading (indicating misfit)
Size_crossloading <- c(0,.3,.6)
# What models to specify and run: within time points/ between time points
Type_crossloading <- c("none","tt","tt1")
# Number of data samples to generate and run model fit on them
N_sim_samples <- 1 # because each core does 1?

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

# 3.0.1 Simulation loop ###########################################################################################
for (i_pers in seq_along(Person_size)) {
  for (i_timep in seq_along(Timepoints)) {
    
    n_p <- Person_size[i_pers] #small n_p because local variable (inside function)
    n_t<- Timepoints[i_timep]
    
    # if condition: match pos i of N_pers to pos i of N_timep. This is needed for workloads
    if(i_pers == i_timep){
      print(paste("Matching position:", i_pers, "-> N_pers:", n_p, "N_timep:", n_t))
      
      for(idk in N_sim_samples){
        
        idk <- rslurm_id
        name_local_SIMULATE_Info <- paste0("results_dfg/results_N",as.character(person_size),"_cluster",as.character(i))
        
        #####################################################################
        #set core specific seed at the start of each gendata
        set.seed(24032024+idk+N_pers)
        
        rnorm(5)
        
          
      }
      # test if sim parallelisation works
      
      # 3.1. For each unique_id: new gendata ###########################################################################
      #TODO
      # for each dataset (tt, tt1/none) reinit the seed
      
      # 3.2. For each unique_id: num_sim: 25, 250 respectively
      #TODO
      
      # 3.3. For each unique_id: run sim and make sure it is correctly parallelised 
      #TODO
      
      # 4.1. In each parallelisation run, save file ####################################################################
      #TODO
      # Find correct folder
      
      #with specific seed, N, Nt, lav/blav, tt_type, tt_size,   
      # TODO
      
    } #end condition for correct workload spec. Sligt redundancy, but O(1)
  } #end loop over N_t
} #end loop over N_p



