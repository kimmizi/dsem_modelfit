# run_par_new.R

# Plan: 
# 1. Initialise everything inside the file so that there are no problems with the parallel workers (this should work indep of part 2)
# 2. Each file gets a seed from the slurm bash file
# 3. For each seed, we have a new gendata, new sim
# 4. Somehow saves each run in a separate file

debugging_mode = "LOCAL"  #"LOCAL" 

# 1.1. Packages and file dependencies #######################################################################

# all required packages should be here and initialised only once at the start
library(mvtnorm)
library(matrixcalc)
library(lavaan)
library(blavaan)
library(rstan)
# I think this is so that rstan and blavaan can operate efficiently. 
# How might this interact with the other part of the simulation?
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# 1.2. DSEM models ###########################################################################################

if(debugging_mode != "LOCAL"){
  # ONLY FOR CLUSTER USE
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
  
}else{
  # ONLY FOR LOCAL USE
  current_dir = "/Users/kimzierahn/PycharmProjects/dsem_modelfit/src"
  
}


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
  ),
  # this last one is a testing one
  sim_11 = list(
    total_workload = 20618,
    combinations = data.frame(
      timepoint = c(1, 2),
      person_size = c(31, 91)
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


#### A first function generating data using within-time points cross-loadings ####
# Inputs:
# N: number of persons
# Nt: number of timepoints
# phi0: covariance matrix for latent factors
# mu0: mean vector for latent factors
# ar0: autoregressive coefficients vector
# ly0: factor loadings (for current time)
# ly1: factor loadings (for cross-loadings) ### WHAT IS DIFFERENCE BETWEEN LY0 AND LY1?
# td: residual variance

# Output:
# a list containing
# is_positive_def: Is cov positive-definite? TRUE / FALSE
# y: observed data over time and individuals
# y0: reshaped data frame
# x_cov: covariance matrix of y0
# x_mean: mean vector of y0

gendata01 <- function(N, Nt, phi0, mu0, ar0, ly0, ly1, td){
  
  # random intercept with var = .3, so icc = .3
  # scale covariance matrix
  covu0 <- phi0*.3
  
  # sample latent factors from multivariate normal distribution
  eta2 <- rmvnorm(N, mu0, covu0) # RANDOM INTERCEPT VALUES?
  
  # initialize empty matrices for 2 latent variables (eta) and 6 observed variables (y)
  eta <- array(NA, c(Nt, N, 2))
  y   <- array(NA, c(Nt, N, 6))
  
  # update residual variance for standardized items
  # adjust the fourth element of matrix
  td[4,4] <- 1 - (ly0[4, 2]^2 + ly1[4, 1]^2 + 2*ly0[4, 2]*ly1[4, 1]*phi0[2, 1])
  
  # latent variables
  # time point 1
  eta[1, , ] <- eta2 + rmvnorm(N, mu0, phi0 - covu0)
  
  # observed variables
  y[1, , ] <- eta[1, , ]%*%t(ly0 + ly1) + rmvnorm(N, sigma = td) # FACTOR MODEL, ly1 is a missspecification lambda 
  #cov(eta[1,,])
  
  # loop through remaining time points
  if(Nt > 1){
    for(j in 2:Nt){
      # sample zeta, residuals for latent variables
      zeta <- rmvnorm(N, c(0, 0), (phi0)*(1 - ar0^2) - covu0) # this is a residual with var (1-phi^2) sp tjat var(eta)==1
      
      # update the two latent variables
      # adds the autoregressive effect (ar0) of the previous time point and residual
      eta[j, , 1] <- eta2[, 1] + ar0[1]*(eta[j - 1, , 1] - eta2[, 1]) + zeta[, 1] # eta2[,1] INTERCEPT
      eta[j, , 2] <- eta2[, 2] + ar0[2]*(eta[j - 1, , 2] - eta2[, 2]) + zeta[, 2]
      
      # generate observed data
      # this does NOT include cross-time effects
      y[j, , ] <- eta[j, , ]%*%t(ly0 + ly1) + rmvnorm(N, sigma = td) # ARRAY 3D, IN VERSION 04 WE ADD THE CROSSTIME ADDITIONALLY 
    }
  }
  # NO NEED TO ADJUST
  #cov(eta[Nt,,]) # check points, is cov what I want it to be? 0.3
  #diag(cov(y[Nt,,]))
  #apply(y[Nt,,],2,mean)
  
  # initialize empty matrix with N rows (individuals) and 6*Nt columns (6 observed variables over all time points)
  y0 <- matrix(NA, N, 6*Nt)
  
  for(i in 1:N){
    for(j in 1:Nt){
      y0[i, (j - 1)*6 + 1:6] <- y[j, i, ]  
    }
  }
  
  # convert observed data y into data frame
  y0 <- data.frame(y0)
  
  # add columns names
  cnom <- paste0("y", 1:6, "t", 1)
  if(Nt > 1){
    for(j in 2:Nt){
      cnom <- c(cnom, paste0("y", 1:6, "t", j))
    }
  }
  colnames(y0) <- cnom
  
  # Compute covariance
  x_cov  <- cov(y0)
  
  # Compute mean
  x_mean <- apply(y0, 2, mean)
  
  # Check if cov is positive-definite
  is_positive_def <- is.positive.definite(x_cov) 
  
  # Output depending on whether cov is pos.-def.
  if(is_positive_def==T){
    out <- list(is_positive_def,y,y0,x_cov,x_mean)
    names(out) <- c("is_positive_def", "y", "y0", "x_cov", "x_mean")
  }else{
    out <- list(is_positive_def)
    names(out) <- "is_positive_def"
  }
  out
}


#### Second function generating data using between-time points cross-loadings ####

# Inputs:
# N: number of persons
# Nt: number of timepoints
# phi0: covariance matrix for latent factors
# mu0: mean vector for latent factors
# ar0: autoregressive coefficients vector
# ly0: factor loadings (for current time)
# ly1: factor loadings (for cross-loadings) ### WHAT IS DIFFERENCE BETWEEN LY0 AND LY1?
# td: residual variance

# Output:
# a list containing
# is_positive_def: Is cov positive-definite? TRUE / FALSE
# y: observed data over time and individuals
# y0: reshaped data frame
# x_cov: covariance matrix of y0
# x_mean: mean vector of y0

gendata02 <- function(N, Nt, phi0, mu0, ar0, ly0, ly1, td){
  
  # scale covariance matrix
  covu0 <- phi0*.3
  
  # sample latent factors from multivariate normal distribution
  eta2 <- rmvnorm(N, mu0, covu0)
  
  # initialize empty matrices for 2 latent variables (eta) and 6 observed variables (y)
  eta <- array(NA, c(Nt, N, 2))
  y <- array(NA, c(Nt, N, 6))
  
  ### WHAT DOES THIS DO EXACTLY?
  # update residual variance for standardized items
  tdtt1 <- td
  # adjust the first element of matrix 
  tdtt1[1,1] <- 1 - (ly0[1, 1]^2 + ly1[1, 1]^2 + 2*ly0[1, 1]*ly1[1, 1]*ar0[1])
  
  # latent variables
  # time point 1
  eta[1,,] <- eta2 + rmvnorm(N, mu0, phi0 - covu0)
  #cov(eta[,1,]) # check cov
  # NOTE: the total variance of the latent factors is currently not 1
  
  # observed variables
  y[1,,] <- eta[1,,]%*%t(ly0) + rmvnorm(N, sigma = td)
  #cov(y[1,,])
  #apply(y[1,,],2,mean)
  
  # loop through remaining time points
  if(Nt > 1){
    for(j in 2:Nt){
      # sample zeta, residuals for latent variables
      zeta <- rmvnorm(N, c(0, 0), (phi0)*(1 - ar0^2) - covu0) # this is a residual with var (1 - phi^2) so that var(eta)==1
      
      # update the two latent variables
      # adds the autoregressive effect (ar0) of the previous time point and residual
      eta[j, , 1] <- eta2[,1] + ar0[1]*(eta[j - 1, , 1] - eta2[,1]) + zeta[,1]
      eta[j, , 2] <- eta2[,2] + ar0[2]*(eta[j - 1, , 2] - eta2[,2]) + zeta[,2]
      
      # generate observed data
      # this includes cross-time effects
      y[j, , ] <- eta[j, , ]%*%t(ly0) + eta[j - 1, , ]%*%t(ly1) + rmvnorm(N, sigma = tdtt1)
    }
  }
  #diag(cov(y[Nt,,]))
  #apply(y[Nt,,],2,mean)
  
  
  # initialize empty matrix with N rows (individuals) and 6*Nt columns (6 observed variables over all time points)
  y0 <- matrix(NA, N, 6*Nt)
  
  for(i in 1:N){ 
    for(j in 1:Nt){ 
      y0[i, (j - 1)*6 + 1:6] <- y[j, i, ]  
    }
  }
  
  # convert observed data y into data frame
  y0 <- data.frame(y0)
  
  # add columns names
  cnom <- paste0("y", 1:6, "t", 1)
  if(Nt > 1){
    for(j in 2:Nt){
      cnom <- c(cnom, paste0("y", 1:6, "t", j))
    }
  }
  colnames(y0) <- cnom
  
  # Compute covariance
  x_cov <- cov(y0)
  
  # Compute mean
  x_mean <- apply(y0, 2, mean)
  
  # Check if cov is positive-definite
  is_positive_def <- is.positive.definite(x_cov)
  
  # Output depending on whether cov is pos.-def.
  if(is_positive_def == T){
    out <- list(is_positive_def, y, y0, x_cov, x_mean)
    names(out) <- c("is_positive_def", "y", "y0", "x_cov", "x_mean")
  }else{
    out <- list(is_positive_def)
    names(out) <- "is_positive_def"
  }
  out
}


# Wrapper function for the previous two to initialize simulation specific data set
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
      cat("Lav model with ", N_timep, "did not converge.")
      return(rep(NA, length(fitnom_lavaan)))
    }
  } 
  else { # if there is a try-error: save NAs so that simulation doesnt get broken
    cat("Lav model NAs due to try error for timep:", N_timep)
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
        cat("Blav Error in fit indices calculation:", conditionMessage(e), "\n")
        return(rep(NA, length(fitnom_blavaan)))
      })
    } else {
      cat("Blav Model did not converge\n")
      return(rep(NA, length(fitnom_blavaan)))
    }
  } else {
    cat("Blav Model fitting failed\n")
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

if(debugging_mode != "LOCAL"){
  #ONLY FOR CLUSTER USAGE
  # Get simulation number from command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1) {
    stop("Usage: Rscript run_par_new.R <simulation_number>")
  }
  sim_num <- as.numeric(args[1])
}else{
  # ONLY FOR LOCAL USAGE
  sim_num = 11 #always use 11 as a test condition
}

sim_num = 10

# Get current simulation workload
current_sim <- workloads[[paste0("sim_", sim_num)]]
if (is.null(current_sim)) {
  stop(paste("No configuration found for simulation", sim_num))
}
#current_sim <- 4


# 2.2. Test and make sure slurm unique_id and worloads are working ################################################

if(debugging_mode != "LOCAL"){
  # ONLY FOR CLUSTER USAGE
  # print rslurm_id
  rslurm_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) 
}else{
  #ONLY FOR LOCAL USAGE
  rslurm_id<-1
}


cat('rslurm_id=', rslurm_id, '\n')

# Workload id
cat("Workload_id =", sim_num, '\n')

# Workload specific conditions
# Number of people / time point?
#N_p <- c(91)
#N_p <- c(31, 61, 91, 121, 151, 181, 211, 501, 1001, 1501, 2001)
N_p <- c(31, 61, 91)


# Number of measurement time points
#N_t <- c(2) # N_timep
#N_t <- c(1:5, 10, 15) # N_timep
N_t <- c(1:3) # N_timep


cat("Workload N_p:", N_p, "\n")
cat("Workload N_t:", N_t, "\n")



# 3.0 Global vars ##############################################################################################

# Setting all hyperparameters excluding N_t and N_p 
# Save as filetype (one of local, RDS, csv)
save_as = 'csv'
# Size of crossloading (indicating misfit)
Size_crossloading <- c(0, .3, .6)
# What models to specify and run: within time points/ between time points
Type_crossloading <- c("none", "tt", "tt1")
# Number of data samples to generate and run model fit on them
N_sim_samples <- 10 # because each core does 1?

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
for (i_pers in seq_along(N_p)) {
  for (i_timep in seq_along(N_t)) {
    
    n_p <- N_p[i_pers] #small n_p because local variable (inside function)
    n_t<- N_t[i_timep]
    cat("Person size: ", n_p, "\n")
    cat("Time points: ", n_t, "\n")
    
    # if condition: match pos i of N_pers to pos i of N_timep. This is needed for workloads
    #if(i_pers == i_timep){
    #  print(paste("Matching position:", i_pers, "-> N_p:", n_p, "N_t:", n_t))
      
      
      # Check theoretical specification requirement
      if (n_t * 6 >= n_p) next
      for (Type_misfit in Type_crossloading) {
        # Picking conditions. Since we are not interested in power, we just want
        # to fit 1 model with no misfit & type none, and then foreach misfit type
        # 1 with misfit size = 0.3, one 0.6.
        #cat(Type_misfit, "\n")
          
        #set core specific seed right before gendata
        #set.seed(24032024+idk+n_p)
        set.seed(131212023 + n_p + n_t * 1000 + as.numeric(as.factor(Type_misfit)) * 10000)
          
        for (Size_misfit in Size_crossloading) {
          # if type tt or tt1, dont run size=0, else run
          if(Size_misfit == 0 && Type_misfit != "none"){
            #cat('skipped over misfit type:', Type_misfit, " and size:", Size_misfit, "\n")
          }else if(Size_misfit != 0 && Type_misfit == "none"){
            #cat('skipped over misfit type:', Type_misfit, " and size:", Size_misfit, "\n")
          }else{ 
            
            #cat(Size_misfit, "\n")
            
            # Initialize factor loadings for the current condition
            ly1 <- initialize_ly1(Type_misfit, Size_misfit)
            
            # Initialize empty matrix to store fit indices
            fitm_lavaan <- matrix(NA, N_sim_samples, length(fitnom_lavaan))
            #fitm_blavaan <- matrix(NA, N_sim_samples, length(fitnom_blavaan))
            #fitm_stan <- matrix(NA, 1, length(fitnom_stan))
            
              
            #print(N_sim_samples)
            for(idk in 1:N_sim_samples){
                
              #idk <- rslurm_id
              
              temp_dat <- generate_temp_data(Type_misfit, n_p, n_t, phi0, mu0, ar0, ly0, ly1, td)
              ydat <- temp_dat$y # use array version of data for stan
              ydat2 <- temp_dat$y0 # use the transformed version in the getfit function
              
                
                # Check if data is positive definite and fit models
                if (temp_dat[["is_positive_def"]]) {
                  y0 <- temp_dat[["y0"]]
                
                  # Lavaan
                  fitm_lavaan[idk, ] <- fit_model_lav(ydat2, n_t) # Fit model and store results
                
                  # Blavaan
                  #fitm_blavaan[1, ] <- fit_model_blav(ydat2, n_t) # Fit model and store results
                
                  # Stan
                  #fitm_stan[1, ] <- as.numeric(fit_model_stan(ydat, ydat2, n_t, n_p)) # Fit model and store results
              } # end if pos definite, fit models
            } #end loop over number of simulations per slurm worker 
            
            # for each misfit type and misfit size
              
              # 4.1. After each parallelisation run, save file ###################################################################
              
            # Colnames?
            colnames(fitm_lavaan) <- fitnom_lavaan
            #colnames(fitm_blavaan) <- fitnom_blavaan
            # #colnames(fitm_stan) <- fitnom_stan
            # Move up one directory level to dsem_modelfit
            parent_dir <- dirname(current_dir)
            # Define the experiment name
            Exp_day <- format(Sys.time(), "%Y-%m-%d") # Format: YYYY-MM-DD_HH-MM-SS
            # Specify the save directory in dsem_modelfit/exp
            save_dir <- file.path(parent_dir, "exp", Exp_day)
            # Create the day directory if it does not exist
            # TODO
            if (!dir.exists(save_dir)) {
              dir.create(save_dir, recursive = TRUE)
            }
            
            # test
            #exp_name <- paste("test", n_p, n_t, idk, "lav_blav", ".csv", sep='_')
            # csv_path_test <- file.path(save_dir, exp_name)
            # write.csv(random_numbers, file=csv_path_test, row.names = FALSE)
            
            # actual
            exp_name_lav <- paste("dsem", n_p, n_t, idk, Type_misfit, Size_misfit, "lav", ".csv", sep='_')
            #exp_name_blav <- paste("dsem", n_p, n_t, idk, Type_misfit, Size_misfit, "blav", ".csv", sep='_')
            
            csv_path_lavaan <- file.path(save_dir, exp_name_lav)
            write.csv(fitm_lavaan, file = csv_path_lavaan, row.names = FALSE)
            #csv_path_blavaan <- file.path(save_dir, exp_name_blav)
            #write.csv(fitm_blavaan, file = csv_path_blavaan, row.names = FALSE)
      
          } # end if selecting 5 model conditions from 9 possible
        } #end loop over misfit size
      } #end loop over misfit type
    #} #end condition for correct workload spec. Sligt redundancy, but O(1)
  } #end loop over N_t
} #end loop over N_p

cat("finished worker number ", rslurm_id, "from Workload", sim_num)


