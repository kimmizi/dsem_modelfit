# This file runs the following
# - set hyperparams for 1 run
# - data generation function
# - fit models
# - get fit indices
# - return results for later analysis



# 0. importing relevant libraries and functions ################################

# 
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
#source(file.path(current_dir,  "gen_data_version04.R"))

# working dir?
#setwd()



# 1. Setting hyperparameters ###################################################

# Simulation params

# Number of cores for parallelisation
cores <- 1

# Number of people / time point?
person_size_SIMULATE <- c(91)
#person_size_SIMULATE <- c(91, 121, 151, 181, 211, 501, 1001, 1501, 2001, 2501, 61, 31)

# Number of measurement time points
time_point_SIMULATE <- 15 # Nt
#time_point_SIMULATE <- c(1:5, 10, 15) # Nt

# Size of crossloading (indicating misfit)
model_TRUE_MISS_SIMULATE <- c(.3,.6)

# What models to specify and run: within time points/ between time points
type_TRUE_MISS_SIMULATE <- c("none","tt","tt1")

# Number of data samples to generate and run model fit on them
run_Samples_SIMULATE <- 3

# DSEM model params
phi0 <- diag(2)*.7+.3 # cov(eta)
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

#summary(res1,standardized=T)



# Parallelise cores ############################################################
# TODO: for core in cores
# temporary core=1 for 1 run
core <- 1

# Generate data: right now only taking y0 from get_data_version03.R
#dat <- y0

# Function to fit a single model, capture fit specific errors for later
# Input: data y0
# Output: model res1
fit_model <- function(y0){ 
  
  res1 <- try(sem(dsem[[time_point]], data = y0, std.lv = TRUE), silent = TRUE)
  
  if(!inherits(res1, "try-error")){ # if there is no try-error
    
    if(res1@optim$converged == T){ # if model converged: save fit indices
      fitmeasures(res1)  
    }
    else
      { # if model did not converge: save NAs
      rep(NA, 46)
      }
  }
  else
    { # if there is a try-error: save NAs
    rep(NA, 46)
    }
}

# what does this line do?
#diag((ly0)%*%phi0%*%t(ly0)+td)



# Huge loop iterating over all simulation combinations

for (person_size in person_size_SIMULATE) { # looping over persons
  
  for (time_point in time_point_SIMULATE) { # looping over time points
    
    if(time_point*6 < person_size){ # checking if N_t+p < Nt? Specification requirement
      
      for(type_MISS in type_TRUE_MISS_SIMULATE){ # check each misfit type
      
      # set core specific seed
      set.seed(131212023 + core + person_size + time_point*1000 + as.numeric(as.factor(type_MISS))*10000)

        
      #####################################################################
      if(type_MISS == "none"){ # core <- 1
        
        # there are no missspecifications
        model_TRUE_MISS <- 0
        
        # save in current directory
        name_local_SIMULATE_Info <- paste(current_dir, as.character(person_size), as.character(time_point), as.character(type_MISS), 
                                          as.character(model_TRUE_MISS), core, "_version03_rand", sep = "_")
        
        
        #######################################
        # empty matrix for results
        #######################################
        #run_Samples_SIMULATE <- SAMPLING <- 1
        
        # Initialize empty matrix that stores fit indices
        fitm1 <- matrix(NA, run_Samples_SIMULATE, 46)
        N <- person_size
        Nt <- time_point
        ly1  <- matrix(c(0, 0, 0, model_TRUE_MISS, 0, 0,
                         0, 0, 0, 0, 0, 0), 6, 2, byrow = F) # factor loadings
        
        
        #######################################
        
        
        # Start looping over number of samplings
        for (SAMPLING in 1:run_Samples_SIMULATE){ # SAMPLING <- 1
          
          # Generate data 
          dat1 <- gendata01(N, Nt, phi0, mu0, ar0, ly0, ly1, td)
          
          # Check if data is positive definite
          if(dat1[["is_positive_def"]] == T){
            
            y0 <- dat1[["y0"]]
            
            # Fit model and store fit indices of this model in previously defined matrix, the row indicating the number of sampling
            fitm1[SAMPLING,] <- fit_model(y0)
            
          }
        } # end looping sampling
        
        # Adjust column names of the fit indices matrix
        colnames(fitm1) <- fitnom
        
        # Save samples run
        saveRDS(fitm1, file = paste0(name_local_SIMULATE_Info,".RDS")) 
        
        
        
        
        #####################################################################
        #####################################################################
      }
        else
          { # type_MISS <- "tt1": between time points
            
        for(model_TRUE_MISS in model_TRUE_MISS_SIMULATE){ # check each misfit type # model_TRUE_MISS <- .3
          
          # save in current directory
          name_local_SIMULATE_Info <- paste(current_dir, as.character(person_size), as.character(time_point), as.character(type_MISS), 
                                            as.character(model_TRUE_MISS), core, "_version03_rand", sep = "_")
          
          
          #######################################
          # empty matrix for results
          #######################################
          
          # Initialize empty matrix that stores fit indices
          fitm1 <- matrix(NA, run_Samples_SIMULATE, 46)
          N <- person_size
          Nt <- time_point
          
          if(type_MISS == "tt"){
            ly1  <- matrix(c(0, 0, 0, model_TRUE_MISS, 0, 0,
                             0, 0, 0, 0, 0, 0), 6, 2, byrow = F) # factor loadings
          }
          else
            { # model_TRUE_MISS <- .3
            ly1  <- matrix(c(model_TRUE_MISS, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0), 6, 2, byrow = F) # factor loadings
            }
          
          
          #######################################
          
          
          # Start looping over number of samplings
          for (SAMPLING in 1:run_Samples_SIMULATE){ # SAMPLING <-1
            
            # Generate data 
            if(type_MISS == "tt"){
              dat1 <- gendata01(N, Nt, phi0, mu0, ar0, ly0, ly1, td)
            }
            else
              {
              dat1 <- gendata02(N, Nt, phi0, mu0, ar0, ly0, ly1, td)
            }
            
            
            # Check if data is positive definite
            if(dat1[["is_positive_def"]] == T){
              
              y0 <- dat1[["y0"]]
              
              # Fit model and store fit indices of this model in previously defined matrix, the row indicating the number of sampling
              fitm1[SAMPLING,] <- fit_model(y0)
              
            } 
          } # end looping sampling
          
          # Adjust column names of the fit indices matrix
          colnames(fitm1) <- fitnom
          
          # Save samples run
          saveRDS(fitm1, file = paste0(name_local_SIMULATE_Info,".RDS")) 
          
        } # end looping over misfit types
            
      } # end else
      
      } # end looping over misfit types
      
    } # end checking specification requirement
    
  } #end looping over time points
  
} #end looping over persons


