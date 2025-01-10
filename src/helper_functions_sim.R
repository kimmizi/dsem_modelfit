library(blavaan)
library(lavaan)
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
      return(rep(NA, length(fitnom_lavaan)))
    }
  } 
  else { # if there is a try-error: save NAs so that simulation doesnt get broken
    return(rep(NA, length(fitnom_lavaan)))
  }
}



#### Bayesian models ####

### BLAVAAN ###
fit_model_blav <- function(y0, N_timep){ 
  
  print("nullmodel_0A")
  nullmodel_0A <- try(bsem(null_model_0A(N_timep), data = y0, 
                           n.chains = 4, burnin = 1000, sample = 1000), silent = TRUE)
  
  print("res1")
  res1 <- try(bsem(dsem[[N_timep]], data = y0, std.lv = TRUE,
                   n.chains = 4, burnin = 1000, sample = 1000), silent = TRUE)
  
  print("nullmodel_0C")
  nullmodel_0C <- try(bsem(null_model_0C(N_timep), data = y0, 
                           n.chains = 4, burnin = 1000, sample = 1000), silent = TRUE)
  
  # if there is no try-error
  if(!inherits(res1, "try-error")){ 
    # if model converged: save fit indices
    if(blavInspect(res1, "converged") == TRUE){ 
      
      # try out first nullmodel
      fit_indices <- blavFitIndices(res1, baseline.model = nullmodel_0C) 
      BCFI_0C <- mean(fit_indices@indices$BCFI)
      BTLI_0C <- mean(fit_indices@indices$BTLI)
      BNFI_0C <- mean(fit_indices@indices$BNFI)
      
      # try out different nullmodel
      fit_indices <- blavFitIndices(res1, baseline.model = nullmodel_0A) 
      BCFI_0A <- mean(fit_indices@indices$BCFI)
      BTLI_0A <- mean(fit_indices@indices$BTLI)
      BNFI_0A <- mean(fit_indices@indices$BNFI)
      
      # remaining fit indices
      BRMSEA <- mean(fit_indices@indices$BRMSEA)
      BGammaHat <- mean(fit_indices@indices$BGammaHat)
      adjBGammaHat <- mean(fit_indices@indices$adjBGammaHat)
      BMc <- mean(fit_indices@indices$BMc)
      
      npar <- blavInspect(res1, "npar")
      marg_loglik <- blavInspect(res1, "test")[[1]]$stat
      ppp <- blavInspect(res1, "test")[[2]]$stat
      chisq <- mean(fit_indices@details$chisq)
      pd <- fit_indices@details$pD
      
      # store
      modelfit <- c(npar, ppp, marg_loglik, 
                    BRMSEA, BGammaHat, adjBGammaHat, BMc, 
                    BCFI_0C, BTLI_0C, BNFI_0C, 
                    BCFI_0A, BTLI_0A, BNFI_0A, 
                    chisq, pd)
      fitresults <- matrix(modelfit, nrow = 1)
      colnames(fitresults) <- fitnom_blavaan
      
      return(fitresults)
    }
    else { # if model did not converge: save NAs so that simulation doesnt get broken
      return(rep(NA, length(fitnom_blavaan)))
    }
  } 
  else { # if there is a try-error: save NAs so that simulation doesnt get broken
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

