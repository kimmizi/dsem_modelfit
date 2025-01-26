#########################################################
#### Nullmodels ####
#########################################################



#########################################################
#### configural invariance model ####
# has the same pattern of fixed and free loadings in the factor loading matrix 
# in each time point but no cross-time invariance constraints in any matrices and
# is thus the least constrained model
#########################################################

configural_invariance_model <- function(timepoints) {
  model_string <- ''
  
  for (t in 1:timepoints) {
    
    # fixed latent structure
    model_string <- paste(model_string, sprintf('
      eta1t%d =~ ly1*y1t%d + ly2*y2t%d + ly3*y3t%d
      eta2t%d =~ ly4*y4t%d + ly5*y5t%d + ly6*y6t%d
    ', t, t, t, t, t, t, t, t), sep="\n")
  }
  
  return(model_string)
}

configural_invariance_model(2)



#########################################################
#### weak factorial invariance model ####
# is identical to the configural invariance model, 
# but invokes cross-time invariance constraints on factor loadings in the factor loading matrix 
#########################################################

#weak_invariance_model <- function(timepoints) {
#  model_string <- configural_invariance_model(timepoints)
#
#  return(model_string)
#}
# not applicable here



#########################################################
#### strong factorial invariance model ####
# adds to the weak factorial invariance model cross-time invariance constraints
# on intercept terms
#########################################################

strong_invariance_model <- function(timepoints) {
  model_string <- configural_invariance_model(timepoints)
  
  # Intercepts Invariance across time points
  for (t in 1:(timepoints-1)) {
    model_string <- paste(model_string, sprintf('
      y1t%d ~ 1*y1t%d
      y2t%d ~ 1*y2t%d
      y3t%d ~ 1*y3t%d
      y4t%d ~ 1*y4t%d
      y5t%d ~ 1*y5t%d
      y6t%d ~ 1*y6t%d
    ', t, t+1, t, t+1, t, t+1, t, t+1, t, t+1, t, t+1), sep="\n")
  }
  
  return(model_string)
}



#########################################################
#### strict factorial invariance model ####
# adds to the strong factorial invariance
# model cross-time constraints on unique factor variances
#########################################################

strict_invariance_model <- function(timepoints) {
  model_string <- strong_invariance_model(timepoints)
  
  # Residual Variance Constraints across time points (Strict Invariance)
  for (t in 1:(timepoints - 1)) {
    model_string <- paste(model_string, sprintf('
      y1t%d ~~ y1t%d
      y2t%d ~~ y2t%d
      y3t%d ~~ y3t%d
      y4t%d ~~ y4t%d
      y5t%d ~~ y5t%d
      y6t%d ~~ y6t%d
    ', t, t+1, t, t+1, t, t+1, t, t+1, t, t+1, t, t+1), sep="\n")
  }
  
  return(model_string)
}




#########################################################
#### null model 0C ####
# variables are uncorrelated across time points
# free intercepts and residual variances
#########################################################

null_model_0C <- function(timepoints) {
  model <- ''
  
  for (t in 1:timepoints) {
    model <- paste0(model, '
      y1t', t, ' ~~ y1t', t, '
      y2t', t, ' ~~ y2t', t, '
      y3t', t, ' ~~ y3t', t, '
      y4t', t, ' ~~ y4t', t, '
      y5t', t, ' ~~ y5t', t, '
      y6t', t, ' ~~ y6t', t, '
      ')
  }
  
  return(model)
}



#########################################################
#### null model 0A ####
# variables are uncorrelated across time points
# invariant intercepts and residual variances
#########################################################

null_model_0A <- function(timepoints) {
  model <- ""
  
  # Define variances and intercepts for all timepoints
  for (t in 1:timepoints) {
    model <- paste0(model, '
      y1t', t, ' ~~ y1t', t, '
      y2t', t, ' ~~ y2t', t, '
      y3t', t, ' ~~ y3t', t, '
      y4t', t, ' ~~ y4t', t, '
      y5t', t, ' ~~ y5t', t, '
      y6t', t, ' ~~ y6t', t, '\n 
      '
                    )
  }
  
  # Add covariance constraints only if timepoints > 1
  if (timepoints > 1) {
    for (t in 1:timepoints) {
      model <- paste0(model, '
      y1t', t, ' ~ 1*1', '
      y2t', t, ' ~ 1*1', '
      y3t', t, ' ~ 1*1', '
      y4t', t, ' ~ 1*1', '
      y5t', t, ' ~ 1*1', '
      y6t', t, ' ~ 1*1', '\n
      '
      )
    }
    
    for (v in 1:6) { # Loop over variables y1 to y6
      var_name <- paste0("y", v, "t") # Append "t" to indicate timepoint
      for (t in 1:(timepoints - 1)) { # Loop over time points
        next_time <- t + 1
        model <- paste0(model, 
                        sprintf("%s%d ~~ 0.5*%s%d\n", var_name, t, var_name, next_time))
      }
    }
  }
  
  # Add intercept constraints only if timepoints > 1
  if (timepoints > 1) {
    for (t in 1:(timepoints - 1)) {
      model <- paste(model, sprintf('
        y1t%d ~ 1*y1t%d
        y2t%d ~ 1*y2t%d
        y3t%d ~ 1*y3t%d
        y4t%d ~ 1*y4t%d
        y5t%d ~ 1*y5t%d
        y6t%d ~ 1*y6t%d
      ', t, t + 1, t, t + 1, t, t + 1, t, t + 1, t, t + 1, t, t + 1), sep = "\n")
    }
  }
  
  return(model)
}

nullmodel_0A <- try(bsem(null_model_0A(2), data = y0, 
                         n.chains = 4, burnin = 1000, sample = 1000), silent = F)



#########################################################
#### null model with correlated residuals ####
#########################################################

null_model_corr_res <- function(timepoints) {
  model <- ''
  
  # Uncorrelated residuals within each timepoint
  for (t in 1:timepoints) {
    model <- paste0(model, '
      y1t', t, ' ~~ y1t', t, '
      y2t', t, ' ~~ y2t', t, '
      y3t', t, ' ~~ y3t', t, '
      y4t', t, ' ~~ y4t', t, '
      y5t', t, ' ~~ y5t', t, '
      y6t', t, ' ~~ y6t', t, '
      ')
  }
  
  # Correlated residuals across timepoints for corresponding variables
  for (t in 1:(timepoints - 1)) {
    model <- paste0(model, '
      y1t', t, ' ~~ y1t', t + 1, '
      y2t', t, ' ~~ y2t', t + 1, '
      y3t', t, ' ~~ y3t', t + 1, '
      y4t', t, ' ~~ y4t', t + 1, '
      y5t', t, ' ~~ y5t', t + 1, '
      y6t', t, ' ~~ y6t', t + 1, '
      ')
  }
  
  return(model)
}


