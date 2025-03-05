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
# intercepts and residual variances constrained to be equal over time
#########################################################

null_model_0A <- function(timepoints) {
  model <- ''
  
  # Define intercepts (constrained to be equal across time)
  model <- paste0(model, '
    # Intercepts
    y1t1 ~ int_y1 * 1
    y2t1 ~ int_y2 * 1
    y3t1 ~ int_y3 * 1
    y4t1 ~ int_y4 * 1
    y5t1 ~ int_y5 * 1
    y6t1 ~ int_y6 * 1
  ')
  
  # Constrain intercepts to be equal across time points
  for (t in 2:timepoints) {
    model <- paste0(model, '
      y1t', t, ' ~ int_y1 * 1
      y2t', t, ' ~ int_y2 * 1
      y3t', t, ' ~ int_y3 * 1
      y4t', t, ' ~ int_y4 * 1
      y5t', t, ' ~ int_y5 * 1
      y6t', t, ' ~ int_y6 * 1
    ')
  }
  
  # Define residual variances (constrained to be equal across time)
  model <- paste0(model, '
    # Residual variances
    y1t1 ~~ res_y1 * y1t1
    y2t1 ~~ res_y2 * y2t1
    y3t1 ~~ res_y3 * y3t1
    y4t1 ~~ res_y4 * y4t1
    y5t1 ~~ res_y5 * y5t1
    y6t1 ~~ res_y6 * y6t1
  ')
  
  # Constrain residual variances to be equal across time points
  for (t in 2:timepoints) {
    model <- paste0(model, '
      y1t', t, ' ~~ res_y1 * y1t', t, '
      y2t', t, ' ~~ res_y2 * y2t', t, '
      y3t', t, ' ~~ res_y3 * y3t', t, '
      y4t', t, ' ~~ res_y4 * y4t', t, '
      y5t', t, ' ~~ res_y5 * y5t', t, '
      y6t', t, ' ~~ res_y6 * y6t', t, '
    ')
  }
  
  return(model)
}


nullmodel_0A <- try(bsem(null_model_0A(Timepoints), ydat2, 
                         n.chains = 4, burnin = 1000, sample = 1000), silent = F)
summary(nullmodel_0A)



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


