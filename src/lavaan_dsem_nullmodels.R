#########################################################
#### Nullmodels ####
#########################################################



#########################################################
#### null model 0C ####
# variables are uncorrelated across time points
# free intercepts and residual variances
#########################################################

null_model_0C <- function(timepoints) {
  model <- ''
  
  # Define variances
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
  
  # Define variances
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
  
  # Define intercepts
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
  
  # Define residual variances
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



# test if it worked:
# nullmodel_0A <- try(bsem(null_model_0A(Timepoints), ydat2, 
#                          n.chains = 4, burnin = 1000, sample = 1000), silent = F)
# summary(nullmodel_0A)




