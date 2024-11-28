#### A function generating data using within-time points cross-loadings ####

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
