# generate data

gendata01 <- function(N, Nt, phi0, mu0, ar0, ly0, ly1, td){
  #N<-10000
  #SAMPLING <- 1
  # random intercept with var=.3, so icc=.3
  covu0 <- phi0*.3
  eta2 <- rmvnorm(N, mu0, covu0) # RANDOM INTERCEPT VALUES
  
  # empty matrices for lvs and ovs
  eta <- array(NA,c(Nt,N,2))
  #matrix[N,6] y[Nt]; declared
  y   <- array(NA,c(Nt,N,6))
  
  # update td for standardized items
  # UPDATED RESIDUAL VARIANCE
  td[4,4] <- 1-(ly0[4,2]^2+ly1[4,1]^2+2*ly0[4,2]*ly1[4,1]*phi0[2,1])
  
  # time point 1
  eta[1,,] <- eta2 + rmvnorm(N,mu0,phi0-covu0)
  y[1,,] <- eta[1,,]%*%t(ly0+ly1)+rmvnorm(N,sigma = td) # FACTOR MODEL, ly1 is a missspecification lambda 
  #cov(eta[1,,])
  #remaining time points
  if(Nt>1){
    for(j in 2:Nt){#j<-2
      zeta <- rmvnorm(N,c(0,0),(phi0)*(1-ar0^2)-covu0) # this is a residual with var (1-phi^2) sp tjat var(eta)==1
      eta[j,,1] <- eta2[,1] + ar0[1]*(eta[j-1,,1]-eta2[,1]) + zeta[,1] # eta2[,1] INTERCEPT
      eta[j,,2] <- eta2[,2] + ar0[2]*(eta[j-1,,2]-eta2[,2]) + zeta[,2]
      y[j,,] <- eta[j,,]%*%t(ly0+ly1)+rmvnorm(N,sigma = td) # ARRAY 3D, IN VERSION 04 WE ADD THE CROSSTIME ADDITIONALLY 
    }
  }
  # NO NEED TO ADJUST
  #cov(eta[Nt,,]) # check points, is cov what I want it to be? 0.3
  #diag(cov(y[Nt,,]))
  #apply(y[Nt,,],2,mean)
  
  y0 <- matrix(NA,N,6*Nt) # N ROWS, N*6 COLUMNS (item 1 - 6 for timepoint 1, then item 1-6 for timepoint 2, ...)
  for(i in 1:N){#i<-1
    for(j in 1:Nt){#j<-1
      y0[i,(j-1)*6+1:6] <- y[j,i,]  
    }
  }
  y0 <- data.frame(y0)
  cnom <- paste0("y",1:6,"t",1)
  if(Nt>1){
    for(j in 2:Nt){
      cnom <- c(cnom,paste0("y",1:6,"t",j))
    }
  }
  colnames(y0)<-cnom
  x_cov  <- cov(y0) # CALC COV MATRIX
  
  #dim(x_cov)
  #round(eigen(x_cov)$values,3)
  x_mean <- apply(y0,2,mean)
  is_positive_def <- is.positive.definite(x_cov) # CHECK IF COV IS POSITIVE DEFINITE, SHOULD WE SINCE WE PUT IF N > N+... IN DATA GEN PROCESS
  #is_positive_def
  if(is_positive_def==T){
    out <- list(is_positive_def,y,y0,x_cov,x_mean)
    names(out) <- c("is_positive_def","y","y0","x_cov","x_mean")
  }else{
    out <- list(is_positive_def)
    names(out) <- "is_positive_def"
  }
  out
}
