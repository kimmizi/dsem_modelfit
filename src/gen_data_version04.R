# generate data

gendata02 <- function(N,Nt,phi0,mu0,ar0,ly0,ly1,td){
  
  #SAMPLING <- 1
  covu0 <- phi0*.3
  eta2 <- rmvnorm(N,mu0,covu0)
  
  # empty matrices for lvs and ovs
  eta <- array(NA,c(Nt,N,2))
  #matrix[N,6] y[Nt]; declared
  y   <- array(NA,c(Nt,N,6))
  
  # update td for standardized items
  tdtt1 <- td
  tdtt1[1,1] <- 1-(ly0[1,1]^2+ly1[1,1]^2+2*ly0[1,1]*ly1[1,1]*ar0[1])
  
  # latent variables
  # time point 1
  eta[1,,] <- eta2 + rmvnorm(N,mu0,phi0-covu0)
  #cov(eta[,1,]) # check cov
  # NOTE: the total variance of the latent factors is currently not 1
  y[1,,] <- eta[1,,]%*%t(ly0)+rmvnorm(N,sigma = td)
  #cov(y[1,,])
  #apply(y[1,,],2,mean)
  if(Nt>1){
    for(j in 2:Nt){#j<-2
      zeta <- rmvnorm(N,c(0,0),(phi0)*(1-ar0^2)-covu0) # this is a residual with var (1-phi^2) sp tjat var(eta)==1
      eta[j,,1] <- eta2[,1] + ar0[1]*(eta[j-1,,1]-eta2[,1]) + zeta[,1]
      eta[j,,2] <- eta2[,2] + ar0[2]*(eta[j-1,,2]-eta2[,2]) + zeta[,2]
      y[j,,] <- eta[j,,]%*%t(ly0)+eta[j-1,,]%*%t(ly1)+rmvnorm(N,sigma = tdtt1)
    }
  }
  #diag(cov(y[Nt,,]))
  #apply(y[Nt,,],2,mean)
  
  y0 <- matrix(NA,N,6*Nt)
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
  x_cov  <- cov(y0)
  
  #dim(x_cov)
  #round(eigen(x_cov)$values,3)
  x_mean <- apply(y0,2,mean)
  is_positive_def <- is.positive.definite(x_cov)
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
