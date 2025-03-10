
getloglik <- function(dat,xmean,xcov){
  #dat <- ydat1
  #xmean <- x_mean
  #xcov <- x_cov
  N <- nrow(dat)
  p <- ncol(dat)
  
  wo.omit <- no.omit <- apply(is.na(dat),1,sum)
  wo.omit[no.omit>0] <- 1
  wo.omit[no.omit==p] <- 2
  #dat <- dat[is.na(no.omit)==F,]
  #no.omit <- na.omit(no.omit)
  
  
  loglik <- rep(NA,N)
  loglik[wo.omit==0] <- dMvn(dat[wo.omit==0,],xmean,xcov,log=T)
  
  which.omit <- which(wo.omit==1)
  
  if(length(which.omit)>1){
    for(i in which.omit){#i<-6
      yi <- dat[i,]
      #yi[1]<-1
      na.yi <- which(is.na(yi))
      loglik[i] <- dMvn(yi[-na.yi],xmean[-na.yi],xcov[-na.yi,-na.yi],log=T)
    }}
  
  loglik
}



dMvn <- function(X,mu,Sigma,log=T) {
  k <- ncol(X)
  rooti <- backsolve(chol(Sigma),diag(k))
  quads <- colSums((crossprod(rooti,(t(X)-mu)))^2)
  return((-(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads))
}


getchisq.mcmc <- function(loglik,logliksat){
  2*(logliksat-loglik)
}


getmcmc <- function(model,jags=T,array=F){
  #model <- fit1
  if(array==F){
    if(jags==T){
      samps <- as.mcmc(model)
      nchains <- length(samps)
      samps2 <- c()
      for(j in 1:nchains){
        samps2 <- rbind(samps2,samps[[j]])  
      }
    }else{
      samps2 <- as.data.frame(model)
    }
  }else{
    samps <- as.mcmc(model)
    nchains <- length(samps)
    niter <- dim(samps[[1]])[1]
    npar <- dim(samps[[1]])[2]
    
    samps2 <- array(NA,c(niter,nchains,npar))
    for(j in 1:nchains){
      samps2[,j,] <- samps[[j]]
    }
    
    dimnames(samps2)[[3]] <- colnames(samps[[j]])
  }
  samps2
}

getmax <- function(x){
  y <- x
  for(i in 1:length(x)){
    y[i] <- max(0,x[i])
  }
  y
}

################################################################
# this includes the whole function in short
################################################################
getfit <- function(fit1,y,y.2){
  samps <- getmcmc(fit1,jags=F) # function extracts samples from jags or stan
  nsamps <- dim(samps)[1]
  #y <- ydat
  N <- dim(y)[2]
  Nt <- dim(y)[1]
  nvar <- dim(y)[3]
  
  # identification constraint
  id0 <- (dim(y)[1]*dim(y)[3]<(dim(y)[2]+1))
  
  #######################################################################
  # parameter extraction
  #######################################################################
  ly1 <- samps[,paste0("ly[",1:4,"]")]
  b1 <- samps[,c(paste0("beta1[",1:2,",1]"),
                 paste0("beta1[",1:2,",2]"))]
  sigma.eps <- samps[,paste0("sigmaeps[",1:6,"]")]
  sigma.eta1 <- samps[,c(paste0("sigmaeta1[",1:2,",1]"),
                         paste0("sigmaeta1[",1:2,",2]"))]
  sigma.eta2 <- samps[,c(paste0("sigmaeta2[",1:2,",1]"),
                         paste0("sigmaeta2[",1:2,",2]"))] ### CHANGED: GOOD
  
  #######################################################################
  # matrices for loglik and model implied covariance matrices
  #######################################################################
  
  #lym.2 <- array(0,c(nvar,2,nsamps)) #HB: added for comparison at the end
  #epsm.2 <- array(0,c(nvar,nvar,nsamps)) #HB: added for comparison at the end
  lym.1 <- array(0,c(nvar*Nt,2*Nt,nsamps))
  epsm.1 <- array(0,c(nvar*Nt,nvar*Nt,nsamps))
  zetam.1 <- etam.1 <-  array(0,c(2*Nt,2*Nt,nsamps))
  
  #etam.2 <- array(0,c(Nt,2,2,nsamps)) # HB: added etam.2.stan for comparison (can be deleted after testing)
  D.1  <- array(0,c(2*Nt,2*Nt,nsamps)) # HB: added D.1.stan for comparison (can be deleted after testing)
  sigmay.1 <- array(NA,c(nvar*Nt,nvar*Nt,nsamps)) # HB: this is overall
  sigmay.2 <- array(NA,c(Nt,nvar,nvar,nsamps)) # HB added: this is time-specific
  m1 <- matrix(1,2,2)
  
  for(j in 1:nsamps){#j<-1
    for(m in 1:Nt){#j<-1
      epsm.1[1:6+(m-1)*nvar,1:6+(m-1)*nvar,j] <- diag(6)*unlist(sigma.eps[j,])^2
    }
  
    for(m in 1:Nt){
      lym.1[c(1,4)+(m-1)*6,(m-1)*2+1:2,j] <- diag(2)
      for(k in 2:3){lym.1[k+(m-1)*6,(m-1)*2+1,j] <- unlist(ly1[j,k-1])}
      for(k in 5:6){lym.1[k+(m-1)*6,(m-1)*2+2,j] <- unlist(ly1[j,k-2])}
    }
    
    for(m in 1:(Nt-1)){#m<-1
      D.1[(m-1)*2+2+1:2,(m-1)*2+1:2,j] = matrix(0, nrow = 2, ncol = 2)-m1*unlist(b1[j,])##diag(2)-m1*unlist(b1[j,])#E)# #### CHANGED. WHY it was commented?
    }
    diag(D.1[,,j]) <- rep(1,2*Nt)
    
    for(m in 1:Nt){#m<-1
      zetam.1[(m-1)*2+1:2,(m-1)*2+1:2,j] = m1*unlist(sigma.eta1[j,])#CHANGED
    }
    
    invD.op <- solve(D.1[,,j])# + t(solve(D.1[,,j])))/2 #HB: CHANGED
    
    etam.1[,,j] <- invD.op %*% zetam.1[,,j] %*% t(invD.op)  #CHANGED
    
    for(m in 1:Nt){### ADDED loop
      etam.1[(m-1)*2+1:2,(m-1)*2+1:2,j] = etam.1[(m-1)*2+1:2,(m-1)*2+1:2,j] + m1*unlist(sigma.eta2[j,])#CHANGED
    }
    
    # HB: final matrices
    sigmay.1[,,j] = lym.1[,,j] %*% etam.1[,,j] %*% t(lym.1[,,j])+ epsm.1[,,j]
    sigmay.1[,,j] = ((sigmay.1[,,j]+ t(sigmay.1[,,j]))/2) #HB: CHANGED
    
    #HB: added for myself
    for(m in 1:Nt){#m<-1
      sigmay.2[m,,,j] <- sigmay.1[1:nvar+(m-1)*nvar,1:nvar+(m-1)*nvar,j]
    }
    
  }
  
  ################################################################################
  # model-implied variance of y. 
  # meanstructure (conditional expectation of y)
  ################################################################################
  tym <- cbind(rep(0,nsamps),samps[,paste0("ty[",1:2,"]")],rep(0,nsamps),samps[,paste0("ty[",3:4,"]")])
  kam <- samps[,paste0("ka[",1:2,"]")]
  muy.1 <- matrix(NA,nsamps,nvar*Nt)
  muy.2 <-  array(NA,c(nvar,Nt,nsamps))
  for(j in 1:nsamps){
    muy.1[j,] <- rep(unlist(tym[j,]),Nt)+(lym.1[,,j])%*%(rep(unlist(kam[j,]),Nt))
    
    for(m in 1:Nt){
      muy.2[,m,j] <- muy.1[j,1:nvar+(m-1)*nvar]
    }
    
  }
  
  
  
  
  ################################################################################################################################
  # Fit indices and chi2
  ################################################################################################################################
  # 1. covariance matrices, mean vectors loglike saturated
  #y.2 <- ydat2
  x_mean.1 <- apply(y.2,2,mean,na.rm=T)
  x_cov.1  <- cov(y.2,use="pair")
  
  x_mean.2 <- array(NA,c(nvar,Nt))
  x_cov.2 <- array(NA,c(nvar,nvar,Nt))
  
  
  ################################################################
  # 2a. global loglikelihood and chisq
  ################################################################
  
  logliki.1 <- array(NA,c(nsamps,N))
  chisqs.1 <- loglik.1 <- rep(NA,nsamps)
  logliksat.1 <- NA
  
  if(id0==T){
    for(j in 1:nsamps){#j<-1
      logliki.1[j,] <- getloglik(y.2,muy.1[j,],sigmay.1[,,j])
    }
    
    loglik.1 <- apply(logliki.1, 1, sum)
    logliksat.1 <- sum(getloglik(y.2,x_mean.1,round(x_cov.1,4))) # loglik sum for all
    chisqs.1 <- getchisq.mcmc(loglik.1,logliksat.1)
  }
  
  
  ################################################################
  # 2b. loglikelihood time-specific (which is extractable from stan)
  ################################################################
  logliki.2 <- array(NA,c(nsamps,Nt,N))
  logliksat.2 <- c()
  loglik.2 <- array(NA,c(nsamps,Nt))
  chisqs.2 <- array(NA,c(nsamps,Nt))
  
  for(j in 1:nsamps){#j<-1
    for(m in 1:Nt){
      logliki.2[j,m,] <- getloglik(y[m,,],muy.2[,m,j],sigmay.2[m,,,j])
      loglik.2[j,m] <- sum(logliki.2[j,m,])
    }
  }
  
  for(m in 1:Nt){#m<-2
    x_cov.2[,,m] <- cov(y[m,,])
    x_mean.2[,m] <- apply(y[m,,],2,mean)
    logliksat.2[m] <- sum(getloglik(y[m,,],x_mean.2[,m],x_cov.2[,,m])) # loglik sum for each Nt separate
    chisqs.2[,m] <- getchisq.mcmc(loglik.2[,m],logliksat.2[m])
  }
  
  ################################################################
  # 3. pd values and other no of vars
  ################################################################
  pd.2 <- pd.2.stan <- c() 
  for(m in 1:Nt){
    loo.2 <- loo(logliki.2[,m,], cores = 4)
    pd.2[m] <- loo.2$estimates["p_loo",1]  ## IMPORTANT!!! FROM R
  }
  
  if(id0==T){
    loo.1 <- loo(logliki.1, cores = 4)## ALSO TRY WITH STAN IMPORT
    pd.1 <- loo.1$estimates["p_loo",1]  ## IMPORTANT!!! FROM R
  }
  
  pstar.1 <- nvar*Nt*((nvar*Nt)+3)/2
  pstar.2 <- nvar*(nvar+3)/2

  p.1 <- nvar*Nt
  p.2 <- nvar

  
  ################################################################
  # 4. actual model fit
  ################################################################
  # version 1 is overall model fit
  if(id0==T){
    BRMSEA.devm.1 <- sqrt(getmax((chisqs.1-pstar.1)/((pstar.1-pd.1)*N)))
    BMc.devm.1       <- exp(-1/(2*N)*(chisqs.1-pstar.1))
    BGamma.devm.1    <- p.1/(p.1+2/N*(chisqs.1-pstar.1))     
    BadjGamma.devm.1 <- 1-pstar.1/(pstar.1-pd.1)*(1-BGamma.devm.1)
  }else{
    BRMSEA.devm.1 <- BMc.devm.1 <- BGamma.devm.1 <- BadjGamma.devm.1 <- rep(NA,nsamps)
  }
  
  # version 2 is time-point specific model fit
  BRMSEA.devm.2 <- BGamma.devm.2 <- BadjGamma.devm.2 <- BMc.devm.2 <- matrix(NA,nsamps,Nt)
  for(m in 1:Nt){
    BRMSEA.devm.2[,m]    <- sqrt(getmax((chisqs.2[,m]-pstar.2)/((pstar.2-pd.2[m])*N)))
    BMc.devm.2[,m]       <- exp(-1/(2*N)*(chisqs.2[,m]-pstar.2))
    BGamma.devm.2[,m]    <- p.2/(p.2+2/N*(chisqs.2[,m]-pstar.2))     
    BadjGamma.devm.2[,m] <- 1-pstar.2/(pstar.2-pd.2[m])*(1-BGamma.devm.2[,m])
  }
  
  
  #### BCFI, BTLI, BNFI 
  BCFI <- 1 - (mean(loglik.1) / logliksat.1)
  
  
  
  out <- list(round(c(mean(BRMSEA.devm.1),mean(BGamma.devm.1),mean(BadjGamma.devm.1),mean(BMc.devm.1)),3),
              round(cbind(apply(BRMSEA.devm.2,2,mean),apply(BGamma.devm.2,2,mean),apply(BadjGamma.devm.2,2,mean),apply(BMc.devm.2,2,mean)),3))
  
  names(out)[[1]] <- c("BRMSEA","BGammahat","adjBGammahat","BMc")
  
  out2 <- list(
    round(mean(BRMSEA.devm.1), 3),        # BRMSEA 
    round(mean(BGamma.devm.1), 3),        # BGamma 
    round(mean(BadjGamma.devm.1), 3),     # BadjGamma 
    round(mean(BMc.devm.1), 3),           # BMc 
    pd.1,                                 # pd value for overall model
    p.1,                                  # number of estimated parameters
    mean(loglik.1),
    logliksat.1,                           # Saturated log likelihood for overall model
    mean(chisqs.1)
  )  
  names(out2) <- c("BRMSEA", "BGammahat", "adjBGammahat", "BMc", "pd_1", "p_1", "loglik", "logliksat_1", "chisqs")
  
  out2

}




