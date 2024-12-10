# this is a short version for the stan and blavaan file that you can use to integrate in your file.

library(mvtnorm)
library(matrixcalc)
library(lavaan)

# more libraries
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


setwd("/Users/kimzierahn/PycharmProjects/dsem_modelfit/src")


source("lavaan_dsem_models_randomintercept.R")
source("lavaan_dsem_nullmodels.R")
source("gen_data_version03.R")
source("fitfunctions_stan.R")



#################################################
# data
#################################################
Person_size <- c(91)
Timepoints <- 5 # N_timep
Size_crossloading <- c(0,.3,.6)
Type_crossloading <- c("none","tt","tt1")
N_sim_samples <- 3

# DSEM model params
phi0 <- diag(2)*.7+.3 # cov(eta)
mu0  <- c(0, 0)        # mean(eta)
ar0  <- c(.3, .3)      # ar(1) structure
ly00 <- .6
ly0  <- matrix(c(ly00, ly00, ly00, 0, 0, 0,
                 0, 0, 0, ly00, ly00, ly00), 6, 2, byrow = F) # factor loadings
td   <- diag(6)*(1 - ly00^2) # cond. var(y) -> res var



# generate data without misfit here
dat1 <- gendata01(N = Person_size, Nt = Timepoints, phi0, mu0, ar0, ly0, ly1 = matrix(0,6,2), td)
length(dat1)

ydat <- dat1$y #use array version of data for stan
dim(ydat) #NtxN#nitems, ok

ydat2 <- dat1$y0 # use the transformed version in the getfit function (calculate covariance matrix for saturated model)

#################################################
# analyze dat with stan
#################################################
rt1 <- stanc("model3_dsem.stan") #load model
sm1 <- stan_model(stanc_ret = rt1, verbose=FALSE) # compile model

params2 <- c("ly","sigmaeps","sigmaeta1","sigmaeta2","ty","ka","beta1") # parameter names

data.stan1 <- list(N=Person_size,Nt=Timepoints,y=ydat)#,x_cov=round(cov(y.2),4),x_mean=round(apply(y.2,2,mean),4))
fit1 <- sampling(sm1, data=data.stan1)

print(fit1,params2) # parameter estimates -> check convergence
modelfit <- getfit(fit1,ydat,ydat2) # estimate bayesian fit (without cfi/nfi at this point)
#modelfit[[1]] # please ignore the second list element for now
fitresults1 <- matrix(modelfit[[1]],nrow=1)
colnames(fitresults1) <- c("BRMSEA","BGammahat","adjBGammahat","BMc")
fitresults1
# -> results are slightly different than blavaan, check if this is sample size or prior specific? Script should be correct.


#################################################
# analyze dat with blavaan
#################################################
library(blavaan)
sem.bl1 <- bsem(dsem[[Timepoints]], ydat2, 
                n.chains = 4, burnin = 1000, sample = 1000) # use similar conditions as stan to make it comparable
summary(sem.bl1,fit=T)

f1.sem.bl1 <- blavFitIndices(sem.bl1) # relevant indices



#################################################
# experimenting with different null models
#################################################

## fit null model to calculate CFI

fit_conv_inv <- bsem(configural_invariance_model(Timepoints), ydat2, 
                 n.chains = 4, burnin = 1000, sample = 1000)

#fit_weak_inv <- bsem(weak_invariance_model(Timepoints), ydat2, 
#                     n.chains = 4, burnin = 1000, sample = 1000)

fit_strong <- bsem(strong_invariance_model(Timepoints), ydat2, 
                     n.chains = 4, burnin = 1000, sample = 1000)

fit_strict <- bsem(strict_invariance_model(Timepoints), ydat2, 
                     n.chains = 4, burnin = 1000, sample = 1000)

fit_null <- bsem(null_model(Timepoints), ydat2, 
                 n.chains = 4, burnin = 1000, sample = 1000)


f1.sem.bl1_inv <- blavFitIndices(sem.bl1, baseline.model = fit_conv_inv) # relevant indices
#f1.sem.bl1_weak <- blavFitIndices(sem.bl1, baseline.model = fit_weak_inv) # relevant indices
f1.sem.bl1_strong <- blavFitIndices(sem.bl1, baseline.model = fit_strong) # relevant indices
f1.sem.bl1_strict <- blavFitIndices(sem.bl1, baseline.model = fit_strong) # relevant indices
f1.sem.bl1_null <- blavFitIndices(sem.bl1, baseline.model = fit_null) # relevant indices


f1.sem.bl1_inv
#f1.sem.bl1_weak
f1.sem.bl1_strong
f1.sem.bl1_strict
f1.sem.bl1_null


#f2.sem.bl1 <- fitmeasures(sem.bl1) # some other fit indices like aic etc., ignore for now


library(lavaan)
sem.l1 <- sem(dsem[[Timepoints]],ydat2) # use similar conditions as stan to make it comparable
summary(sem.l1,fit=T)








#################################################
# some other helpful things for me

# loglikihood from blavaan (you can extract the likelihood within the getfit function in stan)

i <- 1
cbind(unlist(sem.bl1@external$mcmcout@sim$samples[[1]][paste0("log_lik[",i,"]")])[1:1000+1000],
      unlist(sem.bl1@external$mcmcout@sim$samples[[2]][paste0("log_lik[",i,"]")])[1:1000+1000],
      unlist(sem.bl1@external$mcmcout@sim$samples[[3]][paste0("log_lik[",i,"]")])[1:1000+1000],
      unlist(sem.bl1@external$mcmcout@sim$samples[[4]][paste0("log_lik[",i,"]")])[1:1000+1000])


cbind(unlist(sem.bl1@external$mcmcout@sim$samples[[1]][paste0("log_lik_sat[",i,"]")])[1:1000+1000],
      unlist(sem.bl1@external$mcmcout@sim$samples[[2]][paste0("log_lik_sat[",i,"]")])[1:1000+1000],
      unlist(sem.bl1@external$mcmcout@sim$samples[[3]][paste0("log_lik_sat[",i,"]")])[1:1000+1000],
      unlist(sem.bl1@external$mcmcout@sim$samples[[4]][paste0("log_lik_sat[",i,"]")])[1:1000+1000])


















