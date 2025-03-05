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

ydat <- dat1$y # use array version of data for stan
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
print(fit1, params2) # parameter estimates -> check convergence


modelfit <- getfit(fit1, ydat, ydat2) # estimate bayesian fit (without cfi/nfi at this point)
#modelfit[[1]] # please ignore the second list element for now
#fitresults1 <- matrix(modelfit[[1]], nrow = 1)
#colnames(fitresults1) <- c("BRMSEA", "BGammahat", "adjBGammahat", "BMc")

fitresults1 <- matrix(modelfit, nrow = 1)
colnames(fitresults1) <-  c("BRMSEA", "BGammahat", "adjBGammahat", "BMc", "pd_1", "p_1", "loglik", "logliksat_1", "chisqs")

fitresults1
# -> results are slightly different than blavaan, check if this is sample size or prior specific? Script should be correct.


# trying out a function fro run_sim_all_models
fitm_stan <- matrix(NA, N_sim_samples, length(fitnom_stan))
fit_st <- fit_model_stan(ydat, ydat2, Timepoints, Person_size)
fit_st <- as.numeric(fit_st)
fitm_stan[1, ] = fit_st

#modelfit <- getfit(fit1, ydat, ydat2)
#fitresults1 <- matrix(modelfit, nrow = 1)
#colnames(fitresults1) <- fitnom_stan
#fitresults1
                    


#################################################
# analyze dat with blavaan
#################################################
library(blavaan)
sem.bl1 <- bsem(dsem[[Timepoints]], ydat2, 
                n.chains = 4, burnin = 1000, sample = 1000) # use similar conditions as stan to make it comparable
summary(sem.bl1,fit=T)
f1.sem.bl1 <- blavFitIndices(sem.bl1, fit.measures = "all") # relevant indices

# all kinds of fit indices
fit_indices <- blavFitIndices(sem.bl1, baseline.model = fit_null) 
BRMSEA <- mean(fit_indices@indices$BRMSEA)
BGammaHat <- mean(fit_indices@indices$BGammaHat)
adjBGammaHat <- mean(fit_indices@indices$adjBGammaHat)
BMc <- mean(fit_indices@indices$BMc)
BCFI <- mean(fit_indices@indices$BCFI)
BTLI <- mean(fit_indices@indices$BTLI)
BNFI <- mean(fit_indices@indices$BNFI)
npar <- blavInspect(sem.bl1, "npar")
marg_loglik <- blavInspect(sem.bl1, "test")[[1]]$stat
ppp <- blavInspect(sem.bl1, "test")[[2]]$stat
chisq <- mean(fit_indices@details$chisq)
pd <- fit_indices@details$pD

# try out function from run_sim_all_models
fitm_blavaan <- matrix(NA, N_sim_samples, length(fitnom_blavaan))
fit_bl = fit_model_blav(ydat2, Timepoints)
fitm_blavaan[1, ] <- fit_bl



#################################################
# experimenting with different null models
#################################################

## fit null model to calculate CFI

fit_0A <- bsem(null_model_0A(Timepoints), ydat2, 
               n.chains = 4, burnin = 1000, sample = 1000)
summary(fit_0A)

fit_0C <- bsem(null_model_0C(Timepoints), ydat2, 
               n.chains = 4, burnin = 1000, sample = 1000)
summary(fit_0C)



#fit_weak_inv <- bsem(weak_invariance_model(Timepoints), ydat2, 
#                     n.chains = 4, burnin = 1000, sample = 1000)

fit_strong <- bsem(strong_invariance_model(Timepoints), ydat2, 
                     n.chains = 4, burnin = 1000, sample = 1000)

fit_strict <- bsem(strict_invariance_model(Timepoints), ydat2, 
                     n.chains = 4, burnin = 1000, sample = 1000)

fit_null <- bsem(null_model(Timepoints), ydat2, 
                 n.chains = 4, burnin = 1000, sample = 1000)

fit_null_res <- bsem(null_model_corr_res(Timepoints), ydat2, 
                 n.chains = 4, burnin = 1000, sample = 1000)


f1.sem.bl1_inv <- blavFitIndices(sem.bl1, baseline.model = fit_conv_inv) # relevant indices
#f1.sem.bl1_weak <- blavFitIndices(sem.bl1, baseline.model = fit_weak_inv) # relevant indices
f1.sem.bl1_strong <- blavFitIndices(sem.bl1, baseline.model = fit_strong) # relevant indices
f1.sem.bl1_strict <- blavFitIndices(sem.bl1, baseline.model = fit_strong) # relevant indices
f1.sem.bl1_null <- blavFitIndices(sem.bl1, baseline.model = fit_null) # relevant indices
f1.sem.bl1_null_res <- blavFitIndices(sem.bl1, baseline.model = fit_null_res, fit.measures = "all") # relevant indices


f1.sem.bl1_inv
#f1.sem.bl1_weak
f1.sem.bl1_strong
f1.sem.bl1_strict
f1.sem.bl1_null
f1.sem.bl1_null_res

#f2.sem.bl1 <- fitmeasures(sem.bl1) # some other fit indices like aic etc., ignore for now




#################################################
# analyze dat with lavaan
#################################################

library(lavaan)
sem.l1 <- sem(dsem[[Timepoints]], ydat2)
fitmeasures(sem.l1) # use similar conditions as stan to make it comparable
summary(sem.l1,fit=T)

fitm_lavaan <- matrix(NA, N_sim_samples, length(fitnom_lavaan))
fit_la = fit_model_lav(ydat2, Timepoints)
fitm_lavaan[1, ] <- fit_la




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




