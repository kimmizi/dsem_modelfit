setwd("C:\\holger\\SEM\\modelfit\\stanversion")

library(mvtnorm)
library(matrixcalc)
library(lavaan)

type_MISS <- "tt"
i <- 2

#library(foreach)
#library(doParallel)

#registerDoParallel(8) 

#detectCores()
################################################################
#est <- 
#  foreach (i=1:8) %dopar% {
#    library(mvtnorm)
#    library(matrixcalc)
#    library(lavaan)
    
#current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # RStudio
current_dir <- "/Users/kimzierahn/PycharmProjects/dsem_modelfit/"

    source(file.path(current_dir, "src/lavaan_dsem_models_randomintercept_tt.R"))
    source(file.path(current_dir, "src/lavaan_dsem_models_randomintercept_tt1.R"))
    #source("C:/holger/SEM/modelfit/stanversion/funs/lavaan_dsem_models_comp1.R")
    #source("C:/holger/SEM/modelfit/stanversion/funs/lavaan_dsem_models_comp2.R")
    source(file.path(current_dir, "funs/gen_data_version03.R"))
    source(file.path(current_dir, "funs/gen_data_version04.R"))
    
    ##########################
    #person_size_SIMULATE <- c(91,121, 151, 181,211,501,1001,1501,2001,2501,61,31) 
    person_size_SIMULATE <- c(31, 91) 

    #time_point_SIMULATE <- c(1:5,10,15) # Nt
    time_point_SIMULATE <- c(1, 2)
    
    model_TRUE_MISS_SIMULATE <- c(.3,.6)
    type_TRUE_MISS_SIMULATE <- c("tt","tt1") # within time points, between time points #"none",
    #run_Samples_SIMULATE <- 125
    run_Samples_SIMULATE <- 3
    ##########################
    
    
    ##########################
    phi0 <- diag(2)*.7+.3 # cov(eta)
    mu0  <- c(0,0)        #mean(eta)
    ar0  <- c(.3,.3)      # ar(1) structure
    ly00 <- .6
    ly0  <- matrix(c(ly00,ly00,ly00,0,0,0,
                     0,0,0,ly00,ly00,ly00),6,2,byrow=F) # factor loadings
    td   <- diag(6)*(1-ly00^2) # cond. var(y) -> res var
    ##########################
    
    #summary(sem(dsem[[2]], data=y0,std.lv = TRUE),standardized=T)
    fitnom <- c("npar","fmin","chisq","df","pvalue","baseline.chisq","baseline.df",
                "baseline.pvalue","cfi","tli","nnfi","rfi","nfi","pnfi","ifi","rni",
                "logl","unrestricted.logl","aic","bic","ntotal","bic2","rmsea","rmsea.ci.lower",
                "rmsea.ci.upper","rmsea.ci.level","rmsea.pvalue","rmsea.close.h0","rmsea.notclose.pvalue",
                "rmsea.notclose.h0","rmr","rmr_nomean","srmr","srmr_bentler","srmr_bentler_nomean",
                "crmr","crmr_nomean","srmr_mplus","srmr_mplus_nomean","cn_05","cn_01","gfi",
                "agfi","pgfi","mfi","ecvi") 
    
    #summary(res1,standardized=T)
    # new function so that it does not break the loop
    runstuff <- function(y0,type_MISS){
      
      if(type_MISS=="tt"){
        res1 <- try(sem(dsem.tt[[time_point]], data=y0,std.lv = TRUE,se="none"), silent = TRUE)
      }else{
        res1 <- try(sem(dsem.tt1[[time_point]], data=y0,std.lv = TRUE,se="none"), silent = TRUE)
      }
      
      if(!inherits(res1, "try-error")){
        if(res1@optim$converged==T){
          fitmeasures(res1)  
        }else{rep(NA,46)}
      }else{rep(NA,46)}
    }
    
    #diag((ly0)%*%phi0%*%t(ly0)+td)
    
    #i<-1
    #person_size <- 181
    # time_point <- 15
    #N <- 10000
    for (person_size in person_size_SIMULATE) {
      for (time_point in time_point_SIMULATE) {
        if(time_point*6<person_size){
          
          #for(type_MISS in type_TRUE_MISS_SIMULATE){
            
            #####################################################################
            #set core specific seed
            set.seed(131212023+i+person_size+time_point*1000+as.numeric(as.factor(type_MISS))*10000)
            #####################################################################

            #####################################################################
            if(type_MISS=="none"){#i<-1
              model_TRUE_MISS <- 0
              name_local_SIMULATE_Info <- paste(current_dir, as.character(person_size), 
                                                as.character(time_point),as.character(type_MISS),as.character(model_TRUE_MISS),i ,"_version03_rand", sep = "_")
              
              #######################################
              # empty matrix for results
              #######################################
              #run_Samples_SIMULATE <- SAMPLING <- 1
              fitm1  <- matrix(NA,run_Samples_SIMULATE,46)
              N <- person_size  
              Nt <- time_point
              
              ly1  <- matrix(c(0,0,0,model_TRUE_MISS,0,0,
                               0,0,0,0,0,0),6,2,byrow=F) # factor loadings
              #######################################
              # start looping over 
              for (SAMPLING in 1:run_Samples_SIMULATE){#SAMPLING <-1
                dat1 <- gendata01(N,Nt,phi0,mu0,ar0,ly0,ly1,td)
                
                if(dat1[["is_positive_def"]]==T){
                  y0 <- dat1[["y0"]]
                  # loop with error back-up
                  fitm1[SAMPLING,] <- runstuff(y0)
                }#end pos.def
              }# end looping sampling
              
              colnames(fitm1) <- fitnom
              saveRDS(fitm1, file = paste0(name_local_SIMULATE_Info,"_true.RDS")) # Save samples run
              
              #####################################################################
              #####################################################################
            }else{
              #type_MISS <- "tt1"
              for(model_TRUE_MISS in model_TRUE_MISS_SIMULATE){#model_TRUE_MISS<-.3
                name_local_SIMULATE_Info <- paste(current_dir, as.character(person_size), 
                                                  as.character(time_point),as.character(type_MISS),as.character(model_TRUE_MISS),i ,"_version03_rand", sep = "_")
                
                #######################################
                # empty matrix for results
                #######################################
                fitm1  <- matrix(NA,run_Samples_SIMULATE,46)
                N <- person_size #<-2500  
                Nt <- time_point
                
                if(type_MISS=="tt"){
                  ly1  <- matrix(c(0,0,0,model_TRUE_MISS,0,0,
                                   0,0,0,0,0,0),6,2,byrow=F) # factor loadings
                }else{#model_TRUE_MISS<-.3
                  ly1  <- matrix(c(model_TRUE_MISS,0,0,0,0,0,
                                   0,0,0,0,0,0),6,2,byrow=F) # factor loadings
                }
                #######################################
                # start looping over 
                for (SAMPLING in 1:run_Samples_SIMULATE){#SAMPLING <-1
                  if(type_MISS=="tt"){
                    dat1 <- gendata01(N,Nt,phi0,mu0,ar0,ly0,ly1,td)
                  }else{
                    dat1 <- gendata02(N,Nt,phi0,mu0,ar0,ly0,ly1,td)
                  }
                  
                  
                  if(dat1[["is_positive_def"]]==T){
                    y0 <- dat1[["y0"]]
                    # loop with error back-up
                    fitm1[SAMPLING,] <- runstuff(y0,type_MISS)
                    
                  }#end pos.def
                }# end looping sampling
                
                colnames(fitm1) <- fitnom
                saveRDS(fitm1, file = paste0(name_local_SIMULATE_Info,"_true.RDS")) # Save samples run
                
                #######################################
              }# model_true_miss
            }#end else
            
            
            
            
          #}# end misstype
        }# end of if identification
      }#end time points
    }#end person size
  #}#end for each

