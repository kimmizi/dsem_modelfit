setwd("C:\\holger\\SEM\\modelfit\\stanversion")
#setwd("/Users/admin/Desktop/Research Project/CODE/TUUni/SIMULATION")


###################################
###    SIMULATION STARTS HERE
###################################

person_size_SIMULATE <- c(91,121, 151, 181,211,61,31,501,1001,1501,2001,2501) #
time_point_SIMULATE <- c(1:5,10,15,30) # Nt ,30
model_TRUE_MISS_SIMULATE <- c(0,.3,.6)
type_TRUE_MISS_SIMULATE <- c("tt","tt1") # within time points, between time points #"none", HB: taken out because true tt/tt1 model used for type I error rates
run_Samples_SIMULATE <- 125
runpcs <- 2

# FIT INDICES WE ARE ACTUALLY INTERESTED IN
fitnom <- c("cfi","tli","rmsea","srmr","gfi","agfi","cgfi")

cond0 <- data.frame(matrix(NA,1000,4))
res3false <- res3true <- data.frame(matrix(NA,1000,length(fitnom)+2))
res4false <- res4true <- data.frame(matrix(NA,1000,length(fitnom)+2))

# output names
fitnom.res1 <- c("npar","fmin","chisq","df","pvalue","baseline.chisq","baseline.df",
                 "baseline.pvalue","cfi","tli","nnfi","rfi","nfi","pnfi","ifi","rni",
                 "logl","unrestricted.logl","aic","bic","ntotal","bic2","rmsea","rmsea.ci.lower",
                 "rmsea.ci.upper","rmsea.ci.level","rmsea.pvalue","rmsea.close.h0","rmsea.notclose.pvalue",
                 "rmsea.notclose.h0","rmr","rmr_nomean","srmr","srmr_bentler","srmr_bentler_nomean",
                 "crmr","crmr_nomean","srmr_mplus","srmr_mplus_nomean","cn_05","cn_01","gfi",
                 "agfi","pgfi","mfi","ecvi") 

getstuff <- function(name_local_SIMULATE_Info){
  res11 <- try(read.csv(file = name_local_SIMULATE_Info), silent = TRUE)
  if(!inherits(res11, "try-error")){
    res11
  }else{
    matrix(NA,run_Samples_SIMULATE,length(fitnom.res1))    
  }
}

current_dir <- "~/PycharmProjects/dsem_modelfit/src/"

ijk <- 1
for (time_point in time_point_SIMULATE) {
  for (person_size in person_size_SIMULATE) {
    for(type_MISS in type_TRUE_MISS_SIMULATE){
      for (model_TRUE_MISS in model_TRUE_MISS_SIMULATE){
        
        N <- person_size#<-211  # persons
        Nt <- time_point#<-2 #time points
        
        
        if(time_point*6<person_size){
          
          res1 <- res2 <- data.frame(matrix(NA,1,length(fitnom.res1)))
          colnames(res1) <- colnames(res2) <- fitnom.res1
          for(i in 1:runpcs){
            
            # MISSSPECIFICATIONS
            name_local_SIMULATE_Info <- paste(current_dir, as.character(person_size), as.character(time_point), as.character(type_MISS), 
                                              as.character(model_TRUE_MISS), core, "_version03_rand", sep = "_")
            
            # TRUE DGP
            name_local_SIMULATE_Info2 <- paste("C:\\holger\\SEM\\modelfit\\stanversion\\results_lavaan_version03_rand\\local", as.character(person_size), 
                                              as.character(time_point),as.character(type_MISS),as.character(model_TRUE_MISS),i ,"_version03_rand_true.RDS", sep = "_")
            
            
            res11 <- getstuff(name_local_SIMULATE_Info)
            res21 <- getstuff(name_local_SIMULATE_Info2)
            colnames(res21) <- colnames(res11) <- fitnom.res1

            res1 <- rbind(res1,res11)
            res2 <- rbind(res2,res21)
          }
          #colnames(res11) <- colnames(res1m2) <- colnames(res1m2) <- colnames(res0)
          #colnames(res1m12) <- colnames(res1m13) <- c("chi2","chi2_p")
          
          res1 <- data.frame(res1[-1,])
          res2 <- data.frame(res2[-1,])
          
          res1$cgfi <- res1$gfi+(Nt*6+1)*Nt*6/res1$npar/N
          res2$cgfi <- res2$gfi+(Nt*6+1)*Nt*6/res2$npar/N
          
          #dim(res1)
          # average scores m false (except none)
          # FIT INDICES FOR WRONG MODEL
          # AVERAGE SCORES
          res3false[ijk,1] <- mean(res1$chisq/res1$df,na.rm=T) # CHI2
          res3false[ijk,1:length(fitnom)+1] <- apply(res1[fitnom],2,mean,na.rm=T)
          res3false[ijk,length(fitnom)+2] <- length(na.omit(res1$npar)) # LENGHT OF MISSING DATA
          
          # average scores m true
          # FIT INDICES FOR TRUE MODEL
          res3true[ijk,1] <- mean(res2$chisq/res1$df,na.rm=T)
          res3true[ijk,1:length(fitnom)+1] <- apply(res2[fitnom],2,mean,na.rm=T)
          res3true[ijk,length(fitnom)+2] <- length(na.omit(res2$npar))
          
          # CUT OFFS
          #cutoffs good fit (misfit): percentage of rejection
          res4false[ijk,1] <- 1-mean(res1$pvalue>.05,na.rm=T) # if p-value is larger than 5% then its 1 = percentage of rejection = POWER
          res4false[ijk,2] <- 1-mean(res1$chisq < 2*res1$df,na.rm=T)
          res4false[ijk,3] <- 1-mean(res1$cfi > .95,na.rm=T)
          res4false[ijk,4] <- 1-mean(res1$tli > .97,na.rm=T)
          res4false[ijk,5] <- 1-mean(res1$rmsea < .05,na.rm=T)
          res4false[ijk,6] <- 1-mean(res1$srmr < .05,na.rm=T)
          res4false[ijk,7] <- 1-mean(res1$gfi > .95,na.rm=T)
          res4false[ijk,8] <- 1-mean(res1$agfi > .95,na.rm=T)
          res4false[ijk,9] <- 1-mean(res1$cgfi > .95,na.rm=T)
          
          res4true[ijk,1] <- 1-mean(res2$pvalue>.05,na.rm=T) # TYPE 1 ERROR RATES
          res4true[ijk,2] <- 1-mean(res2$chisq < 2*res1$df,na.rm=T)
          res4true[ijk,3] <- 1-mean(res2$cfi > .95,na.rm=T)
          res4true[ijk,4] <- 1-mean(res2$tli > .97,na.rm=T)
          res4true[ijk,5] <- 1-mean(res2$rmsea < .05,na.rm=T)
          res4true[ijk,6] <- 1-mean(res2$srmr < .05,na.rm=T)
          res4true[ijk,7] <- 1-mean(res2$gfi > .95,na.rm=T)
          res4true[ijk,8] <- 1-mean(res2$agfi > .95,na.rm=T)
          res4true[ijk,9] <- 1-mean(res2$cgfi > .95,na.rm=T)
          
          # good fit is enough:
          #cutoffs acceptable fit
          #res4[ijk,1] <- 1-mean(res1$pvalue>.05,na.rm=T)
          #res4[ijk,2] <- 1-mean(res1$chisq < 3*res1$df,na.rm=T)
          #res4[ijk,3] <- 1-mean(res1$cfi > .90,na.rm=T)
          #res4[ijk,4] <- 1-mean(res1$tli > .95,na.rm=T)
          #res4[ijk,5] <- 1-mean(res1$rmsea < .08,na.rm=T)
          #res4[ijk,6] <- 1-mean(res1$srmr < .08,na.rm=T)
          
          
          
          cond0[ijk,] <- c(time_point,person_size,type_MISS,model_TRUE_MISS)
          ijk <- ijk+1 # adding a new line
        }
      }
    }
  }
}

#save all data
res3false2 <- res3false
res4false2 <- res4false
res3true2  <- res3true
res4true2  <- res4true

# omit missing data
res3false1 <- na.omit(cbind(cond0,res3false))[,-(1:4)]
res4false1 <- na.omit(cbind(cond0,res4false))[,-(1:4)]
res3true1 <- na.omit(cbind(cond0,res3true))[,-(1:4)]
res4true1 <- na.omit(cbind(cond0,res4true))[,-(1:4)]
cond1 <- na.omit(cbind(cond0,res4false))[,(1:4)]
cond2 <- na.omit(cbind(cond0,res4true))[,(1:4)]

# check dimensions
dim(cond1)
dim(cond2)

# put col names
res3false <- res3false1
res4false <- res4false1
res3true  <- res3true1
res4true  <- res4true1
cond0 <- cond1

colnames(cond0) <- c("Nt.f","N.f","type","miss.f")

# for plotting:
# as continous measure:
cond0$Nt <- as.numeric(cond0$Nt.f)
cond0$N <- as.numeric(cond0$N.f)
cond0$miss <- as.numeric(cond0$miss.f)

# as factor
cond0$Nt.f <- factor(cond0$Nt,levels = c("1","2","3","4","5","10","15","30"))
cond0$N.f <-  factor(cond0$N,levels=paste0(sort(person_size_SIMULATE)))
cond0$type <- as.factor(cond0$type)


colnames(res3false) <- colnames(res3true) <- c("chi2/df",fitnom,"nconv")
colnames(res4false) <- colnames(res4true) <- c("chi2_p","chi2/df",fitnom)

#cbind(cond0,res2)
cbind(cond0,res3false)
cbind(cond0,res4false)

range(res3false$nconv)
hist(res3false$nconv)
range(res3true$nconv)
hist(res3true$nconv)


##########################################################################
# x-axis is N and separate lines for Nt
##########################################################################

### PLOTTING, developing own plot
plotfun <- function(dat,miss0,misstyp,who,woline=.8,ord=F){ # who selects one of fit indices, wholine = 0.8 for power
  #dat<-res3
  #miss0<-0
  #who<-"chi2_p"
  #misstyp<-"none"
  par(new=F)
  j<-1
  cond11 <- cond0$miss==miss0 & cond0$type==misstyp & cond0$Nt.f==levels(cond0$Nt.f)[j]
  
  if(ord==F){
    plot(sort(cond0$N[cond11]),dat[cond11,who][rank(cond0$N[cond11])],ylim=c(0,1),xlim=c(0,max(cond0$N)),col=j,pch=j,type="b",
         axes=F,xlab="N",ylab="% rejected",main=who)
  }else{
    plot(sort(log(cond0$N[cond11])),dat[cond11,who][rank(cond0$N[cond11])],ylim=c(0,1),xlim=c(log(min(cond0$N)),log(max(cond0$N))),col=j,pch=j,type="b",
         axes=F,xlab="N",ylab="% rejected",main=who)
  }
  
  for(j in 2:length(levels(cond0$Nt.f))){#j <- length(levels(cond0$Nt.f))
    par(new=T)
    cond11 <- cond0$miss==miss0 & cond0$type==misstyp & cond0$Nt.f==levels(cond0$Nt.f)[j]
    
    if(ord==F){
      plot(sort(cond0$N[cond11]),dat[cond11,who][rank(cond0$N[cond11])],ylim=c(0,1),xlim=c(0,max(cond0$N)),col=j,pch=j,type="b",
           axes=F,ylab="",xlab="")
    }else{
      plot(sort(log(cond0$N[cond11])),dat[cond11,who][rank(cond0$N[cond11])],ylim=c(0,1),xlim=c(log(min(cond0$N)),log(max(cond0$N))),col=j,pch=j,type="b",
           axes=F,ylab="",xlab="")
    }
    
  }
  
  if(ord==F){
    axis(1,person_size_SIMULATE)
  }else{
    axis(1,at=log(person_size_SIMULATE),labels=person_size_SIMULATE)
  }
  
  axis(2)
  abline(h=woline,lty=3)
  
  legend("topright",paste0("N_t=",levels(cond0$Nt.f)),lty=1,col=1:length(levels(cond0$Nt.f)),pch=1:length(levels(cond0$Nt.f)))
}


#######################################################################
pdf("type1error01_v1.pdf",height=3*4,width=2*4) # create PDF
par(mfrow=c(3,2))

# LABELING:
plotfun(res3,miss0=0,misstyp="none",who="chi2_p",woline=.05,ord=T)
plotfun(res3,miss0=0,misstyp="none",who="chi2/df",woline=.05,ord=T)
plotfun(res3,miss0=0,misstyp="none",who="rmsea",woline=.05,ord=T)
plotfun(res3,miss0=0,misstyp="none",who="srmr",woline=.05,ord=T)
plotfun(res3,miss0=0,misstyp="none",who="cfi",woline=.05,ord=T)
plotfun(res3,miss0=0,misstyp="none",who="tli",woline=.05,ord=T)
dev.off()

#######################################################################
pdf("power_small_within_v1.pdf",height=3*4,width=2*4)
par(mfrow=c(3,2))
plotfun(res3,miss0=0.3,misstyp="tt",who="chi2_p",woline=.80,ord=T)
plotfun(res3,miss0=0.3,misstyp="tt",who="chi2/df",woline=.80,ord=T)
plotfun(res3,miss0=0.3,misstyp="tt",who="rmsea",woline=.80,ord=T)
plotfun(res3,miss0=0.3,misstyp="tt",who="srmr",woline=.80,ord=T)
plotfun(res3,miss0=0.3,misstyp="tt",who="cfi",woline=.80,ord=T)
plotfun(res3,miss0=0.3,misstyp="tt",who="tli",woline=.80,ord=T)
dev.off()

pdf("power_large_within_v1.pdf",height=3*4,width=2*4)
par(mfrow=c(3,2))
plotfun(res3,miss0=0.6,misstyp="tt",who="chi2_p",woline=.80,ord=T)
plotfun(res3,miss0=0.6,misstyp="tt",who="chi2/df",woline=.80,ord=T)
plotfun(res3,miss0=0.6,misstyp="tt",who="rmsea",woline=.80,ord=T)
plotfun(res3,miss0=0.6,misstyp="tt",who="srmr",woline=.80,ord=T)
plotfun(res3,miss0=0.6,misstyp="tt",who="cfi",woline=.80,ord=T)
plotfun(res3,miss0=0.6,misstyp="tt",who="tli",woline=.80,ord=T)
dev.off()

#######################################################################
pdf("power_small_across_v1.pdf",height=3*4,width=2*4)
# note that N_t=1 is a type I error rate without misspecification
par(mfrow=c(3,2))
plotfun(res3,miss0=0.3,misstyp="tt1",who="chi2_p",woline=.80)
plotfun(res3,miss0=0.3,misstyp="tt1",who="chi2/df",woline=.80)
plotfun(res3,miss0=0.3,misstyp="tt1",who="rmsea",woline=.80)
plotfun(res3,miss0=0.3,misstyp="tt1",who="srmr",woline=.80)
plotfun(res3,miss0=0.3,misstyp="tt1",who="cfi",woline=.80)
plotfun(res3,miss0=0.3,misstyp="tt1",who="tli",woline=.80)
dev.off()

pdf("power_large_across_v1.pdf",height=3*4,width=2*4)
par(mfrow=c(3,2))
plotfun(res3,miss0=0.6,misstyp="tt1",who="chi2_p",woline=.80)
plotfun(res3,miss0=0.6,misstyp="tt1",who="chi2/df",woline=.80)
plotfun(res3,miss0=0.6,misstyp="tt1",who="rmsea",woline=.80)
plotfun(res3,miss0=0.6,misstyp="tt1",who="srmr",woline=.80)
plotfun(res3,miss0=0.6,misstyp="tt1",who="cfi",woline=.80)
plotfun(res3,miss0=0.6,misstyp="tt1",who="tli",woline=.80)
dev.off()

#######################################################################





##########################################################################
# x-axis is Nt and separate lines for N, looks better
##########################################################################
plotfun2 <- function(dat,miss0,misstyp,who,woline=.8,not1=F){
  #dat<-res2
  #miss0<-0
  #who<-"chi2_p"
  par(new=F)
  j<-1
  
  if(not1==T){
    cond11 <- cond0$miss==miss0 & cond0$type==misstyp & cond0$N.f==levels(cond0$N.f)[j] & cond0$Nt!=1
  }else{
    cond11 <- cond0$miss==miss0 & cond0$type==misstyp & cond0$N.f==levels(cond0$N.f)[j]
  }
  
  #l0 <- dat[cond0$miss==miss0&cond0$N1==levels(cond0$N1)[j],who]
  
  plot(cond0$Nt[cond11],dat[cond11,who],
       #l0+runif(length(l0),-.025,.025),
       ylim=c(0,1),xlim=c(0,15),col=j,pch=j,type="b",
       axes=F,xlab="Nt",ylab="% rejected",main=who)
  
  for(j in 2:length(levels(cond0$N.f))){#j<-2
    #l0 <- dat[cond0$miss==miss0&cond0$N1==levels(cond0$N1)[j],who]
    if(not1==T){
      cond11 <- cond0$miss==miss0 & cond0$type==misstyp & cond0$N.f==levels(cond0$N.f)[j] & cond0$Nt!=1
    }else{
      cond11 <- cond0$miss==miss0 & cond0$type==misstyp & cond0$N.f==levels(cond0$N.f)[j]
    }
    
    par(new=T)
    plot(cond0$Nt[cond11],dat[cond11,who],ylim=c(0,1),xlim=c(0,15),col=j,pch=j,type="b",
         axes=F,ylab="",xlab="")
  }
  axis(1,time_point_SIMULATE)
  axis(2)
  abline(h=woline,lty=3)
  
  legend("topleft",paste0("N=",levels(cond0$N.f)),lty=1,col=1:length(levels(cond0$N.f)),pch=1:length(levels(cond0$N.f)))
}

#######################################################################

# THE PLOT HE USED IN THE END:
pdf("type1error01_v2.pdf",height=3*4,width=2*4)
par(mfrow=c(3,2))
plotfun2(res3,miss0=0.0,misstyp="none",who="chi2_p",woline=.05)
plotfun2(res3,miss0=0.0,misstyp="none",who="chi2/df",woline=.05)
plotfun2(res3,miss0=0.0,misstyp="none",who="rmsea",woline=.05)
plotfun2(res3,miss0=0.0,misstyp="none",who="srmr",woline=.05)
plotfun2(res3,miss0=0.0,misstyp="none",who="cfi",woline=.05)
plotfun2(res3,miss0=0.0,misstyp="none",who="tli",woline=.05)
dev.off()

#######################################################################
pdf("power_small_within_v2.pdf",height=3*4,width=2*4)
par(mfrow=c(3,2))
plotfun2(res3,miss0=0.3,misstyp="tt",who="chi2_p")
plotfun2(res3,miss0=0.3,misstyp="tt",who="chi2/df")
plotfun2(res3,miss0=0.3,misstyp="tt",who="rmsea")
plotfun2(res3,miss0=0.3,misstyp="tt",who="srmr")
plotfun2(res3,miss0=0.3,misstyp="tt",who="cfi")
plotfun2(res3,miss0=0.3,misstyp="tt",who="tli")
dev.off()

pdf("power_large_within_v2.pdf",height=3*4,width=2*4)
par(mfrow=c(3,2))
plotfun2(res3,miss0=0.6,misstyp="tt",who="chi2_p")
plotfun2(res3,miss0=0.6,misstyp="tt",who="chi2/df")
plotfun2(res3,miss0=0.6,misstyp="tt",who="rmsea")
plotfun2(res3,miss0=0.6,misstyp="tt",who="srmr")
plotfun2(res3,miss0=0.6,misstyp="tt",who="cfi")
plotfun2(res3,miss0=0.6,misstyp="tt",who="tli")
dev.off()

#######################################################################
pdf("power_small_across_v2.pdf",height=3*4,width=2*4)
par(mfrow=c(3,2))
plotfun2(res3,miss0=0.3,misstyp="tt1",who="chi2_p",not1=T)
plotfun2(res3,miss0=0.3,misstyp="tt1",who="chi2/df",not1=T)
plotfun2(res3,miss0=0.3,misstyp="tt1",who="rmsea",not1=T)
plotfun2(res3,miss0=0.3,misstyp="tt1",who="srmr",not1=T)
plotfun2(res3,miss0=0.3,misstyp="tt1",who="cfi",not1=T)
plotfun2(res3,miss0=0.3,misstyp="tt1",who="tli",not1=T)
dev.off()

pdf("power_large_across_v2.pdf",height=3*4,width=2*4)
par(mfrow=c(3,2))
plotfun2(res3,miss0=0.6,misstyp="tt1",who="chi2_p",not1=T)
plotfun2(res3,miss0=0.6,misstyp="tt1",who="chi2/df",not1=T)
plotfun2(res3,miss0=0.6,misstyp="tt1",who="rmsea",not1=T)
plotfun2(res3,miss0=0.6,misstyp="tt1",who="srmr",not1=T)
plotfun2(res3,miss0=0.6,misstyp="tt1",who="cfi",not1=T)
plotfun2(res3,miss0=0.6,misstyp="tt1",who="tli",not1=T)
dev.off()
#######################################################################

##########################################################################
# x-axis is N and separate lines for Nt
# plots type I error rate and power in the same picture, too crowded
##########################################################################
plotfun3 <- function(dat,miss0,misstyp,who,woline=.8,not1=F){
  #dat<-res2
  #miss0<-0
  #who<-"chi2_p"
  par(new=F)
  j<-1
  
  if(not1==T){
    cond01 <- cond0$miss==0 & cond0$type=="none" & cond0$N.f==levels(cond0$N.f)[j] & cond0$Nt!=1
    cond11 <- cond0$miss==miss0 & cond0$type==misstyp & cond0$N.f==levels(cond0$N.f)[j] & cond0$Nt!=1
  }else{
    cond01 <- cond0$miss==0 & cond0$type=="none" & cond0$N.f==levels(cond0$N.f)[j]
    cond11 <- cond0$miss==miss0 & cond0$type==misstyp & cond0$N.f==levels(cond0$N.f)[j]
  }
  
  #l0 <- dat[cond0$miss==miss0&cond0$N1==levels(cond0$N1)[j],who]
  
  plot(cond0$Nt[cond01],dat[cond01,who],ylim=c(0,1),xlim=c(0,15),col=j,pch=j,type="b",axes=F,xlab="Nt",ylab="% rejected",main=who,lty=2)
  par(new=T)
  plot(cond0$Nt[cond11],dat[cond11,who],ylim=c(0,1),xlim=c(0,15),col=j,pch=j,type="b",axes=F,xlab="",ylab="")
  
  for(j in 2:length(levels(cond0$N.f))){#j<-2
    #l0 <- dat[cond0$miss==miss0&cond0$N1==levels(cond0$N1)[j],who]
    if(not1==T){
      cond01 <- cond0$miss==0 & cond0$type=="none" & cond0$N.f==levels(cond0$N.f)[j] & cond0$Nt!=1
      cond11 <- cond0$miss==miss0 & cond0$type==misstyp & cond0$N.f==levels(cond0$N.f)[j] & cond0$Nt!=1
    }else{
      cond01 <- cond0$miss==0 & cond0$type=="none" & cond0$N.f==levels(cond0$N.f)[j]
      cond11 <- cond0$miss==miss0 & cond0$type==misstyp & cond0$N.f==levels(cond0$N.f)[j]
    }
    
    par(new=T)
    plot(cond0$Nt[cond01],dat[cond01,who],ylim=c(0,1),xlim=c(0,15),col=j,pch=j,type="b",axes=F,xlab="",ylab="",lty=2)
    par(new=T)
    plot(cond0$Nt[cond11],dat[cond11,who],ylim=c(0,1),xlim=c(0,15),col=j,pch=j,type="b",axes=F,xlab="",ylab="")
    
  }
  axis(1,time_point_SIMULATE)
  axis(2)
  abline(h=woline,lty=3)
  
  legend("topleft",paste0("N=",levels(cond0$N.f)),lty=1,col=1:length(levels(cond0$N.f)),pch=1:length(levels(cond0$N.f)))
}

#######################################################################

# TYPE 1 ERROR RATE AND POWER
pdf("both_small_within_v2.pdf",height=3*4,width=2*4)
par(mfrow=c(3,2))
plotfun3(res3,miss0=0.3,misstyp="tt",who="chi2_p")
plotfun3(res3,miss0=0.3,misstyp="tt",who="chi2/df")
plotfun3(res3,miss0=0.3,misstyp="tt",who="rmsea")
plotfun3(res3,miss0=0.3,misstyp="tt",who="srmr")
plotfun3(res3,miss0=0.3,misstyp="tt",who="cfi")
plotfun3(res3,miss0=0.3,misstyp="tt",who="tli")
dev.off()

pdf("both_large_within_v2.pdf",height=3*4,width=2*4)
par(mfrow=c(3,2))
plotfun3(res3,miss0=0.6,misstyp="tt",who="chi2_p")
plotfun3(res3,miss0=0.6,misstyp="tt",who="chi2/df")
plotfun3(res3,miss0=0.6,misstyp="tt",who="rmsea")
plotfun3(res3,miss0=0.6,misstyp="tt",who="srmr")
plotfun3(res3,miss0=0.6,misstyp="tt",who="cfi")
plotfun3(res3,miss0=0.6,misstyp="tt",who="tli")
dev.off()

#######################################################################
pdf("both_small_across_v2.pdf",height=3*4,width=2*4)
par(mfrow=c(3,2))
plotfun3(res3,miss0=0.3,misstyp="tt1",who="chi2_p",not1=T)
plotfun3(res3,miss0=0.3,misstyp="tt1",who="chi2/df",not1=T)
plotfun3(res3,miss0=0.3,misstyp="tt1",who="rmsea",not1=T)
plotfun3(res3,miss0=0.3,misstyp="tt1",who="srmr",not1=T)
plotfun3(res3,miss0=0.3,misstyp="tt1",who="cfi",not1=T)
plotfun3(res3,miss0=0.3,misstyp="tt1",who="tli",not1=T)
dev.off()

pdf("both_large_across_v2.pdf",height=3*4,width=2*4)
par(mfrow=c(3,2))
plotfun3(res3,miss0=0.6,misstyp="tt1",who="chi2_p",not1=T)
plotfun3(res3,miss0=0.6,misstyp="tt1",who="chi2/df",not1=T)
plotfun3(res3,miss0=0.6,misstyp="tt1",who="rmsea",not1=T)
plotfun3(res3,miss0=0.6,misstyp="tt1",who="srmr",not1=T)
plotfun3(res3,miss0=0.6,misstyp="tt1",who="cfi",not1=T)
plotfun3(res3,miss0=0.6,misstyp="tt1",who="tli",not1=T)
dev.off()
#######################################################################

##########################################################################
# x-axis is N and separate lines for Nt
# plots difference between misfit and type I error rate
# use for both average values and rejection rates
##########################################################################


# X-AXIS AS TIME POINTS OR PERSONS?

plotfun4 <- function(dat,miss0,misstyp,who,woline=.0,not1=F,ylim=ylim){
  #dat<-res2
  #miss0<-0.3
  #who<-"chi2/df"
  #misstyp<-"tt"
  par(new=F)
  j<-1
  
  if(not1==T){
    cond01 <- cond0$miss==0 & cond0$type=="none" & cond0$N.f==levels(cond0$N.f)[j] & cond0$Nt!=1
    cond11 <- cond0$miss==miss0 & cond0$type==misstyp & cond0$N.f==levels(cond0$N.f)[j] & cond0$Nt!=1
  }else{
    cond01 <- cond0$miss==0 & cond0$type=="none" & cond0$N.f==levels(cond0$N.f)[j]
    cond11 <- cond0$miss==miss0 & cond0$type==misstyp & cond0$N.f==levels(cond0$N.f)[j]
  }
  
  #l0 <- dat[cond0$miss==miss0&cond0$N1==levels(cond0$N1)[j],who]
  
  plot(cond0$Nt[cond01],dat[cond11,who]-dat[cond01,who],ylim=ylim,xlim=c(0,15),col=j,pch=j,type="b",axes=F,xlab="Nt",ylab="Difference",main=who)
  
  for(j in 2:length(levels(cond0$N.f))){#j<-2
    #l0 <- dat[cond0$miss==miss0&cond0$N1==levels(cond0$N1)[j],who]
    if(not1==T){
      cond01 <- cond0$miss==0 & cond0$type=="none" & cond0$N.f==levels(cond0$N.f)[j] & cond0$Nt!=1
      cond11 <- cond0$miss==miss0 & cond0$type==misstyp & cond0$N.f==levels(cond0$N.f)[j] & cond0$Nt!=1
    }else{
      cond01 <- cond0$miss==0 & cond0$type=="none" & cond0$N.f==levels(cond0$N.f)[j]
      cond11 <- cond0$miss==miss0 & cond0$type==misstyp & cond0$N.f==levels(cond0$N.f)[j]
    }
    
    par(new=T)
    plot(cond0$Nt[cond01],dat[cond11,who]-dat[cond01,who],ylim=ylim,xlim=c(0,15),col=j,pch=j,type="b",axes=F,xlab="",ylab="")
    
  }
  axis(1,time_point_SIMULATE)
  axis(2)
  abline(h=woline,lty=3)
  
  legend("topleft",paste0("N=",levels(cond0$N.f)),lty=1,col=1:length(levels(cond0$N.f)),pch=1:length(levels(cond0$N.f)))
}

# difference between true and misspecified data
pdf("diff_small_across_v2.pdf",height=3*4,width=2*4)
par(mfrow=c(3,2))
plotfun4(res3,miss0=0.3,misstyp="tt1",who="chi2_p",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.3,misstyp="tt1",who="chi2/df",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.3,misstyp="tt1",who="rmsea",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.3,misstyp="tt1",who="srmr",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.3,misstyp="tt1",who="cfi",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.3,misstyp="tt1",who="tli",not1=T,ylim=c(0,1))

plotfun4(res2,miss0=0.3,misstyp="tt1",who="chi2/df",not1=T,ylim=c(0,1))
plotfun4(res2,miss0=0.3,misstyp="tt1",who="rmsea",not1=T,ylim=c(0,.1))
plotfun4(res2,miss0=0.3,misstyp="tt1",who="srmr",not1=T,ylim=c(0,.1))
plotfun4(res2,miss0=0.3,misstyp="tt1",who="cfi",not1=T,ylim=c(-.3,.1))
plotfun4(res2,miss0=0.3,misstyp="tt1",who="tli",not1=T,ylim=c(-.3,.1))
dev.off()

pdf("diff_large_across_v2.pdf",height=3*4,width=2*4)
par(mfrow=c(3,2))
plotfun4(res3,miss0=0.6,misstyp="tt1",who="chi2_p",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.6,misstyp="tt1",who="chi2/df",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.6,misstyp="tt1",who="rmsea",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.6,misstyp="tt1",who="srmr",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.6,misstyp="tt1",who="cfi",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.6,misstyp="tt1",who="tli",not1=T,ylim=c(0,1))

plotfun4(res2,miss0=0.6,misstyp="tt1",who="chi2/df",not1=T,ylim=c(0,1))
plotfun4(res2,miss0=0.6,misstyp="tt1",who="rmsea",not1=T,ylim=c(0,.1))
plotfun4(res2,miss0=0.6,misstyp="tt1",who="srmr",not1=T,ylim=c(0,.1))
plotfun4(res2,miss0=0.6,misstyp="tt1",who="cfi",not1=T,ylim=c(-.3,.1))
plotfun4(res2,miss0=0.6,misstyp="tt1",who="tli",not1=T,ylim=c(-.3,.1))
dev.off()

pdf("diff_small_within_v2.pdf",height=3*4,width=2*4)
par(mfrow=c(3,2))
plotfun4(res3,miss0=0.3,misstyp="tt",who="chi2_p",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.3,misstyp="tt",who="chi2/df",not1=T,ylim=c(0,.1))
plotfun4(res3,miss0=0.3,misstyp="tt",who="rmsea",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.3,misstyp="tt",who="srmr",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.3,misstyp="tt",who="cfi",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.3,misstyp="tt",who="tli",not1=T,ylim=c(0,1))

plotfun4(res2,miss0=0.3,misstyp="tt",who="chi2/df",not1=T,ylim=c(0,1))
plotfun4(res2,miss0=0.3,misstyp="tt",who="rmsea",not1=T,ylim=c(0,.1))
plotfun4(res2,miss0=0.3,misstyp="tt",who="srmr",not1=T,ylim=c(0,.1))
plotfun4(res2,miss0=0.3,misstyp="tt",who="cfi",not1=T,ylim=c(-.3,.1))
plotfun4(res2,miss0=0.3,misstyp="tt",who="tli",not1=T,ylim=c(-.3,.1))
dev.off()

pdf("diff_large_within_v2.pdf",height=3*4,width=2*4)
par(mfrow=c(3,2))
plotfun4(res3,miss0=0.6,misstyp="tt",who="chi2_p",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.6,misstyp="tt",who="chi2/df",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.6,misstyp="tt",who="rmsea",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.6,misstyp="tt",who="srmr",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.6,misstyp="tt",who="cfi",not1=T,ylim=c(0,1))
plotfun4(res3,miss0=0.6,misstyp="tt",who="tli",not1=T,ylim=c(0,1))

plotfun4(res2,miss0=0.6,misstyp="tt",who="chi2/df",not1=T,ylim=c(0,1))
plotfun4(res2,miss0=0.6,misstyp="tt",who="rmsea",not1=T,ylim=c(0,.1))
plotfun4(res2,miss0=0.6,misstyp="tt",who="srmr",not1=T,ylim=c(0,.1))
plotfun4(res2,miss0=0.6,misstyp="tt",who="cfi",not1=T,ylim=c(-.3,.1))
plotfun4(res2,miss0=0.6,misstyp="tt",who="tli",not1=T,ylim=c(-.3,.1))
dev.off()

################################################################################
################################################################################
################################################################################
# tables
################################################################################
################################################################################
################################################################################
library(xtable)

# CREATE TABLES
tabdat1 <- function(res,index,type,miss){
  #type <- "none";miss=0 
  #res <- res3
  #index <- "cfi"
  cond1 <- cond0$type==type & cond0$miss==miss
  dat1 <- data.frame(matrix(NA,length(time_point_SIMULATE),length(person_size_SIMULATE)))
  colnames(dat1) <- sort(person_size_SIMULATE)
  rownames(dat1) <- sort(time_point_SIMULATE)
  
  for(j in 1:length(person_size_SIMULATE)){
    for(k in 1:length(time_point_SIMULATE)){
      if(length(res[cond1&cond0$N==sort(person_size_SIMULATE)[j]&cond0$Nt==sort(time_point_SIMULATE)[k],index])>0){
        dat1[k,j] <- res[cond1&cond0$N==sort(person_size_SIMULATE)[j]&cond0$Nt==sort(time_point_SIMULATE)[k],index]
      }
    }
  }
  dat1
}

tabdat1(res3false,"cfi","tt",0.3)
tabdat1(res3true,"cfi","tt",0.3)

#type I error rate by sample size
print(xtable(tabdat1(res4true,"cfi","tt",0.3)), include.rownames=T)
print(xtable(tabdat1(res4true,"tli","tt",0.3)), include.rownames=T)
print(xtable(tabdat1(res4true,"rmsea","tt",0.3)), include.rownames=T)
print(xtable(tabdat1(res4true,"srmr","tt",0.3)), include.rownames=T)
print(xtable(tabdat1(res4true,"chi2/df","tt",0.3)), include.rownames=T)
print(xtable(tabdat1(res4true,"chi2_p","tt",0.3)), include.rownames=T)

#power tt=.3 rate by sample size
print(xtable(tabdat1(res4false,"cfi","tt",0.3)), include.rownames=T)
print(xtable(tabdat1(res4false,"tli","tt",0.3)), include.rownames=T)
print(xtable(tabdat1(res4false,"rmsea","tt",0.3)), include.rownames=T)
print(xtable(tabdat1(res4false,"srmr","tt",0.3)), include.rownames=T)
print(xtable(tabdat1(res4false,"chi2/df","tt",0.3)), include.rownames=T)
print(xtable(tabdat1(res4false,"chi2_p","tt",0.3)), include.rownames=T)


#type I error rate by sample size
print(xtable(tabdat1(res4true,"cfi","tt1",0.6)), include.rownames=T)
print(xtable(tabdat1(res4true,"tli","tt1",0.6)), include.rownames=T)
print(xtable(tabdat1(res4true,"rmsea","tt1",0.6)), include.rownames=T)
print(xtable(tabdat1(res4true,"srmr","tt1",0.6)), include.rownames=T)
print(xtable(tabdat1(res4true,"chi2/df","tt1",0.6)), include.rownames=T)
print(xtable(tabdat1(res4true,"chi2_p","tt1",0.6)), include.rownames=T)

#power tt=.6 rate by sample size
print(xtable(tabdat1(res4false,"cfi","tt",0.6)), include.rownames=T)
print(xtable(tabdat1(res4false,"tli","tt",0.6)), include.rownames=T)
print(xtable(tabdat1(res4false,"rmsea","tt",0.6)), include.rownames=T)
print(xtable(tabdat1(res4false,"srmr","tt",0.6)), include.rownames=T)
print(xtable(tabdat1(res4false,"chi2/df","tt",0.6)), include.rownames=T)
print(xtable(tabdat1(res4false,"chi2_p","tt",0.6)), include.rownames=T)

########################################
# plots based on tabdat1
########################################
plottab <- function(res,index,type,miss,woline=.8,nom,ylim=c(0,1.2)){
  
  dat <- tabdat1(res,index,type,miss)
  k0 <- 0
  ntname <- c(expression(N[t]*"="*1),expression(N[t]*"="*2),expression(N[t]*"="*3),expression(N[t]*"="*4),
              expression(N[t]*"="*5),expression(N[t]*"="*10),expression(N[t]*"="*15))
  seqname <- c(1,5,2,6,3,7,4)
  
  if(type=="tt1"){
    dat <- dat[-1,]
    k0 <- 1
    ntname <- ntname[-1]
    seqname<- c(2,5,3,6,4,7)-1
  }
  
  j<-1
  plot(log(as.numeric(colnames(dat))),dat[j,],ylim=ylim,xlim=log(c(31,2501)),col=j+k0,pch=j+k0,type="b",
       axes=F,xlab="N",ylab="% rejected",main=nom)
  
  for(j in 2:nrow(dat)){#j<-length(levels(cond0$Nt.f))
    par(new=T)
    
    plot(log(as.numeric(colnames(dat))),dat[j,],ylim=ylim,xlim=log(c(31,2501)),col=j+k0,pch=j+k0,type="b",
         axes=F,ylab="",xlab="")
    
  }
  
  axis(1,at=log(as.numeric(colnames(dat))),labels=as.numeric(colnames(dat)))
  
  axis(2,seq(0,ylim[2]-.2,length.out=5))
  abline(h=woline,lty=3)
  
  #ntname
  
  
  legend("top",
         legend=ntname[seqname],
         #legend=ntname,
         #legend=substitute(expression(N[t],"=",nn), list(nn=rownames(dat))),
         #legend=bquote(.(parse(text=paste(N[t],"~",rownames(dat),sep="")))),
         #legend=bquote(N[t]~"="~.(rownames(dat)) ),
         #expression(N[t],"=",paste(rownames(dat))),
         lty=1,col=(1:nrow(dat)+k0)[seqname],pch=(1:nrow(dat)+k0)[seqname],
         bty="n",ncol=4)#horiz = T,
}


#######################################################################
nomplot1 <- c(fitnom,"chi2/df","chi2_p")
nomplot2 <- c("CFI","TLI","RMSEA","SRMR","GFI","AGFI","CGFI",expression(chi^2/df),expression(p(chi^2)))

pdf("type1error01_v1.pdf",height=3*4,width=3*4)
par(mfrow=c(3,3))
for(j in 1:length(nomplot1)){plottab(res4true,nomplot1[j],"tt",0.3,woline=.05,nomplot2[j])}
for(j in 1:length(nomplot1)){plottab(res4true,nomplot1[j],"tt",0.6,woline=.05,nomplot2[j])}
for(j in 1:length(nomplot1)){plottab(res4true,nomplot1[j],"tt1",0.3,woline=.05,nomplot2[j])}
for(j in 1:length(nomplot1)){plottab(res4true,nomplot1[j],"tt1",0.6,woline=.05,nomplot2[j])}
dev.off()


#######################################################################
pdf("power01_v1.pdf",height=3*4,width=3*4)
par(mfrow=c(3,3))
for(j in 1:length(nomplot1)){plottab(res4false,nomplot1[j],"tt",0.3,woline=.80,nomplot2[j])}
for(j in 1:length(nomplot1)){plottab(res4false,nomplot1[j],"tt",0.6,woline=.80,nomplot2[j])}
for(j in 1:length(nomplot1)){plottab(res4false,nomplot1[j],"tt1",0.3,woline=.80,nomplot2[j])}
for(j in 1:length(nomplot1)){plottab(res4false,nomplot1[j],"tt1",0.6,woline=.80,nomplot2[j])}
dev.off()



#######################################################################
nomplot1b <- c(fitnom,"chi2/df")
nomplot2b <- c("CFI","TLI","RMSEA","SRMR","GFI","AGFI","CGFI",expression(chi^2/df))

wo1 <- c(.95,.97,.05,.05,.95,.95,.95,2)
#whichplot <- c(1:length(nomplot1b))
whichplot <- c(1,2,5,3,4,8)
ylim1 <- matrix(c(rep(0,8),1.2,1.2,0.4,0.4,1.2,1.2,1.2,3),8,2,byrow=F)


ylimi <- matrix(c(-.1,0,
                  -.1,0,
                  0,.1,
                  0,.1,
                  -.1,0,
                  -.1,0,
                  -.1,.3,
                  0,30),8,2,byrow=T)


pdf("meanmodels_v1.pdf",height=2*4,width=6*4)
par(mfrow=c(3,6))
for(j in whichplot){plottab(res3true,nomplot1b[j],"tt",0.3,woline=wo1[j],nomplot2b[j],ylim=ylim1[j,])}
for(j in whichplot){plottab(res3false,nomplot1b[j],"tt",0.3,woline=wo1[j],nomplot2b[j],ylim=ylim1[j,])}
for(j in whichplot){plotdif(res3true,res3false,nomplot1a[j],"tt",0.3,woline=c(-0.02,0,0.02),nomplot2a[j],ylim=ylimi[j,])}

for(j in whichplot){plottab(res3true,nomplot1b[j],"tt",0.6,woline=wo1[j],nomplot2b[j],ylim=ylim1[j,])}
for(j in whichplot){plottab(res3false,nomplot1b[j],"tt",0.6,woline=wo1[j],nomplot2b[j],ylim=ylim1[j,])}
for(j in whichplot){plotdif(res3true,res3false,nomplot1a[j],"tt",0.6,woline=c(-0.02,0,0.02),nomplot2a[j],ylim=ylimi[j,])}

for(j in whichplot){plottab(res3true,nomplot1b[j],"tt1",0.3,woline=wo1[j],nomplot2b[j],ylim=ylim1[j,])}
for(j in whichplot){plottab(res3false,nomplot1b[j],"tt1",0.3,woline=wo1[j],nomplot2b[j],ylim=ylim1[j,])}
for(j in whichplot){plotdif(res3true,res3false,nomplot1a[j],"tt1",0.3,woline=c(-0.02,0,0.02),nomplot2a[j],ylim=ylimi[j,])}

for(j in whichplot){plottab(res3true,nomplot1b[j],"tt1",0.6,woline=wo1[j],nomplot2b[j],ylim=ylim1[j,])}
for(j in whichplot){plottab(res3false,nomplot1b[j],"tt1",0.6,woline=wo1[j],nomplot2b[j],ylim=ylim1[j,])}
for(j in whichplot){plotdif(res3true,res3false,nomplot1a[j],"tt1",0.6,woline=c(-0.02,0,0.02),nomplot2a[j],ylim=ylimi[j,])}
dev.off()


#######################################################################
nomplot1 <- c(fitnom,"chi2/df","chi2_p")
nomplot2 <- c("CFI","TLI","RMSEA","SRMR","GFI","AGFI","CGFI",expression(chi^2/df),expression(p(chi^2)))

wo1 <- rep(.05,9)
wo2 <- rep(.80,9)
whichplot <- c(1:length(nomplot1))
#whichplot <- c(1,2,5,3,4,8)
ylim1 <- matrix(c(rep(0,9),rep(1.2,9)),9,2,byrow=F)



pdf("rejectmodels_v1.pdf",height=3*4,width=9*4)
par(mfrow=c(3,9))
for(j in whichplot){plottab(res4true,nomplot1[j],"tt",0.3,woline=wo1[j],nomplot2[j],ylim=ylim1[j,])}
for(j in whichplot){plottab(res4false,nomplot1[j],"tt",0.3,woline=wo2[j],nomplot2[j],ylim=ylim1[j,])}
for(j in whichplot){plotdif(res4true,res4false,nomplot1[j],"tt",0.3,woline=c(-0.02,0,0.02),nomplot2[j],ylim=ylim1[j,])}

for(j in whichplot){plottab(res4true,nomplot1[j],"tt",0.6,woline=wo1[j],nomplot2[j],ylim=ylim1[j,])}
for(j in whichplot){plottab(res4false,nomplot1[j],"tt",0.6,woline=wo2[j],nomplot2[j],ylim=ylim1[j,])}
for(j in whichplot){plotdif(res4true,res4false,nomplot1[j],"tt",0.6,woline=c(-0.02,0,0.02),nomplot2[j],ylim=ylim1[j,])}

for(j in whichplot){plottab(res4true,nomplot1[j],"tt1",0.3,woline=wo1[j],nomplot2[j],ylim=ylim1[j,])}
for(j in whichplot){plottab(res4false,nomplot1[j],"tt1",0.3,woline=wo2[j],nomplot2[j],ylim=ylim1[j,])}
for(j in whichplot){plotdif(res4true,res4false,nomplot1[j],"tt1",0.3,woline=c(-0.02,0,0.02),nomplot2[j],ylim=ylim1[j,])}

for(j in whichplot){plottab(res4true,nomplot1[j],"tt1",0.6,woline=wo1[j],nomplot2[j],ylim=ylim1[j,])}
for(j in whichplot){plottab(res4false,nomplot1[j],"tt1",0.6,woline=wo2[j],nomplot2[j],ylim=ylim1[j,])}
for(j in whichplot){plotdif(res4true,res4false,nomplot1[j],"tt1",0.6,woline=c(-0.02,0,0.02),nomplot2[j],ylim=ylim1[j,])}
dev.off()

ylimi <- matrix(c(0,1.2,
                  0,1.2,
                  0,1.2,
                  0,1,
                  0,1,
                  0,1,
                  0,1,
                  0,1,
                  0,1),9,2,byrow=T)
#######################################################################



par(mfrow=c(3,2))
plottab(res3,"cfi","none",0,woline=.05,"CFI")
plottab(res2,"cfi","none",0,woline=.95,"CFI")

plottab(res3,"cfi","tt",0.6,woline=.80,"CFI")
plottab(res2,"cfi","tt",0.3,woline=.95,"CFI")

plottab(res3,"cfi","tt",0.3,woline=.80,"CFI")
plottab(res2,"cfi","tt",0.6,woline=.95,"CFI")

print(xtable(tabdat1(res3,"cfi","none",0)), include.rownames=T)
print(xtable(tabdat1(res2,"cfi","none",0)), include.rownames=T)

tabdat1(res2,"cfi","none",0)-tabdat1(res2,"cfi","tt",0.3)
tabdat1(res2,"cfi","none",0)-tabdat1(res2,"cfi","tt",0.6)


plottab(res2,"chi2/df","none",0,woline=1,"CFI")
plottab(res2,"chi2/df","tt",0.6,woline=1,"CFI")


########################################
# plots based on tabdat1
########################################
plotdif <- function(restrue,resfalse,index,type,miss,woline=.8,nom,ylim=ylim){
  
  dat1 <- tabdat1(restrue,index,type,miss)
  dat2 <- tabdat1(resfalse,index,type,miss)
  dat <- dat2-dat1
  
  k0 <- 0
  ntname <- c(expression(N[t]*"="*1),expression(N[t]*"="*2),expression(N[t]*"="*3),expression(N[t]*"="*4),
              expression(N[t]*"="*5),expression(N[t]*"="*10),expression(N[t]*"="*15))
  seqname <- c(1,5,2,6,3,7,4)
  
  if(type=="tt1"){
    dat <- dat[-1,]
    k0 <- 1
    ntname <- ntname[-1]
    seqname<- c(2,5,3,6,4,7)-1
  }
  
  j<-1
  plot(log(as.numeric(colnames(dat))),dat[j,],ylim=ylim,xlim=log(c(31,2501)),col=j+k0,pch=j+k0,type="b",
       axes=F,xlab="N",ylab="Difference",main=nom)
  
  for(j in 2:nrow(dat)){#j<-length(levels(cond0$Nt.f))
    par(new=T)
    
    plot(log(as.numeric(colnames(dat))),dat[j,],ylim=ylim,xlim=log(c(31,2501)),col=j+k0,pch=j+k0,type="b",
         axes=F,ylab="",xlab="")
    
  }
  
  axis(1,at=log(as.numeric(colnames(dat))),labels=as.numeric(colnames(dat)))
  
  axis(2,c(ylim,0))
  abline(h=woline,lty=3)
  
  #ntname
  
  
  legend("top",
         legend=ntname[seqname],
         #legend=ntname,
         #legend=substitute(expression(N[t],"=",nn), list(nn=rownames(dat))),
         #legend=bquote(.(parse(text=paste(N[t],"~",rownames(dat),sep="")))),
         #legend=bquote(N[t]~"="~.(rownames(dat)) ),
         #expression(N[t],"=",paste(rownames(dat))),
         lty=1,col=(1:nrow(dat)+k0)[seqname],pch=(1:nrow(dat)+k0)[seqname],
         bty="n",ncol=4)#horiz = T,
}

par(mfrow=c(2,5))

nomplot1a <- c(fitnom,"chi2/df")
nomplot2a <- c(expression(Delta~"CFI"),expression(Delta~"TLI"),expression(Delta~"RMSEA"),expression(Delta~"SRMR"),expression(Delta~"GFI"),
              expression(Delta~"AGFI"),expression(Delta~"CGFI"),expression(Delta~chi^2/df))

pdf("meandiff_v1.pdf",height=2*4,width=4*4)
par(mfrow=c(2,4))

ylimi <- matrix(c(-.1,0,
                  -.1,0,
                  0,.1,
                  0,.1,
                  -.1,0,
                  -.1,0,
                  -.1,.3,
                  0,30),8,2,byrow=T)

for(j in 1:length(nomplot1a)){plotdif(res3true,res3false,nomplot1a[j],"tt",0.3,woline=c(-0.02,0,0.02),nomplot2a[j],ylim=ylimi[j,])}
for(j in 1:length(nomplot1a)){plotdif(res3true,res3false,nomplot1a[j],"tt",0.6,woline=c(-0.02,0,0.02),nomplot2a[j],ylim=ylimi[j,])}
for(j in 1:length(nomplot1a)){plotdif(res3true,res3false,nomplot1a[j],"tt1",0.3,woline=c(-0.02,0,0.02),nomplot2a[j],ylim=ylimi[j,])}
for(j in 1:length(nomplot1a)){plotdif(res3true,res3false,nomplot1a[j],"tt1",0.6,woline=c(-0.02,0,0.02),nomplot2a[j],ylim=ylimi[j,])}
dev.off()


##########################################################################################
##########################################################################################
plottab3 <- function(restrue,resfalse,index,type,miss,woline=.8,nom,ylim=c(0,1.2)){
  #restrue <- res4true
  #resfalse <- res4true
  #type <- "tt";miss<-0.3;index="cfi";nom<-"CFI"
  dat1 <- tabdat1(restrue,index,type,miss)
  dat2 <- tabdat1(resfalse,index,type,miss)
  dat <- dat2
  
  datgood <- matrix(0,nrow(dat1),ncol(dat1))
  for(k in 1:nrow(dat2)){
    for(j in 1:ncol(dat2)){
      if(is.na(dat1[k,j])==F){
        if(dat1[k,j]<.10){
          datgood[k,j] <- 0
        }else{
          datgood[k,j] <- 1
        }
      }
    }
  }
  
  
  k0 <- 0
  ntname <- c(expression(N[t]*"="*1),expression(N[t]*"="*2),expression(N[t]*"="*3),expression(N[t]*"="*4),
              expression(N[t]*"="*5),expression(N[t]*"="*10),expression(N[t]*"="*15))
  seqname <- c(1,5,2,6,3,7,4)
  
  if(type=="tt1"){
    dat <- dat[-1,]
    k0 <- 1
    ntname <- ntname[-1]
    seqname<- c(2,5,3,6,4,7)-1
  }
  
  j<-1
  par(new=F)
  for(k in 1:(ncol(dat)-1)){
    plot(log(as.numeric(colnames(dat)))[k:(k+1)],dat[j,k:(k+1)],ylim=ylim,xlim=log(c(31,2501)),col=j+k0,pch=j+k0,type="b",lty=datgood[j,k]+1,
         axes=F,xlab="N",ylab="Difference",main=nom)
    par(new=T)
  }
  
  for(j in 2:nrow(dat)){#j<-length(levels(cond0$Nt.f))
    par(new=T)
    
    for(k in 1:(ncol(dat)-1)){
      plot(log(as.numeric(colnames(dat)))[k:(k+1)],dat[j,k:(k+1)],ylim=ylim,xlim=log(c(31,2501)),col=j+k0,pch=j+k0,type="b",lty=datgood[j,k]+1,
         axes=F,ylab="",xlab="")
      par(new=T)
    }
  }
  
  axis(1,at=log(as.numeric(colnames(dat))),labels=as.numeric(colnames(dat)))
  
  axis(2,c(ylim,0))
  abline(h=woline,lty=3)
  
  #ntname
  
  
  legend("top",
         legend=ntname[seqname],
         #legend=ntname,
         #legend=substitute(expression(N[t],"=",nn), list(nn=rownames(dat))),
         #legend=bquote(.(parse(text=paste(N[t],"~",rownames(dat),sep="")))),
         #legend=bquote(N[t]~"="~.(rownames(dat)) ),
         #expression(N[t],"=",paste(rownames(dat))),
         lty=1,col=(1:nrow(dat)+k0)[seqname],pch=(1:nrow(dat)+k0)[seqname],
         bty="n",ncol=4)#horiz = T,
}


pdf("powermodels_v2.pdf",height=4*4,width=9*4)
par(mfrow=c(4,9))
for(j in whichplot){plottab3(res4true,res4false,nomplot1[j],"tt",0.3,woline=wo2[j],nomplot2[j],ylim=ylim1[j,])}
for(j in whichplot){plottab3(res4true,res4false,nomplot1[j],"tt",0.6,woline=wo2[j],nomplot2[j],ylim=ylim1[j,])}
for(j in whichplot){plottab3(res4true,res4false,nomplot1[j],"tt1",0.3,woline=wo2[j],nomplot2[j],ylim=ylim1[j,])}
for(j in whichplot){plottab3(res4true,res4false,nomplot1[j],"tt1",0.6,woline=wo2[j],nomplot2[j],ylim=ylim1[j,])}
dev.off()







