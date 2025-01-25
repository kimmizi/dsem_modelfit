# This file does the following:
# - summarizes information about simulation
# - reads the different csv files from simulation
# - cleans data, calculates cut-off values and stacks them into one big data frame
# - plots data



current_dir = "/Users/kimzierahn/PycharmProjects/dsem_modelfit/exp/2025-01-25"
setwd(current_dir)



################################################################################
# Information about the simulation 
################################################################################

#N_p = c(31, 91, 121, 151, 181, 211, 501, 1001, 1501, 2001, 2501) 
N_p <- c(31, 61, 91, 121, 211, 501, 1001, 1501)
#N_p = c(31, 61, 91) 

#N_t = 15
N_t = 1:5

Size_crossloading = c(0, .3, .6)

Type_crossloading = c("none", "tt", "tt1")

N_sim_samples = 10 

# frequentist fit indices
fitnom_lavaan = c("npar","fmin","chisq","df","pvalue","baseline.chisq","baseline.df",
                   "baseline.pvalue","cfi","tli","nnfi","rfi","nfi","pnfi","ifi","rni",
                   "logl","unrestricted.logl","aic","bic","ntotal","bic2","rmsea","rmsea.ci.lower",
                   "rmsea.ci.upper","rmsea.ci.level","rmsea.pvalue","rmsea.close.h0","rmsea.notclose.pvalue",
                   "rmsea.notclose.h0","rmr","rmr_nomean","srmr","srmr_bentler","srmr_bentler_nomean",
                   "crmr","crmr_nomean","srmr_mplus","srmr_mplus_nomean","cn_05","cn_01","gfi",
                   "agfi","pgfi","mfi","ecvi") 

fitnom_lavaan_red = c(
  "pvalue_rejection_rate",               # Percentage of p-values <= 0.05 (rejection rate)
  "chisq_df_rejection_rate",             # Percentage of chi-sq > 2 * df
  "cfi_rejection_rate",                  # Percentage of CFI <= 0.95
  "tli_rejection_rate",                  # Percentage of TLI <= 0.97
  "rmsea_rejection_rate",                # Percentage of RMSEA >= 0.05
  "srmr_rejection_rate",                 # Percentage of SRMR >= 0.05
  "gfi_rejection_rate",                  # Percentage of GFI <= 0.95
  "agfi_rejection_rate",                 # Percentage of AGFI <= 0.95
  "cgfi_rejection_rate"                 # Percentage of CGFI <= 0.95
)

# bayesian fit indices
fitnom_blavaan = c("npar", "PPP", "MargLogLik", 
                    "BRMSEA", "BGammaHat", "adjBGammaHat", "BMc", 
                    "BCFI_0C", "BTLI_0C", "BNFI_0C",
                    "BCFI_0A", "BTLI_0A", "BNFI_0A",
                    "chisq", "pd")



################################################################################
# Reading data 
################################################################################

# define function that reads csv file
getdata_lavaan <- function(file_name){
  
  data <- try(read.csv(file_name), silent = TRUE)

  if(!inherits(data, "try-error")){
    data
  }else{
    matrix(NA, N_sim_samples, length(fitnom_lavaan))    
  }
}

getdata_blavaan <- function(file_name){
  
  data <- try(read.csv(file_name), silent = TRUE)
  
  if(!inherits(data, "try-error")){
    data
  }else{
    matrix(NA, N_sim_samples, length(fitnom_blavaan))    
  }
}




################################################################################
# Cleaning data 
################################################################################

# Initializing final dataframes
condition = data.frame(matrix(data = NA, nrow = 100, ncol = 4))
colnames(condition) = c("n_t", "n_p", "Type_misfit", "Size_misfit")

data_all_lav = data.frame(matrix(NA, 100, length(fitnom_lavaan_red)))
colnames(data_all_lav) = fitnom_lavaan_red

index = 1

### Combining data frames to one big data frame & calculate cut-offs
# define loop that loops over number of persons, timepoints, type of missfit and size of missfit
for (i_pers in seq_along(N_p)) {
  for (i_timep in seq_along(N_t)) {
    
    n_p = N_p[i_pers]
    n_t = N_t[i_timep]
    cat("Person size: ", n_p, "\n")
    cat("Time points: ", n_t, "\n")
    
    for (Type_misfit in Type_crossloading) {
      for (Size_misfit in Size_crossloading) {
        if(Size_misfit == 0 && Type_misfit != "none"){
        }else if(Size_misfit != 0 && Type_misfit == "none"){
        }else{ 
          
          # Initialize empty matrix that will store data
          data_lav_1 = data.frame(matrix(data = NA, nrow = 1, ncol = length(fitnom_lavaan)))
          colnames(data_lav_1) = fitnom_lavaan
          
          #data_blav_1 = data.frame(matrix(data = NA, nrow = 1, ncol = length(fitnom_blavaan)))
          #colnames(data_blav_1) = fitnom_lavaan
          
          
          # Get name of csv file & read data
          # lavaan
          file_name_lav = paste("dsem", n_p, n_t, N_sim_samples, Type_misfit,
                            Size_misfit, "lav_.csv", sep = "_")
          data_temp_lav = getdata_lavaan(file_name_lav)
          
          data_lav = rbind(data_lav_1, data_temp_lav)
          data_lav = data.frame(data_lav[-1,])
          
          
          # blavaan
          #file_name_blav = paste("dsem", n_p, n_t, N_sim_samples, Type_misfit,
          #                      Size_misfit, "blav_.csv", sep = "_")
          #data_temp_blav = getdata_blavaan(file_name_blav)
          
          #data_blav = rbind(data_blav_1, data_temp_blav)
          #data_blav = data.frame(data_blav[-1,])
          
          
          
          #### FIT INDICES ####
          data_lav$cgfi = data_lav$gfi + (n_t*6 + 1) * n_t*6 / data_lav$npar / n_p
          #data_blav$cgfi = data_blav$gfi + (n_t*6 + 1) * n_t*6 / data_blav$npar / n_p
          
          
          # FIT INDICES
          #data_all_lav[index] = mean(data_lav$chisq/data_lav$df, na.rm = T) # CHI2
          #data_all_lav[index, 1:length(fitnom_lavaan) + 1] = apply(data_lav[fitnom_lavaan], 2, mean, na.rm = T)
          #data_all_lav[index, length(fitnom_lavaan) + 2] = length(na.omit(data_lav$npar)) # LENGHT OF MISSING DATA
          
          # CUT OFFS
          # cutoffs good fit (misfit): percentage of rejection
          
          # lavaan
          data_all_lav[index, 1] = 1 - mean(data_lav$pvalue > .05, na.rm = T) # if p-value is larger than 5% then its 1 = percentage of rejection = POWER
          data_all_lav[index, 2] = 1 - mean(data_lav$chisq < 2*data_lav$df, na.rm = T)
          data_all_lav[index, 3] = 1 - mean(data_lav$cfi > .95, na.rm = T)
          data_all_lav[index, 4] = 1 - mean(data_lav$tli > .97, na.rm = T)
          data_all_lav[index, 5] = 1 - mean(data_lav$rmsea < .05, na.rm = T)
          data_all_lav[index, 6] = 1 - mean(data_lav$srmr < .05, na.rm = T)
          data_all_lav[index, 7] = 1 - mean(data_lav$gfi > .95, na.rm = T)
          data_all_lav[index, 8] = 1 - mean(data_lav$agfi > .95, na.rm = T)
          data_all_lav[index, 9] = 1 - mean(data_lav$cgfi > .95, na.rm = T)
          
          # blavaan
          #data_all_blav[index, 1] = 1 - mean(data_blav$BRMSEA > .05, na.rm = T)
          #data_all_blav[index, 2] = 1 - mean(resb1$BGammaHat < .95, na.rm = T)
          #data_all_blav[index, 3] = 1 - mean(resb1$adjBGammaHat < .95, na.rm = T)
          #data_all_blav[index, 4] = 1 - mean(resb1$BMc < .95, na.rm = T)
          #data_all_blav[index, 5] = 1 - mean(resb1$BCFI_0C < .95, na.rm = T)
          #data_all_blav[index, 6] = 1 - mean(resb1$BTLI_0C < .95, na.rm = T)
          #data_all_blav[index, 7] = 1 - mean(resb1$BNFI_0C < .95, na.rm = T)
          #data_all_blav[index, 8] = 1 - mean(resb1$BCFI_0A < .95, na.rm = T)
          #data_all_blav[index, 9] = 1 - mean(resb1$BTLI_0A < .95, na.rm = T)
          #data_all_blav[index, 10] = 1 - mean(resb1$BNFI_0A < .95, na.rm = T)


          condition[index, ] = c(n_t, n_p, Type_misfit, Size_misfit)
          index = index + 1 
          
        }      
      }
    }
  }
}


# Combine data and omit missing data
data_no_na_lav = na.omit(cbind(condition, data_all_lav))
#data_no_na_blav = na.omit(cbind(condition, data_all_blav))


# Convert columns into numeric and factors for plotting
# as continous measure:
data_no_na_lav$n_t = as.numeric(data_no_na_lav$n_t)
data_no_na_lav$n_p = as.numeric(data_no_na_lav$n_p)
data_no_na_lav$Size_misfit = as.numeric(data_no_na_lav$Size_misfit)

#data_no_na_blav$n_t = as.numeric(data_no_na_blav$n_t)
#data_no_na_blav$n_p = as.numeric(data_no_na_blav$n_p)
#data_no_na_blav$Size_misfit = as.numeric(data_no_na_blav$Size_misfit)


# as factor
#data_no_na_lav$n_t_f = factor(data_no_na_lav$n_t, levels = c("1","2","3","4","5","10","15","30"))
data_no_na_lav$n_t_f = factor(data_no_na_lav$n_t, levels = c("1", "2", "3", "4", "5"))
data_no_na_lav$n_p_f = factor(data_no_na_lav$n_p, levels = paste0(sort(N_p)))
data_no_na_lav$Type_misfit = as.factor(data_no_na_lav$Type_misfit)

#data_no_na_blav$n_t_f = factor(data_no_na_blav$n_t, levels = c("1","2","3","4","5","10","15","30"))
#data_no_na_blav$n_t_f = factor(data_no_na_blav$n_t, levels = c("1","2","3"))
#data_no_na_blav$n_p_f = factor(data_no_na_blav$n_p, levels = paste0(sort(N_p)))
#data_no_na_blav$Type_misfit = as.factor(data_no_na_blav$Type_misfit)



################################################################################
# Plotting 
################################################################################

#### Plot version 1 ####
# X-Axis: N_p (log-scaled)
# Y-Axis: % rejected
# Different lines for different number of time points
plot_results = function(df, miss_size, miss_type, fit_index, power = 0.8){
  
  # Initialize plot
  par(new = FALSE)
  
  index = 1
  
  condition = df$Type_misfit == miss_type & df$Size_misfit == miss_size & df$n_t_f == levels(df$n_t_f)[index]
  
  plot(sort(log(df$n_p[condition])), df[condition, fit_index][rank(df$n_p[condition])], ylim = c(0, 1), 
       xlim = c(log(min(df$n_p)), log(max(df$n_p))), col = index, pch = index, type = "b",
       axes = F, xlab = "N", ylab = "% rejected", main = fit_index)
  
  for (index in 2:length(levels(df$n_t_f))) {
    par(new = TRUE)
    condition = df$Type_misfit == miss_type & df$Size_misfit == miss_size & df$n_t_f == levels(df$n_t_f)[index]
    
    plot(sort(log(df$n_p[condition])), df[condition, fit_index][rank(df$n_p[condition])], ylim = c(0, 1), 
        xlim = c(log(min(df$n_p)), log(max(df$n_p))), col = index, pch = index, type = "b",
        axes = F, ylab = "", xlab = "")
  }
  
  axis(1, at = log(N_p), labels = N_p)
  axis(2)
  abline(h = power, lty = 3)
  
  legend("topright", paste0("N_t = ", levels(df$n_t_f)), lty = 1, col = 1:length(levels(df$n_t_f)), pch = 1:length(levels(df$n_t_f)))
}


pdf("type1error01_v1_test.pdf", height = 3*4, width = 2*4) # create PDF
par(mfrow = c(3, 2))

# LABELING:
plot_results(data_no_na_lav, miss_size = 0, miss_type = "none", fit_index = "pvalue_rejection_rate", power = .05)
plot_results(data_no_na_lav, miss_size = 0, miss_type = "none", fit_index = "chisq_df_rejection_rate", power = .05)
plot_results(data_no_na_lav, miss_size = 0, miss_type = "none", fit_index = "rmsea_rejection_rate", power = .05)
plot_results(data_no_na_lav, miss_size = 0, miss_type = "none", fit_index = "srmr_rejection_rate", power = .05)
plot_results(data_no_na_lav, miss_size = 0, miss_type = "none", fit_index = "cfi_rejection_rate", power = .05)
plot_results(data_no_na_lav, miss_size = 0, miss_type = "none", fit_index = "tli_rejection_rate", power = .05)
dev.off()




#### Plot version 2 ####
# Plots both power and type 1
# X-Axis: N_p (log-scaled)
# Y-Axis: % rejected
# Different lines for different number of time points

plot_results_pow_typ1 = function(df, miss_size, miss_type, fit_index, timepoint, power = 0.8, type1 = 0.05){
  
  # Initialize plot
  par(new = FALSE)
  
  condition_none = df$Type_misfit == "none" & df$Size_misfit == 0 & df$n_t_f == timepoint
  condition_missfit = df$Type_misfit == miss_type & df$Size_misfit == miss_size & df$n_t_f == timepoint
  
  
  plot(sort(log(df$n_p[condition_missfit])), df[condition_missfit , fit_index][rank(df$n_p[condition_missfit])], ylim = c(0, 1), 
       xlim = log(c(min(df$n_p), max(df$n_p))), col = "black", pch = 1, type = "b",
       axes = F, xlab = "N", ylab = "% rejected", main = fit_index)
  
  par(new = T)
  plot(sort(log(df$n_p[condition_none])), df[condition_none , fit_index][rank(df$n_p[condition_none])], ylim = c(0, 1), 
       xlim = log(c(min(df$n_p), max(df$n_p))), col = "red", pch = 2, type = "b",
       axes = F, xlab = "N", ylab = "% rejected", main = fit_index)
  
  axis(1, at = log(N_p), labels = N_p)
  axis(2)
  abline(h = c(type1, power), lty = 3)
  
  legend("right", c("Power", "Type I error"), lty = 1, col = c("black", "red"), pch = 1:2, bty = "n")
  
}


pdf("fitmean_dfg_test.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2,3)) 

plot_results_pow_typ1(df = data_no_na_lav, miss_size = 0.6, miss_type = "tt1", fit_index = "pvalue_rejection_rate", timepoint = 3)
plot_results_pow_typ1(df = data_no_na_lav, miss_size = 0.6, miss_type = "tt1", fit_index = "chisq_df_rejection_rate", timepoint = 3)
plot_results_pow_typ1(df = data_no_na_lav, miss_size = 0.6, miss_type = "tt1", fit_index = "cfi_rejection_rate", timepoint = 3)
plot_results_pow_typ1(df = data_no_na_lav, miss_size = 0.6, miss_type = "tt1", fit_index = "tli_rejection_rate", timepoint = 3)
plot_results_pow_typ1(df = data_no_na_lav, miss_size = 0.6, miss_type = "tt1", fit_index = "rmsea_rejection_rate", timepoint = 3)
#plot_results_pow_typ1(df = data_no_na_lav, miss_size = 0.6, miss_type = "tt1", fit_index = "gfi_rejection_rate", timepoint = 3)
plot_results_pow_typ1(df = data_no_na_lav, miss_size = 0.6, miss_type = "tt1", fit_index = "srmr_rejection_rate", timepoint = 3)
#plot_results_pow_typ1(df = data_no_na_lav, miss_size = 0.6, miss_type = "tt1", fit_index = "cgfi_rejection_rate", timepoint = 3)

dev.off()



#### Plot version 3 ####
#### Plotting lavaan vs. blavaan ####

plot_results_pow_typ1 = function(df_lav, df_blav, miss_size, miss_type, fit_index, timepoint, power = 0.8, type1 = 0.05){
  
  # Initialize plot
  par(new = FALSE)
  
  # lavaan
  condition_none = df_lav$Type_misfit == "none" & df_lav$Size_misfit == 0 & df_lav$n_t_f == timepoint
  condition_missfit = df_lav$Type_misfit == miss_type & df_lav$Size_misfit == miss_size & df_lav$n_t_f == timepoints
  
  plot(sort(log(df_lav$n_p[condition_missfit])), df_lav[condition_missfit , fit_index][rank(df_lav$n_p[condition_missfit])], ylim = c(0, 1), 
       xlim = log(c(min(df_lav$n_p), max(df_lav$n_p))), col = "black", pch = 1, type = "b",
       axes = F, xlab = "N", ylab = "% rejected", main = fit_index)
  
  par(new = T)
  plot(sort(log(df_lav$n_p[condition_none])), df_lav[condition_none , fit_index][rank(df_lav$n_p[condition_none])], ylim = c(0, 1), 
       xlim = log(c(min(df_lav$n_p), max(df_lav$n_p))), col = "red", pch = 2, type = "b",
       axes = F, xlab = "", ylab = "", main = fit_index)
  
  
  # blavaan
  condition_none = df_blav$Type_misfit == "none" & df_blav$Size_misfit == 0 & df_blav$n_t_f == timepoint
  condition_missfit = df_blav$Type_misfit == miss_type & df_blav$Size_misfit == miss_size & df_blav$n_t_f == timepoints
  
  par(new = T)
  plot(sort(log(df_blav$n_p[condition_missfit])), df_blav[condition_missfit , fit_index][rank(df_blav$n_p[condition_missfit])], ylim = c(0, 1), 
       xlim = log(c(min(df_blav$n_p), max(df_blav$n_p))), col = "blue", pch = 3, type = "b",
       axes = F, xlab = "", ylab = "", main = fit_index)
  
  par(new = T)
  plot(sort(log(df_blav$n_p[condition_none])), df_blav[condition_none , fit_index][rank(df_blav$n_p[condition_none])], ylim = c(0, 1), 
       xlim = log(c(min(df_blav$n_p), max(df_blav$n_p))), col = "green", pch = 4, type = "b",
       axes = F, xlab = "", ylab = "", main = fit_index)
  
  
  
  
  axis(1, log(person_size_SIMULATE), person_size_SIMULATE)
  axis(2)
  abline(h = c(type1, power), lty = 3)
  
  legend("right", 
         c("Lavaan (Power)", "Lavaan (Type I)", "Blavaan (Power)", "Blavaan (Type I)"), 
         lty = 1, col = c("black", "red", "blue", "green"), pch = 1:4, bty = "n")  
}


