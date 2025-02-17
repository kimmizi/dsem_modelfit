# This file does the following:
# - summarizes information about simulation
# - reads the different csv files from simulation
# - cleans data, calculates cut-off values and stacks them into one big data frame
# - plots data


library(ggplot2)
library(dplyr)
library(tidyr)


# # Kim's dir
# #current_dir = "/Users/kimzierahn/PycharmProjects/dsem_modelfit/exp/2025-01-25-final"
# #current_dir = "/Users/kimzierahn/PycharmProjects/dsem_modelfit/exp/2025-01-27"
# current_dir = "~/Desktop/runs_until_28_morning.nosync"

# #Mihai's dir
current_dir = "C:/Users/mihai/Documents/Faculta/Research Project/dsem_modelfit/exp"
folders = c("2025-01-26", "2025-01-27", "2025-01-28")

setwd(current_dir)

################################################################################
# Information about the simulation 
################################################################################

N_p = c(31, 61, 91, 121, 151, 181, 211, 501, 1001, 1501, 2001) 

N_t = c(1:5, 10, 15)

Size_crossloading = c(0, .3, .6)

Type_crossloading = c("none", "tt", "tt1")

#N_sim_samples = 50 
N_sim_samples = 550:700


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
  "cgfi_rejection_rate"                  # Percentage of CGFI <= 0.95
)

# bayesian fit indices
fitnom_blavaan = c("npar", "PPP", "MargLogLik", 
                    "BRMSEA", "BGammaHat", "adjBGammaHat", "BMc", 
                    "BCFI_0C", "BTLI_0C", "BNFI_0C",
                    "BCFI_0A", "BTLI_0A", "BNFI_0A",
                    "chisq", "pd")

fitnom_blavaan_red = c(
  "cfi_rejection_rate_0A",               # Percentage of CFI <= 0.95
  "cfi_rejection_rate_0C",               # Percentage of CFI <= 0.95
  "tli_rejection_rate_0A",               # Percentage of TLI <= 0.97
  "tli_rejection_rate_0C",               # Percentage of TLI <= 0.97
  "nfi_rejection_rate_0A",               # Percentage of NFI <= 0.97
  "nfi_rejection_rate_0C",               # Percentage of NFI <= 0.97
  "brmsea_rejection_rate",               # Percentage of BRMSEA >= 0.05
  "bgammahat_rejection_rate",            # Percentage of BGammaHat <= 0.95
  "adjbgammahat_rejection_rate",         # Percentage of adjBGammaHat <= 0.95
  "bmc_rejection_rate"                  # Percentage of BMc <= 0.95
)



################################################################################
# Reading data 
################################################################################

# Initializing dataframes
condition = data.frame(matrix(data = NA, 1, ncol = 6))
colnames(condition) = c("Sim_num", "n_t", "n_p", "Type_misfit", "Size_misfit", "Model")

condition2 = data.frame(matrix(data = NA, nrow = 5*length(N_p)*length(N_t), ncol = 4))
colnames(condition2) = c("n_t", "n_p", "Type_misfit", "Size_misfit")

# data_all_lav = data.frame(matrix(NA, 5*length(N_p)*length(N_t), length(fitnom_lavaan_red)))
# colnames(data_all_lav) = fitnom_lavaan_red
# 
# data_all_blav = data.frame(matrix(NA, 5*length(N_p)*length(N_t), length(fitnom_blavaan_red)))
# colnames(data_all_blav) = fitnom_blavaan_red

data_lav = data.frame(matrix(data = NA, nrow = 1, ncol = length(fitnom_lavaan)))
colnames(data_lav) = fitnom_lavaan

data_blav = data.frame(matrix(data = NA, nrow = 1, ncol = length(fitnom_blavaan)))
colnames(data_blav) = fitnom_blavaan

index = 1

for(folder in folders){
  
  directory = paste0(current_dir, "/", folder)
  setwd(directory)
  
  file_names = list.files(path = directory, full.names = TRUE)
  
  extracted_data <- data.frame(
    Sim_num = as.numeric(sub(".*dsem_(\\d+)_.*", "\\1", file_names)),  # Simulation number
    n_t = as.numeric(sub(".*_\\d+_(\\d+)_.*", "\\1", file_names)),  # Sample size
    n_p = as.numeric(sub(".*_(\\d+)_.*_(none|tt|tt1)_.*", "\\1", file_names)),  # Time point
    Type_misfit = sub(".*_(none|tt|tt1)_.*", "\\1", file_names),  # Condition
    Size_misfit = as.numeric(sub(".*_(\\d+\\.?\\d*)_(blav|lav)_.*\\.csv", "\\1", file_names)),  # Misfit size
    Model = sub(".*_(blav|lav)_.*\\.csv", "\\1", file_names)  # Model type (blav or lav)
  )
  
  for (file in file_names) {
    data_temp = read.csv(file) 
    
    if (length(data_temp) == length(fitnom_lavaan)){
      colnames(data_temp) = fitnom_lavaan
      data_lav = rbind(data_lav, data_temp)
    }
    if (length(data_temp) == length(fitnom_blavaan)){
      colnames(data_temp) = fitnom_blavaan
      data_blav = rbind(data_blav, data_temp)
    }
  }

  condition = rbind(condition, extracted_data)
}


# remove first row of NAs
data_lav = data_lav[-1, ]
data_blav = data_blav[-1, ]

condition = condition[-1, ]

# check dimensions
nrow(data_blav)
nrow(data_lav)
nrow(data_blav) + nrow(data_lav)

nrow(condition)


# filter for lavaan & blavaan
condition_lav = condition[condition$Model == "lav", ]
condition_blav = condition[condition$Model == "blav", ]

nrow(condition_lav)
nrow(condition_blav)

data_lav = cbind(condition_lav, data_lav)
data_blav = cbind(condition_blav, data_blav)


################################################################################
# Calculating Cut Offs 
################################################################################

#### Lavaan ####
data_lav$cgfi = data_lav$gfi + (data_lav$n_t * 6 + 1) * data_lav$n_t * 6 / data_lav$npar / data_lav$n_p

data_lav_fit = data_lav %>%
  group_by(n_t, n_p, Type_misfit, Size_misfit) %>%
  summarize(
    power = 1 - mean(pvalue > 0.05, na.rm = TRUE),
    chi_squared_power = 1 - mean(chisq < 2 * df, na.rm = TRUE),
    cfi_power = 1 - mean(cfi > 0.95, na.rm = TRUE),
    tli_power = 1 - mean(tli > 0.97, na.rm = TRUE),
    rmsea_power = 1 - mean(rmsea < 0.05, na.rm = TRUE),
    srmr_power = 1 - mean(srmr < 0.05, na.rm = TRUE),
    gfi_power = 1 - mean(gfi > 0.95, na.rm = TRUE),
    agfi_power = 1 - mean(agfi > 0.95, na.rm = TRUE),
    cgfi_power = 1 - mean(cgfi > 0.95, na.rm = TRUE)
  ) %>%
  ungroup()

data_lav_fit = as.data.frame(data_lav_fit)
col_names = append(c("n_t", "n_p", "Type_misfit", "Size_misfit"), fitnom_lavaan_red)
colnames(data_lav_fit) = col_names


#### Blavaan ####
data_blav_fit = data_blav %>%
  group_by(n_t, n_p, Type_misfit, Size_misfit) %>%
  summarize(
    power_bcf1_0A = 1 - mean(BCFI_0A > 0.95, na.rm = TRUE),
    power_bcf1_0C = 1 - mean(BCFI_0C > 0.95, na.rm = TRUE),
    power_btli_0A = 1 - mean(BTLI_0A > 0.97, na.rm = TRUE),
    power_btli_0C = 1 - mean(BTLI_0C > 0.97, na.rm = TRUE),
    power_bnfi_0A = 1 - mean(BNFI_0A > 0.97, na.rm = TRUE),
    power_bnfi_0C = 1 - mean(BNFI_0C > 0.97, na.rm = TRUE),
    power_brmsea = 1 - mean(BRMSEA < 0.05, na.rm = TRUE),
    power_bgamma = 1 - mean(BGammaHat > 0.95, na.rm = TRUE),
    power_adj_bgamma = 1 - mean(adjBGammaHat > 0.95, na.rm = TRUE),
    power_bmc = 1 - mean(BMc > 0.95, na.rm = TRUE)
  ) %>%
  ungroup()

data_blav_fit = as.data.frame(data_blav_fit)
col_names = append(c("n_t", "n_p", "Type_misfit", "Size_misfit"), fitnom_blavaan_red)
colnames(data_blav_fit) = col_names

data_lav_fit
data_blav_fit



################################################################################
# Descriptives 
################################################################################


#### Number NAs ####
# na_count_lav
# #na_count_blav
# 
# number_entries_lav = row_count_lav*length(data_temp_lav)
# #number_entries_blav = row_count_blav*length(data_temp_blav)
# 
# perc_na_lav = (na_count_lav/number_entries_lav)*100
# cat(round(perc_na_lav, 2), "%")
# 
# #perc_na_blav = (na_count_blav/row_count_blav)*100
# cat(round(perc_na_blav, 2), "%")


# Convert columns into numeric and factors for plotting
# as continous measure:
data_lav_fit$n_t = as.numeric(data_lav_fit$n_t)
data_lav_fit$n_p = as.numeric(data_lav_fit$n_p)
data_lav_fit$Size_misfit = as.numeric(data_lav_fit$Size_misfit)

data_blav_fit$n_t = as.numeric(data_blav_fit$n_t)
data_blav_fit$n_p = as.numeric(data_blav_fit$n_p)
data_blav_fit$Size_misfit = as.numeric(data_blav_fit$Size_misfit)


# as factor
data_lav_fit$n_t_f = factor(data_lav_fit$n_t, levels = c("1","2","3","4","5","10","15"))
#data_lav_fit$n_t_f = factor(data_lav_fit$n_t, levels = c("1", "2", "3", "4", "5"))
# data_lav_fit$n_t_f = factor(data_lav_fit$n_t, levels = c("1", "2"))
data_lav_fit$n_p_f = factor(data_lav_fit$n_p, levels = paste0(sort(N_p)))
data_lav_fit$Type_misfit = as.factor(data_lav_fit$Type_misfit)

data_blav_fit$n_t_f = factor(data_blav_fit$n_t, levels = c("1", "2", "3", "4", "5", "10", "15"))
# data_blav_fit$n_t_f = factor(data_blav_fit$n_t, levels = c("1", "2"))
data_blav_fit$n_p_f = factor(data_blav_fit$n_p, levels = paste0(sort(N_p)))
data_blav_fit$Type_misfit = as.factor(data_blav_fit$Type_misfit)



################################################################################
# Plotting 
################################################################################

#### Plot version 1 ####
# X-Axis: N_p (log-scaled)
# Y-Axis: % rejected
# Different lines for different number of time points
plot_results = function(df, miss_size, miss_type, fit_index, power = 0.8, title){
  
  # Initialize plot
  par(new = FALSE)
  
  if(miss_type == "tt1"){
    index = 2
  }else{
    index = 1
  }
  
  condition = df$Type_misfit == miss_type & df$Size_misfit == miss_size & df$n_t_f == levels(df$n_t_f)[index] & !is.na(df[[fit_index]])
  
  plot(sort(log(df$n_p[condition])), df[condition, fit_index][rank(df$n_p[condition])], ylim = c(0, 1), 
      xlim = c(log(min(df$n_p)), log(max(df$n_p))), col = index, pch = index, type = "b",
      axes = F, xlab = "N", ylab = "% rejected", main = paste(title, miss_type, miss_size))
  
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
  if(miss_type == "tt1"){
    legend("topright", paste0("N_t = ", levels(df$n_t_f)[-1]), lty = 1, col = 2:length(levels(df$n_t_f)), pch = 2:length(levels(df$n_t_f)), border = "white")
  }else{
    legend("topright", paste0("N_t = ", levels(df$n_t_f)), lty = 1, col = 1:length(levels(df$n_t_f)), pch = 1:length(levels(df$n_t_f)), border = "white")
  }
}

data_lav_fit$Type_misfit

#### LAVAAN ####
pdf("type1error01_lav.pdf", height = 3*4, width = 2*4) # create PDF
par(mfrow = c(3, 2))
plot_results(data_lav_fit, miss_size = 0, miss_type = "none", fit_index = "pvalue_rejection_rate", power = .05, title = "Type 1 Error: p-value")
plot_results(data_lav_fit, miss_size = 0, miss_type = "none", fit_index = "chisq_df_rejection_rate", power = .05, title = "Type 1 Error: chi-squared")
plot_results(data_lav_fit, miss_size = 0, miss_type = "none", fit_index = "rmsea_rejection_rate", power = .05, title = "Type 1 Error: RMSEA")
plot_results(data_lav_fit, miss_size = 0, miss_type = "none", fit_index = "srmr_rejection_rate", power = .05, title = "Type 1 Error: SRMR")
plot_results(data_lav_fit, miss_size = 0, miss_type = "none", fit_index = "cfi_rejection_rate", power = .05, title = "Type 1 Error: CFI")
plot_results(data_lav_fit, miss_size = 0, miss_type = "none", fit_index = "tli_rejection_rate", power = .05, title = "Type 1 Error: TLI")
dev.off()

pdf("power_tt1_0.6_lav.pdf", height = 3*4, width = 2*4) # create PDF
par(mfrow = c(3, 2))
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "pvalue_rejection_rate", power = .80, title = "Power: p-value")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "chisq_df_rejection_rate", power = .80, title = "Power: chi-squared")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "rmsea_rejection_rate", power = .80, title = "Power: RMSEA")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "srmr_rejection_rate", power = .80, title = "Power: SRMR")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "cfi_rejection_rate", power = .80, title = "Power: CFI")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "tli_rejection_rate", power = .80, title = "Power: TLI")
dev.off()

pdf("power_tt_0.6_lav.pdf", height = 3*4, width = 2*4) # create PDF
par(mfrow = c(3, 2))
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "pvalue_rejection_rate", power = .80, title = "Power: p-value")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "chisq_df_rejection_rate", power = .80, title = "Power: chi-squared")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "rmsea_rejection_rate", power = .80, title = "Power: RMSEA")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "srmr_rejection_rate", power = .80, title = "Power: SRMR")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "cfi_rejection_rate", power = .80, title = "Power: CFI")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate", power = .80, title = "Power: TLI")
dev.off()

pdf("power_tt1_0.3_lav.pdf", height = 3*4, width = 2*4) # create PDF
par(mfrow = c(3, 2))
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "pvalue_rejection_rate", power = .80, title = "Power: p-value")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "chisq_df_rejection_rate", power = .80, title = "Power: chi-squared")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "rmsea_rejection_rate", power = .80, title = "Power: RMSEA")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "srmr_rejection_rate", power = .80, title = "Power: SRMR")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "cfi_rejection_rate", power = .80, title = "Power: CFI")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "tli_rejection_rate", power = .80, title = "Power: TLI")
dev.off()

pdf("power_tt_0.3_lav.pdf", height = 3*4, width = 2*4) # create PDF
par(mfrow = c(3, 2))
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "pvalue_rejection_rate", power = .80, title = "Power: p-value")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "chisq_df_rejection_rate", power = .80, title = "Power: chi-squared")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "rmsea_rejection_rate", power = .80, title = "Power: RMSEA")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "srmr_rejection_rate", power = .80, title = "Power: SRMR")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "cfi_rejection_rate", power = .80, title = "Power: CFI")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "tli_rejection_rate", power = .80, title = "Power: TLI")
dev.off()


#### Plot version 1 ####
# X-Axis: N_p (log-scaled)
# Y-Axis: % rejected
# Different lines for different number of time points
plot_results = function(df, miss_size, miss_type, fit_index, power = 0.8, title){
  
  # Initialize plot
  par(new = FALSE)
  
  if(miss_type == "tt1"){
    index = 2
  }else{
    index = 1
  }
  
  condition = df$Type_misfit == miss_type & df$Size_misfit == miss_size & df$n_t_f == levels(df$n_t_f)[index] & !is.na(df[[fit_index]])
  
  plot(sort(log(df$n_p[condition])), df[condition, fit_index][rank(df$n_p[condition])], ylim = c(0, 1), 
       xlim = c(log(min(df$n_p)), log(max(df$n_p))), col = index, pch = index, type = "b",
       axes = F, xlab = "N", ylab = "% rejected", main = paste(title, miss_type, miss_size))
  
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
  if(miss_type == "tt1"){
    legend("topright", paste0("N_t = ", levels(df$n_t_f)[-1]), lty = 1, col = 2:length(levels(df$n_t_f)), pch = 2:length(levels(df$n_t_f)), border = "white")
  }else{
    legend("topright", paste0("N_t = ", levels(df$n_t_f)), lty = 1, col = 1:length(levels(df$n_t_f)), pch = 1:length(levels(df$n_t_f)), border = "white")
  }
}



#### LAVAAN ####
pdf("type1error01_lav.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 3))
plot_results(data_lav_fit, miss_size = 0, miss_type = "none", fit_index = "pvalue_rejection_rate", power = .05, title = "Type 1 Error: p-value")
plot_results(data_lav_fit, miss_size = 0, miss_type = "none", fit_index = "chisq_df_rejection_rate", power = .05, title = "Type 1 Error: chi-squared")
plot_results(data_lav_fit, miss_size = 0, miss_type = "none", fit_index = "rmsea_rejection_rate", power = .05, title = "Type 1 Error: RMSEA")
plot_results(data_lav_fit, miss_size = 0, miss_type = "none", fit_index = "srmr_rejection_rate", power = .05, title = "Type 1 Error: SRMR")
plot_results(data_lav_fit, miss_size = 0, miss_type = "none", fit_index = "cfi_rejection_rate", power = .05, title = "Type 1 Error: CFI")
plot_results(data_lav_fit, miss_size = 0, miss_type = "none", fit_index = "tli_rejection_rate", power = .05, title = "Type 1 Error: TLI")
dev.off()

pdf("power_tt1_0.6_lav.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 3))
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "pvalue_rejection_rate", power = .80, title = "Power: p-value")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "chisq_df_rejection_rate", power = .80, title = "Power: chi-squared")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "rmsea_rejection_rate", power = .80, title = "Power: RMSEA")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "srmr_rejection_rate", power = .80, title = "Power: SRMR")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "cfi_rejection_rate", power = .80, title = "Power: CFI")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "tli_rejection_rate", power = .80, title = "Power: TLI")
dev.off()

pdf("power_tt_0.6_lav.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 3))
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "pvalue_rejection_rate", power = .80, title = "Power: p-value")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "chisq_df_rejection_rate", power = .80, title = "Power: chi-squared")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "rmsea_rejection_rate", power = .80, title = "Power: RMSEA")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "srmr_rejection_rate", power = .80, title = "Power: SRMR")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "cfi_rejection_rate", power = .80, title = "Power: CFI")
plot_results(data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate", power = .80, title = "Power: TLI")
dev.off()

pdf("power_tt1_0.3_lav.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 3))
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "pvalue_rejection_rate", power = .80, title = "Power: p-value")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "chisq_df_rejection_rate", power = .80, title = "Power: chi-squared")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "rmsea_rejection_rate", power = .80, title = "Power: RMSEA")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "srmr_rejection_rate", power = .80, title = "Power: SRMR")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "cfi_rejection_rate", power = .80, title = "Power: CFI")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "tli_rejection_rate", power = .80, title = "Power: TLI")
dev.off()

pdf("power_tt_0.3_lav.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 3))
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "pvalue_rejection_rate", power = .80, title = "Power: p-value")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "chisq_df_rejection_rate", power = .80, title = "Power: chi-squared")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "rmsea_rejection_rate", power = .80, title = "Power: RMSEA")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "srmr_rejection_rate", power = .80, title = "Power: SRMR")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "cfi_rejection_rate", power = .80, title = "Power: CFI")
plot_results(data_lav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "tli_rejection_rate", power = .80, title = "Power: TLI")
dev.off()



#### BLAVAAN ####
pdf("type1error01_blav.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 3))
plot_results(data_blav_fit, miss_size = 0, miss_type = "none", fit_index = "cfi_rejection_rate_0A", power = .05, title = "Type 1 Error: CFI_0A")
plot_results(data_blav_fit, miss_size = 0, miss_type = "none", fit_index = "cfi_rejection_rate_0C", power = .05, title = "Type 1 Error: CFI_0C")
plot_results(data_blav_fit, miss_size = 0, miss_type = "none", fit_index = "tli_rejection_rate_0A", power = .05, title = "Type 1 Error: TLI_0A")
plot_results(data_blav_fit, miss_size = 0, miss_type = "none", fit_index = "tli_rejection_rate_0C", power = .05, title = "Type 1 Error: TLI_0C")
plot_results(data_blav_fit, miss_size = 0, miss_type = "none", fit_index = "nfi_rejection_rate_0A", power = .05, title = "Type 1 Error: NFI_0A")
plot_results(data_blav_fit, miss_size = 0, miss_type = "none", fit_index = "nfi_rejection_rate_0C", power = .05, title = "Type 1 Error: NFI_0C")
dev.off()

pdf("power_tt1_0.6_blav.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 3))
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "cfi_rejection_rate_0A", power = .80, title = "Power: CFI_0A")
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "cfi_rejection_rate_0C", power = .80, title = "Power: CFI_0C")
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "tli_rejection_rate_0A", power = .80, title = "Power: TLI_0A")
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "tli_rejection_rate_0C", power = .80, title = "Power: TLI_0C")
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "nfi_rejection_rate_0A", power = .80, title = "Power: NFI_0A")
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "nfi_rejection_rate_0C", power = .80, title = "Power: NFI_0C")
dev.off()

pdf("power_tt_0.6_blav.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 3))
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "cfi_rejection_rate_0A", power = .80, title = "Power: CFI_0A")
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "cfi_rejection_rate_0C", power = .80, title = "Power: CFI_0C")
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate_0A", power = .80, title = "Power: TLI_0A")
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate_0C", power = .80, title = "Power: TLI_0C")
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "nfi_rejection_rate_0A", power = .80, title = "Power: NFI_0A")
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "nfi_rejection_rate_0C", power = .80, title = "Power: NFI_0C")
dev.off()

pdf("power_tt1_0.3_blav.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 3))
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "cfi_rejection_rate_0A", power = .80, title = "Power: CFI_0A")
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "cfi_rejection_rate_0C", power = .80, title = "Power: CFI_0C")
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "tli_rejection_rate_0A", power = .80, title = "Power: TLI_0A")
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "tli_rejection_rate_0C", power = .80, title = "Power: TLI_0C")
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "nfi_rejection_rate_0A", power = .80, title = "Power: NFI_0A")
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "nfi_rejection_rate_0C", power = .80, title = "Power: NFI_0C")
dev.off()

pdf("power_tt_0.3_blav.pdf", height = 2*4, width = 3*44) # create PDF
par(mfrow = c(2, 3))
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "cfi_rejection_rate_0A", power = .80, title = "Power: CFI_0A")
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "cfi_rejection_rate_0C", power = .80, title = "Power: CFI_0C")
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "tli_rejection_rate_0A", power = .80, title = "Power: TLI_0A")
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "tli_rejection_rate_0C", power = .80, title = "Power: TLI_0C")
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "nfi_rejection_rate_0A", power = .80, title = "Power: NFI_0A")
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "nfi_rejection_rate_0C", power = .80, title = "Power: NFI_0C")
dev.off()



pdf("type1error01_blav_2.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 2))
plot_results(data_blav_fit, miss_size = 0, miss_type = "none", fit_index = "brmsea_rejection_rate", power = .05, title = "Type 1 Error: BRMSEA")
plot_results(data_blav_fit, miss_size = 0, miss_type = "none", fit_index = "bgammahat_rejection_rate", power = .05, title = "Type 1 Error: BGammaHat")
plot_results(data_blav_fit, miss_size = 0, miss_type = "none", fit_index = "adjbgammahat_rejection_rate", power = .05, title = "Type 1 Error: adj BGammaHat")
plot_results(data_blav_fit, miss_size = 0, miss_type = "none", fit_index = "bmc_rejection_rate", power = .05, title = "Type 1 Error: BMc")
dev.off()

pdf("power_tt1_0.6_blav_2.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 2))
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "brmsea_rejection_rate", power = .80, title = "Power: BRMSEA")
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "bgammahat_rejection_rate", power = .80, title = "Power: BGammaHat")
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "adjbgammahat_rejection_rate", power = .80, title = "Power: adj BGammaHat")
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index = "bmc_rejection_rate", power = .80, title = "Power: BMc")
dev.off()

pdf("power_tt_0.6_blav_2.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 2))
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "brmsea_rejection_rate", power = .80, title = "Power: BRMSEA")
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "bgammahat_rejection_rate", power = .80, title = "Power: BGammaHat")
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "adjbgammahat_rejection_rate", power = .80, title = "Power: adj BGammaHat")
plot_results(data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "bmc_rejection_rate", power = .80, title = "Power: BMc")
dev.off()

pdf("power_tt1_0.3_blav_2.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 2))
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "brmsea_rejection_rate", power = .80, title = "Power: BRMSEA")
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "bgammahat_rejection_rate", power = .80, title = "Power: BGammaHat")
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "adjbgammahat_rejection_rate", power = .80, title = "Power: adj BGammaHat")
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index = "bmc_rejection_rate", power = .80, title = "Power: BMc")
dev.off()

pdf("power_tt_0.3_blav_2.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 2))
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "brmsea_rejection_rate", power = .80, title = "Power: BRMSEA")
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "bgammahat_rejection_rate", power = .80, title = "Power: BGammaHat")
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "adjbgammahat_rejection_rate", power = .80, title = "Power: adj BGammaHat")
plot_results(data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "bmc_rejection_rate", power = .80, title = "Power: BMc")
dev.off()





#### Plot version 2 ####
# Plots both power and type 1
# X-Axis: N_p (log-scaled)
# Y-Axis: % rejected
# Different lines for different number of time points

plot_results_pow_typ1 = function(df, miss_size, miss_type, fit_index, timepoint, power = 0.8, type1 = 0.05, title){
  
  # Initialize plot
  par(new = FALSE)
  
  condition_none = df$Type_misfit == "none" & df$Size_misfit == 0 & df$n_t_f == timepoint
  condition_missfit = df$Type_misfit == miss_type & df$Size_misfit == miss_size & df$n_t_f == timepoint
  
  
  plot(sort(log(df$n_p[condition_missfit])), df[condition_missfit , fit_index][rank(df$n_p[condition_missfit])], ylim = c(0, 1), 
       xlim = log(c(min(df$n_p), max(df$n_p))), col = "royalblue", pch = 1, type = "b",
       axes = F, xlab = "N", ylab = "% rejected", main = paste(title, "time point:", timepoint, miss_type, miss_size))
  
  par(new = T)
  plot(sort(log(df$n_p[condition_none])), df[condition_none , fit_index][rank(df$n_p[condition_none])], ylim = c(0, 1), 
       xlim = log(c(min(df$n_p), max(df$n_p))), col = "darkorange", pch = 2, type = "b",
       axes = F, xlab = "N", ylab = "% rejected")
  
  axis(1, at = log(N_p), labels = N_p)
  axis(2)
  abline(h = c(type1, power), lty = 3)
  
  legend("right", c("Power", "Type I error"), lty = 1, col = c("royalblue", "darkorange"), pch = 1:2, bty = "n")
  
}


#### Lavaan ####
pdf("power_type1_lav.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2,3)) 
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "pvalue_rejection_rate", timepoint = 2, title = "RMSEA")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "chisq_df_rejection_rate", timepoint = 2, title = "chi-squared")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "cfi_rejection_rate", timepoint = 2, title = "CFI")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate", timepoint = 2, title = "TLI")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "rmsea_rejection_rate", timepoint = 2, title = "RMSEA")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "srmr_rejection_rate", timepoint = 2, title = "SRMR")
dev.off()

pdf("power_type1_lav_pval.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2,3)) 
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "pvalue_rejection_rate", timepoint = 1, title = "p-value")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "pvalue_rejection_rate", timepoint = 2, title = "p-value")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "pvalue_rejection_rate", timepoint = 3, title = "p-value")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "pvalue_rejection_rate", timepoint = 5, title = "p-value")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "pvalue_rejection_rate", timepoint = 10, title = "p-value")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "pvalue_rejection_rate", timepoint = 15, title = "p-value")
dev.off()

pdf("power_type1_lav_rmsea.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2,3)) 
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "rmsea_rejection_rate", timepoint = 1, title = "RMSEA")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "rmsea_rejection_rate", timepoint = 2, title = "RMSEA")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "rmsea_rejection_rate", timepoint = 3, title = "RMSEA")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "rmsea_rejection_rate", timepoint = 5, title = "RMSEA")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "rmsea_rejection_rate", timepoint = 10, title = "RMSEA")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "rmsea_rejection_rate", timepoint = 15, title = "RMSEA")
dev.off()

pdf("power_type1_lav_tli.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2,3)) 
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate", timepoint = 1, title = "TLI")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate", timepoint = 2, title = "TLI")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate", timepoint = 3, title = "TLI")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate", timepoint = 5, title = "TLI")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate", timepoint = 10, title = "TLI")
plot_results_pow_typ1(df = data_lav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate", timepoint = 15, title = "TLI")
dev.off()



#### Blavaan ####
pdf("power_type1_blav.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2,3))
plot_results_pow_typ1(df = data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "brmsea_rejection_rate", timepoint = 2, title = "BRMSEA")
plot_results_pow_typ1(df = data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "bgammahat_rejection_rate", timepoint = 2, title = "BGammaHat")
plot_results_pow_typ1(df = data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "adjbgammahat_rejection_rate", timepoint = 2, title = "adj BGammaHat")
plot_results_pow_typ1(df = data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "bmc_rejection_rate", timepoint = 2, title = "BMc")
dev.off()

pdf("power_type1_blav_tli.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2,3))
plot_results_pow_typ1(df = data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate_0C", timepoint = 1, title = "TLI")
plot_results_pow_typ1(df = data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate_0C", timepoint = 2, title = "TLI")
plot_results_pow_typ1(df = data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate_0C", timepoint = 3, title = "TLI")
plot_results_pow_typ1(df = data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate_0C", timepoint = 5, title = "TLI")
plot_results_pow_typ1(df = data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate_0C", timepoint = 10, title = "TLI")
plot_results_pow_typ1(df = data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index = "tli_rejection_rate_0C", timepoint = 15, title = "TLI")
dev.off()



plot_results_pow_typ1_all = function(df, fit_index, timepoint, power = 0.8, type1 = 0.05, title){
  
  # Initialize plot
  par(new = FALSE)
  
  condition_none = df$Type_misfit == "none" & df$Size_misfit == 0 & df$n_t_f == timepoint
  condition_tt_0.3 = df$Type_misfit == "tt" & df$Size_misfit == 0.3 & df$n_t_f == timepoint
  condition_tt_0.6 = df$Type_misfit == "tt" & df$Size_misfit == 0.6 & df$n_t_f == timepoint
  condition_tt1_0.3 = df$Type_misfit == "tt1" & df$Size_misfit == 0.3 & df$n_t_f == timepoint
  condition_tt1_0.6 = df$Type_misfit == "tt1" & df$Size_misfit == 0.6 & df$n_t_f == timepoint
  
  
  plot(sort(log(df$n_p[condition_tt_0.3])), df[condition_tt_0.3 , fit_index][rank(df$n_p[condition_tt_0.3])], ylim = c(0, 1), 
       xlim = log(c(min(df$n_p), max(df$n_p))), col = "royalblue", pch = 1, type = "b",
       axes = F, xlab = "N", ylab = "% rejected", main = paste(title, "time points:", timepoint))
  
  par(new = T)
  plot(sort(log(df$n_p[condition_tt_0.6])), df[condition_tt_0.6 , fit_index][rank(df$n_p[condition_tt_0.6])], ylim = c(0, 1), 
       xlim = log(c(min(df$n_p), max(df$n_p))), col = "hotpink", pch = 2, type = "b",
       axes = F, xlab = "N", ylab = "% rejected")
  
  if (timepoint != 1){
    par(new = T)
    plot(sort(log(df$n_p[condition_tt1_0.3])), df[condition_tt1_0.3 , fit_index][rank(df$n_p[condition_tt1_0.3])], ylim = c(0, 1), 
         xlim = log(c(min(df$n_p), max(df$n_p))), col = "palegreen4", pch = 3, type = "b",
         axes = F, xlab = "N", ylab = "% rejected")
    
    par(new = T)
    plot(sort(log(df$n_p[condition_tt1_0.6])), df[condition_tt1_0.6 , fit_index][rank(df$n_p[condition_tt1_0.6])], ylim = c(0, 1), 
         xlim = log(c(min(df$n_p), max(df$n_p))), col = "darkorange", pch = 4, type = "b",
         axes = F, xlab = "N", ylab = "% rejected")
    
  }
  
  par(new = T)
  plot(sort(log(df$n_p[condition_none])), df[condition_none , fit_index][rank(df$n_p[condition_none])], ylim = c(0, 1), 
       xlim = log(c(min(df$n_p), max(df$n_p))), col = "black", pch = 5, type = "b",
       axes = F, xlab = "N", ylab = "% rejected")
  
  
  axis(1, at = log(N_p), labels = N_p)
  axis(2)
  abline(h = c(type1, power), lty = 3)
  
  legend("right", c("Power tt 0.3", "Power tt 0.6", "Power tt1 0.3", "Power tt1 0.6", "Type I error"), lty = 1, col = c("royalblue", "hotpink", "palegreen4", "darkorange", "black"), pch = 1:2, bty = "n")
}

pdf("power_type1_lav_pval_all.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2,3)) 
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "pvalue_rejection_rate", timepoint = 1, title = "p-value")
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "pvalue_rejection_rate", timepoint = 2, title = "p-value")
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "pvalue_rejection_rate", timepoint = 3, title = "p-value")
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "pvalue_rejection_rate", timepoint = 5, title = "p-value")
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "pvalue_rejection_rate", timepoint = 10, title = "p-value")
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "pvalue_rejection_rate", timepoint = 15, title = "p-value")
dev.off()

pdf("power_type1_lav_rmsea_all.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2,3)) 
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "rmsea_rejection_rate", timepoint = 1, title = "RMSEA")
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "rmsea_rejection_rate", timepoint = 2, title = "RMSEA")
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "rmsea_rejection_rate", timepoint = 3, title = "RMSEA")
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "rmsea_rejection_rate", timepoint = 5, title = "RMSEA")
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "rmsea_rejection_rate", timepoint = 10, title = "RMSEA")
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "rmsea_rejection_rate", timepoint = 15, title = "RMSEA")
dev.off()

pdf("power_type1_lav_TLI_all.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2,3)) 
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "tli_rejection_rate", timepoint = 1, title = "TLI")
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "tli_rejection_rate", timepoint = 2, title = "TLI")
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "tli_rejection_rate", timepoint = 3, title = "TLI")
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "tli_rejection_rate", timepoint = 5, title = "TLI")
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "tli_rejection_rate", timepoint = 10, title = "TLI")
plot_results_pow_typ1_all(df = data_lav_fit, fit_index = "tli_rejection_rate", timepoint = 15, title = "TLI")
dev.off()


pdf("power_type1_blav_TLI_all.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2,3)) 
plot_results_pow_typ1_all(df = data_blav_fit, fit_index = "tli_rejection_rate_0C", timepoint = 1, title = "TLI")
plot_results_pow_typ1_all(df = data_blav_fit, fit_index = "tli_rejection_rate_0C", timepoint = 2, title = "TLI")
plot_results_pow_typ1_all(df = data_blav_fit, fit_index = "tli_rejection_rate_0C", timepoint = 3, title = "TLI")
plot_results_pow_typ1_all(df = data_blav_fit, fit_index = "tli_rejection_rate_0C", timepoint = 5, title = "TLI")
plot_results_pow_typ1_all(df = data_blav_fit, fit_index = "tli_rejection_rate_0C", timepoint = 10, title = "TLI")
plot_results_pow_typ1_all(df = data_blav_fit, fit_index = "tli_rejection_rate_0C", timepoint = 15, title = "TLI")
dev.off()

pdf("power_type1_blav_BMc_all.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2,3)) 
plot_results_pow_typ1_all(df = data_blav_fit, fit_index = "bmc_rejection_rate", timepoint = 1, title = "BMc")
plot_results_pow_typ1_all(df = data_blav_fit, fit_index = "bmc_rejection_rate", timepoint = 2, title = "BMc")
plot_results_pow_typ1_all(df = data_blav_fit, fit_index = "bmc_rejection_rate", timepoint = 3, title = "BMc")
plot_results_pow_typ1_all(df = data_blav_fit, fit_index = "bmc_rejection_rate", timepoint = 5, title = "BMc")
plot_results_pow_typ1_all(df = data_blav_fit, fit_index = "bmc_rejection_rate", timepoint = 10, title = "BMc")
plot_results_pow_typ1_all(df = data_blav_fit, fit_index = "bmc_rejection_rate", timepoint = 15, title = "BMc")
dev.off()




#### Plot version 3 ####
#### Plotting lavaan vs. blavaan ####

plot_results_lav_blav = function(df_lav, df_blav, miss_size, miss_type, fit_index_lav, fit_index_blav, timepoint, power = 0.8, type1 = 0.05, title){
  
  # Initialize plot
  par(new = FALSE)
  
  # lavaan
  condition_none = df_lav$Type_misfit == "none" & df_lav$Size_misfit == 0 & df_lav$n_t_f == timepoint
  condition_missfit = df_lav$Type_misfit == miss_type & df_lav$Size_misfit == miss_size & df_lav$n_t_f == timepoint
  
  plot(sort(log(df_lav$n_p[condition_missfit])), df_lav[condition_missfit , fit_index_lav][rank(df_lav$n_p[condition_missfit])], ylim = c(0, 1), 
       xlim = log(c(min(df_lav$n_p), max(df_lav$n_p))), col = "lightskyblue", pch = 1, type = "b",
       axes = F, xlab = "N", ylab = "% rejected", main = paste(title, miss_type, miss_size, "time points:", timepoint, sep = " "))
  
  par(new = T)
  plot(sort(log(df_lav$n_p[condition_none])), df_lav[condition_none , fit_index_lav][rank(df_lav$n_p[condition_none])], ylim = c(0, 1), 
       xlim = log(c(min(df_lav$n_p), max(df_lav$n_p))), col = "darkorange2", pch = 2, type = "b",
       axes = F, xlab = "", ylab = "", main = paste(title, miss_type, miss_size, "time points:", timepoint, sep = " "))
  
  
  # blavaan
  condition_none = df_blav$Type_misfit == "none" & df_blav$Size_misfit == 0 & df_blav$n_t_f == timepoint
  condition_missfit = df_blav$Type_misfit == miss_type & df_blav$Size_misfit == miss_size & df_blav$n_t_f == timepoint
  
  par(new = T)
  plot(sort(log(df_blav$n_p[condition_missfit])), df_blav[condition_missfit , fit_index_blav][rank(df_blav$n_p[condition_missfit])], ylim = c(0, 1), 
       xlim = log(c(min(df_blav$n_p), max(df_blav$n_p))), col = "royalblue", pch = 3, type = "b",
       axes = F, xlab = "", ylab = "", main = paste(title, miss_type, miss_size, "time points:", timepoint, sep = " "))
  
  par(new = T)
  plot(sort(log(df_blav$n_p[condition_none])), df_blav[condition_none , fit_index_blav][rank(df_blav$n_p[condition_none])], ylim = c(0, 1), 
       xlim = log(c(min(df_blav$n_p), max(df_blav$n_p))), col = "goldenrod1", pch = 4, type = "b",
       axes = F, xlab = "", ylab = "", main = paste(title, miss_type, miss_size, "time points:", timepoint, sep = " "))
  

  axis(1, log(N_p), N_p)
  axis(2)
  abline(h = c(type1, power), lty = 3)
  
  legend("right", 
         c("Blavaan (Power)", "Lavaan (Power)", "Blavaan (Type I)", "Lavaan (Type I)"), 
         lty = 1, col = c("royalblue", "lightskyblue", "goldenrod1", "darkorange2"), pch = 1:4, bty = "n")  
}



pdf("plot_results_lav_blav_rmsea.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2, 2)) 
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index_lav = "rmsea_rejection_rate", 
                      fit_index_blav = "brmsea_rejection_rate", timepoint = 2, power = 0.8, type1 = 0.05, title = "RMSEA")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index_lav = "rmsea_rejection_rate", 
                      fit_index_blav = "brmsea_rejection_rate", timepoint = 2, power = 0.8, type1 = 0.05, title = "RMSEA")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index_lav = "rmsea_rejection_rate", 
                      fit_index_blav = "brmsea_rejection_rate", timepoint = 2, power = 0.8, type1 = 0.05, title = "RMSEA")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index_lav = "rmsea_rejection_rate", 
                      fit_index_blav = "brmsea_rejection_rate", timepoint = 2, power = 0.8, type1 = 0.05, title = "RMSEA")
dev.off()


pdf("plot_results_lav_blav_CFI.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2, 2)) 
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index_lav = "cfi_rejection_rate", 
                      fit_index_blav = "cfi_rejection_rate_0C", timepoint = 2, power = 0.8, type1 = 0.05, title = "CFI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index_lav = "cfi_rejection_rate", 
                      fit_index_blav = "cfi_rejection_rate_0C", timepoint = 2, power = 0.8, type1 = 0.05, title = "CFI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index_lav = "cfi_rejection_rate", 
                      fit_index_blav = "cfi_rejection_rate_0C", timepoint = 2, power = 0.8, type1 = 0.05, title = "CFI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index_lav = "cfi_rejection_rate", 
                      fit_index_blav = "cfi_rejection_rate_0C", timepoint = 2, power = 0.8, type1 = 0.05, title = "CFI")
dev.off()
  
pdf("plot_results_lav_blav_TLI_2.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2, 2)) 
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 2, power = 0.8, type1 = 0.05, title = "TLI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 2, power = 0.8, type1 = 0.05, title = "TLI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 2, power = 0.8, type1 = 0.05, title = "TLI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 2, power = 0.8, type1 = 0.05, title = "TLI")
dev.off()

pdf("plot_results_lav_blav_TLI_1.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2, 2)) 
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 1, power = 0.8, type1 = 0.05, title = "TLI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 1, power = 0.8, type1 = 0.05, title = "TLI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 1, power = 0.8, type1 = 0.05, title = "TLI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 1, power = 0.8, type1 = 0.05, title = "TLI")
dev.off()

pdf("plot_results_lav_blav_TLI_5.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2, 2)) 
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 5, power = 0.8, type1 = 0.05, title = "TLI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 5, power = 0.8, type1 = 0.05, title = "TLI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 5, power = 0.8, type1 = 0.05, title = "TLI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 5, power = 0.8, type1 = 0.05, title = "TLI")
dev.off()

pdf("plot_results_lav_blav_TLI_10.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2, 2)) 
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 10, power = 0.8, type1 = 0.05, title = "TLI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 10, power = 0.8, type1 = 0.05, title = "TLI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 10, power = 0.8, type1 = 0.05, title = "TLI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 10, power = 0.8, type1 = 0.05, title = "TLI")
dev.off()

pdf("plot_results_lav_blav_TLI_15.pdf", height = 2*4, width = 3*4)
par(mfrow = c(2, 2)) 
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 15, power = 0.8, type1 = 0.05, title = "TLI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 15, power = 0.8, type1 = 0.05, title = "TLI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 15, power = 0.8, type1 = 0.05, title = "TLI")
plot_results_lav_blav(data_lav_fit, data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index_lav = "tli_rejection_rate", 
                      fit_index_blav = "tli_rejection_rate_0C", timepoint = 15, power = 0.8, type1 = 0.05, title = "TLI")
dev.off()





#### Plot version 4 ####
#### Plotting results for different null models ####

# Plots null model 0A vs null model 0C
# X-Axis: N_p (log-scaled)
# Y-Axis: % rejected
# Different lines for different number of time points

plot_results_null_models = function(df, miss_size, miss_type, fit_index_0A, fit_index_0C, timepoint, power = 0.8, type1 = 0.05, title){
  
  # Initialize plot
  par(new = F)
  
  condition_none = df$Type_misfit == "none" & df$Size_misfit == 0 & df$n_t_f == timepoint
  condition_missfit = df$Type_misfit == miss_type & df$Size_misfit == miss_size & df$n_t_f == timepoint
  
  
  plot(sort(log(df$n_p[condition_missfit])), df[condition_missfit , fit_index_0A][rank(df$n_p[condition_missfit])], ylim = c(0, 1), 
       xlim = log(c(min(df$n_p), max(df$n_p))), col = "lightskyblue", pch = 1, type = "b",
       axes = F, xlab = "N", ylab = "", main = paste(title, miss_type, miss_size, "time points:", timepoint))
  
  par(new = T)
  plot(sort(log(df$n_p[condition_none])), df[condition_none , fit_index_0A][rank(df$n_p[condition_none])], ylim = c(0, 1), 
       xlim = log(c(min(df$n_p), max(df$n_p))), col = "darkorange2", pch = 2, type = "b",
       axes = F, xlab = "", ylab = "")
  
  
  par(new = T)
  plot(sort(log(df$n_p[condition_missfit])), df[condition_missfit , fit_index_0C][rank(df$n_p[condition_missfit])], ylim = c(0, 1), 
       xlim = log(c(min(df$n_p), max(df$n_p))), col = "royalblue", pch = 3, type = "b",
       axes = F, xlab = "", ylab = "")
  
  par(new = T)
  plot(sort(log(df$n_p[condition_none])), df[condition_none , fit_index_0C][rank(df$n_p[condition_none])], ylim = c(0, 1), 
       xlim = log(c(min(df$n_p), max(df$n_p))), col = "goldenrod1", pch = 4, type = "b",
       axes = F, xlab = "", ylab = "")
  
  
  axis(1, at = log(N_p), labels = N_p)
  axis(2)
  abline(h = c(power, type1), lty = 3)
  
  legend("right", c("Power Null Model 0C", "Power Null Model 0A", "Type 1 Error Model 0C", "Type 1 Error Model 0A"), lty = 1, col = c("royalblue", "lightskyblue", "goldenrod1", "darkorange2"), pch = 1:4, bty = "n")
  
}


pdf("power_type1_0.6_null_models_blav.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 3))
plot_results_null_models(data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index_0A = "cfi_rejection_rate_0A", fit_index_0C = "cfi_rejection_rate_0C",  timepoint = 2, power = .80, title = "CFI")
plot_results_null_models(data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index_0A = "tli_rejection_rate_0A", fit_index_0C = "tli_rejection_rate_0C", timepoint = 2, power = .80, title = "TLI")
plot_results_null_models(data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index_0A = "nfi_rejection_rate_0A", fit_index_0C = "nfi_rejection_rate_0C", timepoint = 2, power = .80, title = "NFI")
plot_results_null_models(data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index_0A = "cfi_rejection_rate_0A", fit_index_0C = "cfi_rejection_rate_0C", timepoint = 2, power = .80, title = "CFI")
plot_results_null_models(data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index_0A = "tli_rejection_rate_0A", fit_index_0C = "tli_rejection_rate_0C", timepoint = 2, power = .80, title = "TLI")
plot_results_null_models(data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index_0A = "nfi_rejection_rate_0A", fit_index_0C = "nfi_rejection_rate_0C", timepoint = 2, power = .80, title = "NFI")
dev.off()

pdf("power_type1_0.3_null_models_blav.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 3))
plot_results_null_models(data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index_0A = "cfi_rejection_rate_0A", fit_index_0C = "cfi_rejection_rate_0C",  timepoint = 2, power = .80, title = "CFI")
plot_results_null_models(data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index_0A = "tli_rejection_rate_0A", fit_index_0C = "tli_rejection_rate_0C", timepoint = 2, power = .80, title = "TLI")
plot_results_null_models(data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index_0A = "nfi_rejection_rate_0A", fit_index_0C = "nfi_rejection_rate_0C", timepoint = 2, power = .80, title = "NFI")
plot_results_null_models(data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index_0A = "cfi_rejection_rate_0A", fit_index_0C = "cfi_rejection_rate_0C", timepoint = 2, power = .80, title = "CFI")
plot_results_null_models(data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index_0A = "tli_rejection_rate_0A", fit_index_0C = "tli_rejection_rate_0C", timepoint = 2, power = .80, title = "TLI")
plot_results_null_models(data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index_0A = "nfi_rejection_rate_0A", fit_index_0C = "nfi_rejection_rate_0C", timepoint = 2, power = .80, title = "NFI")
dev.off()

pdf("power_type1_TLI_null_models_blav.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 3))
plot_results_null_models(data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index_0A = "tli_rejection_rate_0A", fit_index_0C = "tli_rejection_rate_0C",  timepoint = 1, power = .80, title = "TLI")
plot_results_null_models(data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index_0A = "tli_rejection_rate_0A", fit_index_0C = "tli_rejection_rate_0C", timepoint = 2, power = .80, title = "TLI")
plot_results_null_models(data_blav_fit, miss_size = 0.3, miss_type = "tt1", fit_index_0A = "tli_rejection_rate_0A", fit_index_0C = "tli_rejection_rate_0C", timepoint = 3, power = .80, title = "TLI")
plot_results_null_models(data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index_0A = "tli_rejection_rate_0A", fit_index_0C = "tli_rejection_rate_0C", timepoint = 5, power = .80, title = "TLI")
plot_results_null_models(data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index_0A = "tli_rejection_rate_0A", fit_index_0C = "tli_rejection_rate_0C", timepoint = 10, power = .80, title = "TLI")
plot_results_null_models(data_blav_fit, miss_size = 0.3, miss_type = "tt", fit_index_0A = "tli_rejection_rate_0A", fit_index_0C = "tli_rejection_rate_0C", timepoint = 15, power = .80, title = "TLI")
dev.off()



pdf("power_type1_0.6_null_models_blav_3.pdf", height = 2*4, width = 3*4) # create PDF
par(mfrow = c(2, 3))
plot_results_null_models(data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index_0A = "cfi_rejection_rate_0A", fit_index_0C = "cfi_rejection_rate_0C",  timepoint = 3, power = .80, title = "CFI")
plot_results_null_models(data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index_0A = "tli_rejection_rate_0A", fit_index_0C = "tli_rejection_rate_0C", timepoint = 3, power = .80, title = "TLI")
plot_results_null_models(data_blav_fit, miss_size = 0.6, miss_type = "tt1", fit_index_0A = "nfi_rejection_rate_0A", fit_index_0C = "nfi_rejection_rate_0C", timepoint = 3, power = .80, title = "NFI")
plot_results_null_models(data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index_0A = "cfi_rejection_rate_0A", fit_index_0C = "cfi_rejection_rate_0C", timepoint = 3, power = .80, title = "CFI")
plot_results_null_models(data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index_0A = "tli_rejection_rate_0A", fit_index_0C = "tli_rejection_rate_0C", timepoint = 3, power = .80, title = "TLI")
plot_results_null_models(data_blav_fit, miss_size = 0.6, miss_type = "tt", fit_index_0A = "nfi_rejection_rate_0A", fit_index_0C = "nfi_rejection_rate_0C", timepoint = 3, power = .80, title = "NFI")
dev.off()




################################################################################
# Tables 
################################################################################

#### Table 1 ####
# Table with minimum number of sample size to have minimum power of 0.8 and max type 1 error of 0.05
# for each number of timepoints and type of misfit


# function only considering power
table_min_sample_power = function(df, fit_index, power = 0.8) {
  # Create an empty matrix to store the results
  results <- matrix(NA, nrow = length(N_t), ncol = 4)
  colnames(results) <- c("tt_0.3", "tt_0.6", "tt1_0.3", "tt1_0.6")
  rownames(results) <- N_t
  
  # Loop through each time point
  for (timepoint in N_t) {
    # Filter data for the current time point
    df_timepoint <- df %>% filter(n_t == timepoint)
    
    # Define conditions for each misfit type and size
    condition_tt_0.3 <- df_timepoint$Type_misfit == "tt" & df_timepoint$Size_misfit == 0.3
    condition_tt_0.6 <- df_timepoint$Type_misfit == "tt" & df_timepoint$Size_misfit == 0.6
    condition_tt1_0.3 <- df_timepoint$Type_misfit == "tt1" & df_timepoint$Size_misfit == 0.3
    condition_tt1_0.6 <- df_timepoint$Type_misfit == "tt1" & df_timepoint$Size_misfit == 0.6
    
    # Filter data for each condition
    data_tt_0.3 <- df_timepoint %>% filter(condition_tt_0.3)
    data_tt_0.6 <- df_timepoint %>% filter(condition_tt_0.6)
    data_tt1_0.3 <- df_timepoint %>% filter(condition_tt1_0.3)
    data_tt1_0.6 <- df_timepoint %>% filter(condition_tt1_0.6)
    
    # Function to calculate minimum sample size or return NA if no valid values
    calculate_min_sample <- function(data, fit_index, power) {
      valid_samples <- data$n_p[data[[fit_index]] >= power]
      if (length(valid_samples) > 0) {
        return(min(valid_samples))
      } else {
        return(NA)  # Return NA if no valid sample sizes exist
      }
    }
    
    # Calculate minimum sample size for each condition
    min_n_tt_0.3 <- calculate_min_sample(data_tt_0.3, fit_index, power)
    min_n_tt_0.6 <- calculate_min_sample(data_tt_0.6, fit_index, power)
    min_n_tt1_0.3 <- calculate_min_sample(data_tt1_0.3, fit_index, power)
    min_n_tt1_0.6 <- calculate_min_sample(data_tt1_0.6, fit_index, power)
    
    # Store the results in the matrix
    results[timepoint, ] <- c(min_n_tt_0.3, min_n_tt_0.6, min_n_tt1_0.3, min_n_tt1_0.6)
  }
  
  return(results)
}

result_table_power_pval = table_min_sample_power(data_lav_fit, fit_index = "pvalue_rejection_rate")
result_table_power_rmsea = table_min_sample_power(data_lav_fit, fit_index = "rmsea_rejection_rate")
result_table_power_tli = table_min_sample_power(data_lav_fit, fit_index = "tli_rejection_rate")
result_table_power_chisq = table_min_sample_power(data_lav_fit, fit_index = "chisq_df_rejection_rate")

print(result_table_power_pval)
print(result_table_power_rmsea)
print(result_table_power_tli)
print(result_table_power_chisq)


# function considering both power and type 1 error
table_min_sample_pow_typ1 = function(df, fit_index, power = 0.8, type1 = 0.05) {
  # Create an empty matrix to store the results
  results <- matrix(NA, nrow = length(unique(df$n_t)), ncol = 4)
  colnames(results) <- c("tt_0.3", "tt_0.6", "tt1_0.3", "tt1_0.6")
  rownames(results) <- unique(df$n_t)
  
  # Loop through each time point
  for (timepoint in unique(df$n_t)) {
    # Filter data for the current time point
    df_timepoint <- df %>% filter(n_t == timepoint)
    
    # Define conditions for each misfit type and size
    condition_none <- df_timepoint$Type_misfit == "none"
    condition_tt_0.3 <- df_timepoint$Type_misfit == "tt" & df_timepoint$Size_misfit == 0.3
    condition_tt_0.6 <- df_timepoint$Type_misfit == "tt" & df_timepoint$Size_misfit == 0.6
    condition_tt1_0.3 <- df_timepoint$Type_misfit == "tt1" & df_timepoint$Size_misfit == 0.3
    condition_tt1_0.6 <- df_timepoint$Type_misfit == "tt1" & df_timepoint$Size_misfit == 0.6
    
    # Filter data for each condition
    data_none <- df_timepoint %>% filter(condition_none)
    data_tt_0.3 <- df_timepoint %>% filter(condition_tt_0.3)
    data_tt_0.6 <- df_timepoint %>% filter(condition_tt_0.6)
    data_tt1_0.3 <- df_timepoint %>% filter(condition_tt1_0.3)
    data_tt1_0.6 <- df_timepoint %>% filter(condition_tt1_0.6)
    
    # Function to calculate minimum sample size or return NA if no valid values
    calculate_min_sample <- function(data, data_typ1, fit_index, power, type1) {
      valid_samples <- data$n_p[data[[fit_index]] >= power & data_typ1[[fit_index]] <= type1]
      if (length(valid_samples) > 0) {
        return(min(valid_samples))
      } else {
        return(NA)  # Return NA if no valid sample sizes exist
      }
    }
    
    # Calculate minimum sample size for each condition
    min_n_tt_0.3 <- calculate_min_sample(data_tt_0.3, data_none, fit_index, power, type1)
    min_n_tt_0.6 <- calculate_min_sample(data_tt_0.6, data_none, fit_index, power, type1)
    min_n_tt1_0.3 <- calculate_min_sample(data_tt1_0.3, data_none, fit_index, power, type1)
    min_n_tt1_0.6 <- calculate_min_sample(data_tt1_0.6, data_none, fit_index, power, type1)
    
    # Store the results in the matrix
    # results[timepoint, ] <- c(min_n_tt_0.3, min_n_tt_0.6, min_n_tt1_0.3, min_n_tt1_0.6)
    
    timepoint_index <- which(rownames(results) == timepoint)
    
    # Store the results in the matrix using timepoint_index
    if (length(timepoint_index) > 0) {
      results[timepoint_index, ] <- c(min_n_tt_0.3, min_n_tt_0.6, min_n_tt1_0.3, min_n_tt1_0.6)
    } else {
      warning(paste("Timepoint", timepoint, "not found in the results matrix"))
    }
  }
  
  return(results)
}

  

result_table_power_type1_pval_lav = table_min_sample_pow_typ1(data_lav_fit, fit_index = "pvalue_rejection_rate")
result_table_power_type1_rmsea_lav = table_min_sample_pow_typ1(data_lav_fit, fit_index = "rmsea_rejection_rate")
result_table_power_type1_tli_lav = table_min_sample_pow_typ1(data_lav_fit, fit_index = "tli_rejection_rate")
result_table_power_type1_chisq_lav = table_min_sample_pow_typ1(data_lav_fit, fit_index = "chisq_df_rejection_rate")

# print(result_table_power_type1_pval_lav)
# print(result_table_power_type1_rmsea_lav)
write.csv(result_table_power_type1_tli_lav, "result_table_power_type1_tli_lav.txt")
# print(result_table_power_type1_chisq_lav)


result_table_power_type1_brmsea_blav = table_min_sample_pow_typ1(data_blav_fit, fit_index = "brmsea_rejection_rate")
result_table_power_type1_tli_blav = table_min_sample_pow_typ1(data_blav_fit, fit_index = "tli_rejection_rate_0C")
result_table_power_type1_cfi_blav = table_min_sample_pow_typ1(data_blav_fit, fit_index = "cfi_rejection_rate_0C")

# print(result_table_power_type1_brmsea_blav)
write.csv(result_table_power_type1_tli_blav, "result_table_power_type1_tli_blav.txt")
# print(result_table_power_type1_cfi_blav)




################################################################################
# Heatmaps 
################################################################################


#### Heatmap ####
heatmap_power_type1 = function(df, miss_size, miss_type, fit_index) {
  df_filtered = df %>% filter(Type_misfit == miss_type & Size_misfit == miss_size)
  ggplot(df_filtered, aes(x = n_p, y = n_t_f, fill = .data[[fit_index]])) +
    geom_tile() +
    scale_x_log10() +
    labs(x = "N_p (log scale)", y = "N_t", fill = fit_index, title = paste("Heatmap of", fit_index)) +
    theme_minimal()
}

heatmap = heatmap_power_type1(data_lav_fit, 0, "none", "pvalue_rejection_rate")



# Function to create a heatmap of power or type 1 error
heatmap_power_type1 = function(df, miss_size, miss_type, fit_index, type = "power") {
  # Filter data for the specific misfit type and size
  df_filtered = df %>%
    filter(Type_misfit == miss_type & Size_misfit == miss_size)
  
  # Aggregate data to calculate power or type 1 error for each N_p and N_t
  df_aggregated = df_filtered %>%
    group_by(n_p, n_t_f) %>%
    summarise(
      power = mean(.data[[fit_index]] >= 0.8, na.rm = TRUE),  # Power: % of rejections >= 0.8
      type1_error = mean(.data[[fit_index]], na.rm = TRUE)    # Type 1 error: mean rejection rate
    )
  
  # Choose the metric to plot (power or type 1 error)
  if (type == "power") {
    fill_var = df_aggregated$power
    fill_label = "Power"
  } else if (type == "type1") {
    fill_var = df_aggregated$type1_error
    fill_label = "Type 1 Error"
  } else {
    stop("Invalid type. Choose 'power' or 'type1'.")
  }
  
  # Create the heatmap
  ggplot(df_aggregated, aes(x = n_p, y = n_t_f, fill = fill_var)) +
    geom_tile() +  # Create tiles for the heatmap
    scale_x_log10() +  # Log-scale for sample size (N_p)
    scale_fill_gradient(low = "white", high = "steelblue") +  # Color gradient
    labs(
      x = "Sample Size (N_p, log scale)",
      y = "Number of Time Points (N_t)",
      fill = fill_label,
      title = paste("Heatmap of", fit_index, "for", miss_type, "Misfit (Size =", miss_size, ")")
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
}

# Example usage
heatmap_power_type1(data_lav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "pvalue_rejection_rate", type = "power")
heatmap_power_type1(data_lav_fit, miss_size = 0.3, miss_type = "tt", fit_index = "pvalue_rejection_rate", type = "type1")






#### HEATMAP FOR CORRELATION BETWEEN FIT INDICES ####





