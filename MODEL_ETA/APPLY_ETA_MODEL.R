#********************************************************
#1. INDIVIDUAL R0 - SIMULATE DATA & SAMPLE USING MCMC                         
#********************************************************
library(MASS)
source("MCMC_ADAPTIVE_ETA.R")
source("PLOT_MCMC_ETA_GRID.R")
source("../helper_functions.R")
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_individual_nu/v1"
ifelse(!dir.exists(file.path(OUTER_FOLDER)), dir.create(file.path(OUTER_FOLDER), recursive = TRUE), FALSE)

#LOAD DATA

#**************************************
#* APPLY + SIMULATIONS 
#***************************************

#DATA I
seedX = 4
set.seed(seedX)
# simX1 = SIMULATE_ETA()
# dataI = simX1$epidemic_data
# plot.ts(dataI) #DATA II LOOKS GOOD; SEED = 7. seed 4 (data I)
# saveRDS(simX1, file = paste0(OUTER_FOLDER, '/simX1_', seedX, '.rds' ))
# 
# #LIKELIHOOD
# loglike1 = LOG_LIKELIHOOD_ETA(dataI, c(1.2,0.16), simX1$eta_vec)
# loglike1

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcI = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcI$time_elap = time_elap
saveRDS(mcmcI, file = paste0(OUTER_FOLDER, '/mcmcI', seedX, '.rds' ))

#PLOT
dfI = PLOT_MCMC_ETA_GRID(dataI, mcmcI, seedX, simX1$eta_vec, loglike1)

# #ETA plots
# eta_start = 1; eta_step = 18
# eta1 = mcmcI$eta_matrix
# PLOT_ETA(dataI, mcmcI$eta_matrix, simX1$eta_vec, seedX) 
# 
# #Plot eta II
# eta_start = eta_start + eta_step
# eta_start = eta_start + 1
# print(eta_start)
# PLOT_ETA(dataI, mcmcI$eta_matrix, simX1$eta_vec, seedX, eta_start = eta_start)

#BAYES FACTOR
#bf1 = get_bayes_factor(mcmcI$log_like_vec, mcmcIB$log_like_vec)
#bf1

#****************
#*DATA II
#****************
seedX = 7
set.seed(seedX)
# simX2 = SIMULATE_ETA()
# dataII = simX2$epidemic_data
# plot.ts(dataII, ylab = 'Daily Infection count',
#         main = 'Individual reproduction number model') #DATA II LOOKS GOOD; SEED = 7. seed 4 (data I)
# saveRDS(simX2, file = paste0(OUTER_FOLDER, '/simX2_', seedX, '.rds' ))
# 
# #LIKELIHOOD
# loglike2 = LOG_LIKELIHOOD_ETA(dataII, c(1.2,0.16), simX2$eta_vec)
# loglike2

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcII = MCMC_ADAPTIVE_ETA(dataII, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcII$time_elap = time_elap
saveRDS(mcmcII, file = paste0(OUTER_FOLDER, '/mcmcII', seedX, '.rds' ))

#PLOT
dfII = PLOT_MCMC_ETA_GRID(dataII, mcmcII, seedX, simX2$eta_vec, loglike2)
# 
# #ETA plots
# eta_start = 1; eta_step = 18
# eta2 = mcmcII$eta_matrix
# PLOT_ETA(dataII, mcmcII$eta_matrix, simX2$eta_vec, seedX) 
# 
# #Plot eta II
# eta_start = eta_start + eta_step
# print(eta_start)
# PLOT_ETA(dataII, mcmcII$eta_matrix, simX2$eta_vec, seedX, eta_start = eta_start)

#****************
#*DATA III (alpha = 0.9, k = 0.16)
#****************
seedX = 3
set.seed(seedX)
# simX3 = SIMULATE_ETA(alphaX = 0.9)
# dataIII = simX3$epidemic_data
# plot.ts(dataIII, ylab = 'Daily Infection count',
#         main = 'Individual reproduction number model') #DATA II LOOKS GOOD; SEED = 7. seed 4 (data I)
# saveRDS(simX3, file = paste0(OUTER_FOLDER, '/simX3_', seedX, '.rds' ))
# 
# #LIKELIHOOD
# loglike3 = LOG_LIKELIHOOD_ETA(dataIII, c(0.9,0.16), simX3$eta_vec)
# loglike3

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcIII = MCMC_ADAPTIVE_ETA(dataIII, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcIII$time_elap = time_elap
saveRDS(mcmcIII, file = paste0(OUTER_FOLDER, '/mcmcIII', seedX, '.rds' ))

#PLOT 
dfIII = PLOT_MCMC_ETA_GRID(dataIII, mcmcIII, seedX, simX3$eta_vec, loglike3)

#****************
#*DATA IV
#****************
seedX = 4
set.seed(seedX)
# simX4 = SIMULATE_ETA(k = 1.0)
# dataIV = simX4$epidemic_data
# plot.ts(dataIV, ylab = 'Daily Infection count',
#         main = 'Individual reproduction number model') #DATA II LOOKS GOOD; SEED = 7. seed 4 (data I)
# saveRDS(simX4, file = paste0(OUTER_FOLDER, '/simX4_', seedX, '.rds' ))
# 
# #LIKELIHOOD
# loglike4 = LOG_LIKELIHOOD_ETA(dataIV, c(1.2,1.0), simX4$eta_vec)
# loglike4

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcIV = MCMC_ADAPTIVE_ETA(dataIV, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcIV$time_elap = time_elap
saveRDS(mcmcIV, file = paste0(OUTER_FOLDER, '/mcmcIV', seedX, '.rds' ))

#PLOT 
dfIV = PLOT_MCMC_ETA_GRID(dataIV, mcmcIV, seedX, simX4$eta_vec, loglike4)

#**************
#* ETA CREDIBLE INTERVALS PRACTICE
library(RChronoModel)
library(plotrix)

ETA_CREDIBLE_INTERVALS <- function(eta_matrix, eta_true){
  
  #Create a vector of means across columns
  eta_means = colMeans(eta_matrix)
  #Upper & lower limits
  ci = get_ci_matrix(eta_matrix) 
  print(ci$vec_lower[10])
  print(ci$vec_upper[10])
  print(mean(ci$vec_upper - ci$vec_lower))
  #Plot
  plotCI(eta_true, eta_means, ui = ci$vec_upper, li = ci$vec_lower,
         xlab = 'True eta', ylab = 'Eta', main = 'Eta MCMC Posterior Mean &
       95 % Credible intervals. N MCMC = 100k', lwd = 2) #xlim = c(min(vec_alpha), max(vec_alpha)))
  lines(eta_true, eta_true, col = 'red', lwd = 2)
  
  
}

#PLOT ETA CREDIBLE INTERVALS
ETA_CREDIBLE_INTERVALS <- function(eta_matrix, eta_true, lwdX = 2){
  
  #Create a vector of means across columns
  eta_means = colMeans(eta_matrix)
  #Upper & lower limits
  ci = get_ci_matrix(eta_matrix) 
  
  #Plot
  plotCI(seq_along(eta_true), eta_means, ui = ci$vec_upper, li = ci$vec_lower,
         xlab = 'Day of Epidemic', ylab = 'Eta', main = 'Eta MCMC Posterior Mean &
       95 % Credible intervals. Red (True) ', lwd = lwdX, pch = 16) #xlim = c(min(vec_alpha), max(vec_alpha)))
  lines(eta_true, col = 'red', lwd = lwdX, pch = 16)
  points(eta_true, col = 'red', lwd = lwdX, pch = 16)
  
  
}

#Apply
ETA_CREDIBLE_INTERVALS(mcmcI$eta_matrix, simX1$eta_vec)

plot.new(); par(mfrow = c(1,1))
eta_matrix_2 = mcmcII$eta_matrix
ETA_CREDIBLE_INTERVALS(mcmcII$eta_matrix, simX2$eta_vec)

plotCI(c)