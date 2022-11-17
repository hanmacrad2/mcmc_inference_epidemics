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

#ETA CREDIBLE INTERVALS
plot.new(); par(mfrow=c(1,1))
ETA_CREDIBLE_INTERVALS(mcmcI$eta_matrix, simX1$eta_vec, lwdX = 2)

#ETA plots
eta_start = 1; eta_step = 17
PLOT_ETA(dataI, mcmcI$eta_matrix, simX1$eta_vec, seedX)

#Plot eta II
eta_start = eta_start + eta_step; print(eta_start)
eta_start = 94
PLOT_ETA(dataI, mcmcI$eta_matrix, simX1$eta_vec, seedX, eta_start = eta_start)

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

#ETA CREDIBLE INTERVALS
plot.new(); par(mfrow=c(1,1))
ETA_CREDIBLE_INTERVALS(mcmcII$eta_matrix, simX2$eta_vec, lwdX = 2)

#ETA plots
eta_start = 1; eta_step = 17
PLOT_ETA(dataII, mcmcII$eta_matrix, simX2$eta_vec, seedX)

#Plot eta II
eta_start = eta_start + eta_step
print(eta_start)
eta_start = 94
PLOT_ETA(dataII, mcmcII$eta_matrix, simX2$eta_vec, seedX, eta_start = eta_start)

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

#ETA CREDIBLE INTERVALS
plot.new(); par(mfrow=c(1,1))
ETA_CREDIBLE_INTERVALS(mcmcIII$eta_matrix, simX3$eta_vec, lwdX = 2)

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

#ETA CREDIBLE INTERVALS
plot.new(); par(mfrow=c(1,1))
ETA_CREDIBLE_INTERVALS(mcmcIV$eta_matrix, simX4$eta_vec, lwdX = 2)


#****************
#*DATA V (alpha = 0.9, k = 0.16)
#****************
seedX = 5
set.seed(seedX)
simX5 = SIMULATE_ETA(alphaX = 1.0)
dataV = simX5$epidemic_data
plot.ts(dataV, ylab = 'Daily Infection count',
        main = 'Individual reproduction number model') #DATA II LOOKS GOOD; SEED = 7. seed 4 (data I)
saveRDS(simX5, file = paste0(OUTER_FOLDER, '/simX5_', seedX, '.rds' ))

#LIKELIHOOD
loglike5 = LOG_LIKELIHOOD_ETA(dataV, c(1.0,0.16), simX5$eta_vec)
loglike5

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcV = MCMC_ADAPTIVE_ETA(dataV, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcV$time_elap = time_elap
saveRDS(mcmcV, file = paste0(OUTER_FOLDER, '/mcmcV', seedX, '.rds' ))

#PLOT 
dfV = PLOT_MCMC_ETA_GRID(dataV, mcmcV, seedX, simX5$eta_vec, loglike5)

#ETA CREDIBLE INTERVALS
plot.new(); par(mfrow=c(1,1))
ETA_CREDIBLE_INTERVALS(mcmcV$eta_matrix, simX5$eta_vec, lwdX = 2)

