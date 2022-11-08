#PLOTS
plot.new()
par(mfrow=c(1,1))

vec_alpha = c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7)
num_iters = 10
vec_means = vector("numeric", length = num_iters)

#DATA 1
seedX = 1
set.seed(seedX)
sim11 = SIMULATE_NU(alphaX = vec_alpha[seedX])
data11 = sim11$epidemic_data
plot.ts(data11)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc111 = MCMC_ADAPTIVE_ETA(data11, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc111$time_elap = time_elap

print(paste0('mean 1', mean(mcmc111$nu_params_matrix[, 1])))
vec_means[seedX] = mean(mcmc111$nu_params_matrix[, 1])

#*********************
#DATA 2
seedX = seedX + 1
#set.seed(seedX)
#sim22 = SIMULATE_NU(alphaX = vec_alpha[seedX])
#data22 = sim22$epidemic_data
plot.ts(data22)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc222 = MCMC_ADAPTIVE_ETA(data22, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc222$time_elap = time_elap
vec_means[seedX] = mean(mcmc222$nu_params_matrix[, 1])

#****************
#DATA 3
#seedX = seedX + 1
seedX = 3
set.seed(seedX)
#sim33 = SIMULATE_NU(alphaX = vec_alpha[seedX])
#data33 = sim33$epidemic_data
#plot.ts(data33)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc333 = MCMC_ADAPTIVE_ETA(data33, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc333$time_elap = time_elap
vec_means[seedX] = mean(mcmc333$nu_params_matrix[, 1])

#*****************************
#DATA: ALPHA = 1.2
seedX = seedX + 1
set.seed(seedX)
#simX1 = SIMULATE_NU()
#dataj = simX1$epidemic_data
#plot.ts(dataj) #DATA II LOOKS GOOD; SEED = 7. seed 4 (data I)

#LIKELIHOOD
#loglike1 = LOG_LIKELIHOOD_ETA(dataI, c(1.2,0.16), simX1$eta_vec)
#loglike1

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc44 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc44$time_elap = time_elap
vec_means[seedX] = mean(mcmc44$nu_params_matrix[, 1])

#****************
#DATA 5
seedX = seedX + 1
set.seed(seedX)
#sim55 = SIMULATE_NU(alphaX = vec_alpha[seedX])
#data55 = sim55$epidemic_data
#plot.ts(data55)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc555 = MCMC_ADAPTIVE_ETA(data55, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc555$time_elap = time_elap
vec_means[seedX] = mean(mcmc555$nu_params_matrix[, 1]) 

#*****************************
#* DATA 6
seedX = seedX + 1
set.seed(seedX)
#sim66 = SIMULATE_NU(alphaX = vec_alpha[seedX])
#data66 = sim66$epidemic_data
#plot.ts(data66)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc666 = MCMC_ADAPTIVE_ETA(data66, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc666$time_elap = time_elap
vec_means[seedX] = mean(mcmc666$nu_params_matrix[, 1]) 

#*****************************
#* DATA 7
seedX = seedX + 1
set.seed(seedX)
#sim77 = SIMULATE_NU(alphaX = vec_alpha[seedX])
#data77 = sim77$epidemic_data
#plot.ts(data77)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc777 = MCMC_ADAPTIVE_ETA(data77, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc777$time_elap = time_elap
vec_means[seedX] = mean(mcmc777$nu_params_matrix[, 1]) 

#*****************************
#* DATA 8
seedX = seedX + 1
set.seed(seedX)
#sim88 = SIMULATE_NU(alphaX = vec_alpha[seedX])
#data88 = sim88$epidemic_data
#plot.ts(data88)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc888 = MCMC_ADAPTIVE_ETA(data88, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc888$time_elap = time_elap
vec_means[seedX] = mean(mcmc888$nu_params_matrix[, 1]) 

#*****************************
#* DATA 9
seedX = seedX + 1
set.seed(seedX)
#sim99 = SIMULATE_NU(alphaX = vec_alpha[seedX])
#data99 = sim99$epidemic_data
#plot.ts(data99)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc999 = MCMC_ADAPTIVE_ETA(data99, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc999$time_elap = time_elap
vec_means[seedX] = mean(mcmc999$nu_params_matrix[, 1]) 

#*****************************
#* DATA 10
seedX = seedX + 1
set.seed(seedX)
#sim10 = SIMULATE_NU(alphaX = vec_alpha[seedX])
#data10 = sim10$epidemic_data
#plot.ts(data10)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc10 = MCMC_ADAPTIVE_ETA(data10, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc10$time_elap = time_elap
vec_means[seedX] = mean(mcmc10$nu_params_matrix[, 1]) 

#************************************************************************
#************************************************************************
#************************************************************************

#PLOT
#dfI = PLOT_MCMC_ETA_GRID(dataI, mcmcI, seedX, simX1$eta_vec, loglike1)

#BAYES FACTOR
bf1 = get_bayes_factor(mcmcI$log_like_vec, mcmcIB$log_like_vec)
bf1

#******************
#******************
#******************

#STORAGE
num_iters = 10
vec_lme = vector("numeric", length = num_iters)

#MCMC 1 
i = 1
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc1 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc1$time_elap = time_elap

log_model_evidence1 = LOG_MODEL_EVIDENCE(mcmc1$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence1))
vec_lme[i] = log_model_evidence1

df11 = PLOT_MCMC_ETA_GRID(dataI, mcmc1, seedX, simX1$eta_vec, loglike1)

#MCMC 2
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc2 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc2$time_elap = time_elap

log_model_evidence2 = LOG_MODEL_EVIDENCE(mcmc2$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence2))
vec_lme[i] = log_model_evidence2

#MCMC 3
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc3 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc3$time_elap = time_elap

log_model_evidence3 = LOG_MODEL_EVIDENCE(mcmc3$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence3))
vec_lme[i] = log_model_evidence3
vec_lme

#MCMC 4
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc4 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc4$time_elap = time_elap

log_model_evidence4 = LOG_MODEL_EVIDENCE(mcmc4$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence4))
vec_lme[i] = log_model_evidence4
vec_lme

#MCMC 5
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc5 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc5$time_elap = time_elap

log_model_evidence5 = LOG_MODEL_EVIDENCE(mcmc5$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence5))
vec_lme[i] = log_model_evidence5

#MCMC 6
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc6 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc6$time_elap = time_elap

log_model_evidence6 = LOG_MODEL_EVIDENCE(mcmc6$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence6))
vec_lme[i] = log_model_evidence6

#MCMC 7
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc7 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc7$time_elap = time_elap

log_model_evidence7 = LOG_MODEL_EVIDENCE(mcmc7$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence7))
vec_lme[i] = log_model_evidence7

#MCMC 8
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc8 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc8$time_elap = time_elap

log_model_evidence8 = LOG_MODEL_EVIDENCE(mcmc8$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence8))
vec_lme[i] = log_model_evidence8

#MCMC 9
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc9 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc9$time_elap = time_elap

log_model_evidence9 = LOG_MODEL_EVIDENCE(mcmc9$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence9))
vec_lme[i] = log_model_evidence9

df9 = PLOT_MCMC_ETA_GRID(dataI, mcmc9, seedX, simX1$eta_vec, loglike1)

#MCMC 10
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc10 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc10$time_elap = time_elap
