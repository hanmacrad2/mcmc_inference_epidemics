#RUN MODEL EVIDENCE

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

log_model_evidence10 = LOG_MODEL_EVIDENCE(mcmc10$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence10))
vec_lme[i] = log_model_evidence10

#Plot log model evidence
plot.new()
par(mfrow = c(1,1))
plot(vec_lme,  pch = 19,
     ylab = 'Log of Model Evidence',
     main = 'Log of Model Evidence. MCMC 10 runs (100k) - Dataset I')


#**********************
#* DATASET II
#* 

#STORAGE
num_iters = 3
vec_lme2 = vector("numeric", length = num_iters)

#MCMC 1 
i = 1
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc11 = MCMC_ADAPTIVE_ETA(dataII, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc11$time_elap = time_elap

log_model_evidence11 = LOG_MODEL_EVIDENCE(mcmc11$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence11))
vec_lme2[i] = log_model_evidence11
vec_lme2

#MCMC 2
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc22 = MCMC_ADAPTIVE_ETA(dataII, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc22$time_elap = time_elap

log_model_evidence22 = LOG_MODEL_EVIDENCE(mcmc22$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence22))
vec_lme2[i] = log_model_evidence22
vec_lme2

#MCMC 3
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc33 = MCMC_ADAPTIVE_ETA(dataII, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc33$time_elap = time_elap

log_model_evidence33 = LOG_MODEL_EVIDENCE(mcmc33$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence33))
vec_lme2[i] = log_model_evidence33
vec_lme2

#***********
#* TAKE 2

#STORAGE
num_iters = 5
vec_lme2 = vector("numeric", length = num_iters)
#vec_lme2[i] = log_model_evidence33
#vec_lme2

#MCMC 1 
i = 1
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc11 = MCMC_ADAPTIVE_ETA(dataII, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc11$time_elap = time_elap

log_model_evidence11 = LOG_MODEL_EVIDENCE(mcmc11$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence11))
vec_lme2[i] = log_model_evidence11
vec_lme2

#MCMC 2
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc22 = MCMC_ADAPTIVE_ETA(dataII, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc22$time_elap = time_elap

log_model_evidence22 = LOG_MODEL_EVIDENCE(mcmc22$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence22))
vec_lme2[i] = log_model_evidence22
vec_lme2

#MCMC 3
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc33 = MCMC_ADAPTIVE_ETA(dataII, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc33$time_elap = time_elap

log_model_evidence33 = LOG_MODEL_EVIDENCE(mcmc33$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence33))
vec_lme2[i] = log_model_evidence33
vec_lme2

#MCMC 4
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc44 = MCMC_ADAPTIVE_ETA(dataII, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc44$time_elap = time_elap

log_model_evidence44 = LOG_MODEL_EVIDENCE(mcmc44$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence44))
vec_lme2[i] = log_model_evidence44
vec_lme2

#MCMC 5
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc55 = MCMC_ADAPTIVE_ETA(dataII, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc55$time_elap = time_elap

log_model_evidence55 = LOG_MODEL_EVIDENCE(mcmc55$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence55))
vec_lme2[i] = log_model_evidence55
vec_lme2

#MCMC 6
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc66 = MCMC_ADAPTIVE_ETA(dataII, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc66$time_elap = time_elap

log_model_evidence66 = LOG_MODEL_EVIDENCE(mcmc66$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence66))
vec_lme2[i] = log_model_evidence66
vec_lme2

#MCMC 7
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc77 = MCMC_ADAPTIVE_ETA(dataII, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc77$time_elap = time_elap

log_model_evidence77 = LOG_MODEL_EVIDENCE(mcmc77$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence77))
vec_lme2[i] = log_model_evidence77
vec_lme2

#MCMC 8
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc88 = MCMC_ADAPTIVE_ETA(dataII, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc88$time_elap = time_elap

log_model_evidence88 = LOG_MODEL_EVIDENCE(mcmc88$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence88))
vec_lme2[i] = log_model_evidence88
vec_lme2