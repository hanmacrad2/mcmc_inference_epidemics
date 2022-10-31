#RUN MODEL EVIDENCE

#STORAGE
num_iters = 10
vec_lme = vector(numeric, length = num_iters)

#MCMC 1
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc1 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc1$time_elap = time_elap

log_model_evidence = log_model_evidence(mcmc1$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence))
vec_lme[i] = log_model_evidence

#MCMC 2
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc2 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc2$time_elap = time_elap

log_model_evidence = log_model_evidence(mcmc2$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence))
vec_lme[i] = log_model_evidence

#MCMC 3
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc3 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc3$time_elap = time_elap

log_model_evidence = log_model_evidence(mcmc3$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence))
vec_lme[i] = log_model_evidence

#MCMC 4
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc3 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc4$time_elap = time_elap

log_model_evidence = log_model_evidence(mcmc4$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence))
vec_lme[i] = log_model_evidence

#MCMC 5
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc5 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc5$time_elap = time_elap

log_model_evidence = log_model_evidence(mcmc5$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence))
vec_lme[i] = log_model_evidence

#MCMC 6
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc6 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc6$time_elap = time_elap

log_model_evidence = log_model_evidence(mcmc6$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence))
vec_lme[i] = log_model_evidence

#MCMC 7
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc7 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc7$time_elap = time_elap

log_model_evidence = log_model_evidence(mcmc7$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence))
vec_lme[i] = log_model_evidence

#MCMC 8
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc8 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc8$time_elap = time_elap

log_model_evidence = log_model_evidence(mcmc8$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence))
vec_lme[i] = log_model_evidence

#MCMC 9
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc9 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc9$time_elap = time_elap

log_model_evidence = log_model_evidence(mcmc9$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence))
vec_lme[i] = log_model_evidence

#MCMC 10
i = i + 1; print(paste0('i:', i))
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc10 = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc10$time_elap = time_elap

log_model_evidence = log_model_evidence(mcmc10$log_like_vec)
print(paste0('log_model_evidence =  ', log_model_evidence))
vec_lme[i] = log_model_evidence