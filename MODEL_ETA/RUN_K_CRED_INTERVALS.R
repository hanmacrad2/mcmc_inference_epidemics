#Librarys
library(RChronoModel)
library(plotrix)

#CREDIBLE INTERVAL
get_lower_ci <- function(mcmc_output){
  
  lower_interval = CredibleInterval(mcmc_output$nu_params_matrix[, 1], level = 0.95)[[2]]
  
  return(lower_interval)
}

get_upper_ci <- function(mcmc_output){
  
  upper_interval = CredibleInterval(mcmc_output$nu_params_matrix[, 1], level = 0.95)[[3]]
  
  return(upper_interval)
}

#PLOT CREDIBLE INTERVAL 

#VARIABLES
vec_k = c(0.1, 0.16, 0.25, 0.5, 1, 1.5, 2, 4, 8, 10, 100)
num_iters = length(vec_k)
vec_k_means = vector("numeric", length = num_iters)
vec_k_lower = vector("numeric", length = num_iters); vec_k_upper = vector("numeric", length = num_iters)

#PLOTS
plot.new()
par(mfrow=c(1,1))

#DATA 1
seedX = 1
set.seed(seedX)
#simX112 = SIMULATE_NU(k = vec_k[seedX])
#dataX112 = simX112$epidemic_data
#plot.ts(dataX112)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcX111 = MCMC_ADAPTIVE_ETA(dataX112, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcX111$time_elap = time_elap

print(paste0('mean 1', mean(mcmcX111$nu_params_matrix[, 1])))
vec_k_means[seedX] = mean(mcmcX111$nu_params_matrix[, 1])
vec_k_lower[seedX] = get_lower_ci(mcmcX111) 
vec_k_upper[seedX] = get_upper_ci(mcmcX111) 

#*********************
#DATA 2
seedX = seedX + 1
set.seed(seedX)
# simX22 = SIMULATE_NU(k = vec_k[seedX])
# dataX22 = simX22$epidemic_data
# plot.ts(dataX22)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcX222 = MCMC_ADAPTIVE_ETA(dataX22, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcX222$time_elap = time_elap
vec_k_means[seedX] = mean(mcmcX222$nu_params_matrix[, 1])
vec_k_lower[seedX] = get_lower_ci(mcmcX222) 
vec_k_upper[seedX] = get_upper_ci(mcmcX222) 

#****************
#DATA 3
seedX = seedX + 1
set.seed(seedX)
# simX33 = SIMULATE_NU(k = vec_k[seedX])
# dataX33 = simX33$epidemic_data
# plot.ts(dataX33)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcX333 = MCMC_ADAPTIVE_ETA(dataX33, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcX333$time_elap = time_elap

#Results
vec_k_means[seedX] = mean(mcmcX333$nu_params_matrix[, 1])
vec_k_lower[seedX] = get_lower_ci(mcmcX333) 
vec_k_upper[seedX] = get_upper_ci(mcmcX333) 

#*****************************
#DATA: 
seedX = seedX + 1
set.seed(seedX)
# simX44 = SIMULATE_NU(k = vec_k[seedX])
# dataX44 = simX44$epidemic_data
# plot.ts(dataX44)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcX444 = MCMC_ADAPTIVE_ETA(dataX44, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcX444$time_elap = time_elap

#Results
vec_k_means[seedX] = mean(mcmcX444$nu_params_matrix[, 1])
vec_k_lower[seedX] = get_lower_ci(mcmcX444) 
vec_k_upper[seedX] = get_upper_ci(mcmcX444)

#****************
#DATA 5
seedX = seedX + 1
set.seed(seedX)
# simX55 = SIMULATE_NU(k = vec_k[seedX])
# dataX55 = simX55$epidemic_data
# plot.ts(dataX55)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcX555 = MCMC_ADAPTIVE_ETA(dataX55, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcX555$time_elap = time_elap

#Results
vec_k_means[seedX] = mean(mcmcX555$nu_params_matrix[, 1])
vec_k_lower[seedX] = get_lower_ci(mcmcX555) 
vec_k_upper[seedX] = get_upper_ci(mcmcX555)

#*****************************
#* DATA 6
seedX = seedX + 1
set.seed(seedX)
# simX66 = SIMULATE_NU(k = vec_k[seedX])
# dataX66 = simX66$epidemic_data
# plot.ts(dataX66)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcX666 = MCMC_ADAPTIVE_ETA(dataX66, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcX666$time_elap = time_elap

#Results
vec_k_means[seedX] = mean(mcmcX666$nu_params_matrix[, 1]) 
vec_k_lower[seedX] = get_lower_ci(mcmcX666) 
vec_k_upper[seedX] = get_upper_ci(mcmcX666)

#*****************************
#* DATA 7
seedX = seedX + 1
set.seed(seedX)
# simX77 = SIMULATE_NU(k = vec_k[seedX])
# dataX77 = simX77$epidemic_data
# plot.ts(dataX77)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcX777 = MCMC_ADAPTIVE_ETA(dataX77, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcX777$time_elap = time_elap

#Results
vec_k_means[seedX] = mean(mcmcX777$nu_params_matrix[, 1]) 
vec_k_lower[seedX] = get_lower_ci(mcmcX777) 
vec_k_upper[seedX] = get_upper_ci(mcmcX777)

#*****************************
#* DATA 8
seedX = seedX + 1
set.seed(seedX)
# simX88 = SIMULATE_NU(k = vec_k[seedX])
# dataX88 = simX88$epidemic_data
# plot.ts(dataX88)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcX888 = MCMC_ADAPTIVE_ETA(dataX88, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcX888$time_elap = time_elap

#Results
vec_k_means[seedX] = mean(mcmcX888$nu_params_matrix[, 1]) 
vec_k_lower[seedX] = get_lower_ci(mcmcX888) 
vec_k_upper[seedX] = get_upper_ci(mcmcX888)

#*****************************
#* DATA 9
seedX = seedX + 1
set.seed(seedX)
# simX99 = SIMULATE_NU(k = vec_k[seedX])
# dataX99 = simX99$epidemic_data
# plot.ts(dataX99)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcX999 = MCMC_ADAPTIVE_ETA(dataX99, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcX999$time_elap = time_elap

#Results
vec_k_means[seedX] = mean(mcmcX999$nu_params_matrix[, 1]) 
vec_k_lower[seedX] = get_lower_ci(mcmcX999) 
vec_k_upper[seedX] = get_upper_ci(mcmcX999)

#*****************************
#* DATA 10
seedX = seedX + 1
set.seed(seedX)
# simX10 = SIMULATE_NU(k = vec_k[seedX])
# dataX10 = simX10$epidemic_data
# plot.ts(dataX10)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcX10 = MCMC_ADAPTIVE_ETA(dataX10, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcX10$time_elap = time_elap

#Results
vec_k_means[seedX] = mean(mcmcX10$nu_params_matrix[, 1]) 
vec_k_lower[seedX] = get_lower_ci(mcmcX10) 
vec_k_upper[seedX] = get_upper_ci(mcmcX10)

#*****************************
#* DATA 11
seedX = seedX + 1
set.seed(seedX)
# simX11 = SIMULATE_NU(k = vec_k[seedX])
# dataX11 = simX11$epidemic_data
# plot.ts(dataX11)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcX11 = MCMC_ADAPTIVE_ETA(dataX11, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcX11$time_elap = time_elap

#Results
vec_k_means[seedX] = mean(mcmcX11$nu_params_matrix[, 1]) 
vec_k_lower[seedX] = get_lower_ci(mcmcX11) 
vec_k_upper[seedX] = get_upper_ci(mcmcX11)


#********************************
#PLOT CREDIBLE INTERVALS
plotCI(vec_k, vec_k_means, ui= vec_k_upper, li= vec_k_lower,
       xlab = 'k', ylab = 'k mean', main = 'k MCMC Posterior Mean &
       95 % Credible intervals. N MCMC = 100k',
       xlim = c(min(vec_alpha), max(vec_alpha)), lwd = 2)
lines(vec_k, vec_k, col = 'red', lwd = 2)