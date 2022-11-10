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
vec_alpha = c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7)
num_iters = 10
vec_means = vector("numeric", length = num_iters)
vec_lower = vector("numeric", length = num_iters); vec_upper = vector("numeric", length = num_iters)

#PLOTS
plot.new()
par(mfrow=c(1,1))

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
vec_lower[seedX] = get_lower_ci(mcmc111) 
vec_upper[seedX] = get_upper_ci(mcmc111) 

#*********************
#DATA 2
seedX = seedX + 1
set.seed(seedX)
sim22 = SIMULATE_NU(alphaX = vec_alpha[seedX])
data22 = sim22$epidemic_data
plot.ts(data22)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc222 = MCMC_ADAPTIVE_ETA(data22, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc222$time_elap = time_elap
vec_means[seedX] = mean(mcmc222$nu_params_matrix[, 1])
vec_lower[seedX] = get_lower_ci(mcmc222) 
vec_upper[seedX] = get_upper_ci(mcmc222) 

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

#Results
vec_means[seedX] = mean(mcmc333$nu_params_matrix[, 1])
vec_lower[seedX] = get_lower_ci(mcmc333) 
vec_upper[seedX] = get_upper_ci(mcmc333) 

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

#Results
vec_means[seedX] = mean(mcmc44$nu_params_matrix[, 1])
vec_lower[seedX] = get_lower_ci(mcmc44) 
vec_upper[seedX] = get_upper_ci(mcmc44)

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

#Results
vec_means[seedX] = mean(mcmc555$nu_params_matrix[, 1])
vec_lower[seedX] = get_lower_ci(mcmc555) 
vec_upper[seedX] = get_upper_ci(mcmc555)

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

#Results
vec_means[seedX] = mean(mcmc666$nu_params_matrix[, 1]) 
vec_lower[seedX] = get_lower_ci(mcmc666) 
vec_upper[seedX] = get_upper_ci(mcmc666)

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

#Results
vec_means[seedX] = mean(mcmc777$nu_params_matrix[, 1]) 
vec_lower[seedX] = get_lower_ci(mcmc777) 
vec_upper[seedX] = get_upper_ci(mcmc777)

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

#Results
vec_means[seedX] = mean(mcmc888$nu_params_matrix[, 1]) 
vec_lower[seedX] = get_lower_ci(mcmc888) 
vec_upper[seedX] = get_upper_ci(mcmc888)

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

#Results
vec_means[seedX] = mean(mcmc999$nu_params_matrix[, 1]) 
vec_lower[seedX] = get_lower_ci(mcmc999) 
vec_upper[seedX] = get_upper_ci(mcmc999)

#*****************************
#* DATA 10
seedX = seedX + 1
seedX = 10
set.seed(seedX)
sim10 = SIMULATE_NU(alphaX = vec_alpha[seedX])
data10 = sim10$epidemic_data
plot.ts(data10)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc10 = MCMC_ADAPTIVE_ETA(data10, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc10$time_elap = time_elap

#Results
vec_means[seedX] = mean(mcmc10$nu_params_matrix[, 1]) 
vec_lower[seedX] = get_lower_ci(mcmc10) 
vec_upper[seedX] = get_upper_ci(mcmc10)


#********************************
#PLOT CREDIBLE INTERVALS
require(plotrix)
plotCI(vec_alpha, vec_means, ui= vec_upper, li= vec_lower,
       xlab = 'alpha', ylab = 'alpha mean', main = 'Alpha MCMC Posterior Mean &
       95 % Credible intervals. N MCMC = 100k',
       xlim = c(min(vec_alpha), max(vec_alpha)), lwd = 2)
lines(vec_alpha, vec_alpha, col = 'red', lwd = 2)