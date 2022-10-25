#********************************************************
#1. INDIVIDUAL R0 - SIMULATE DATA & SAMPLE USING MCMC                         
#********************************************************
source("~/Github/epidemic_modelling/helper_functions.R")

#SIMULATION FUNCTION
SIMULATE_NU = function(num_days = 110, alphaX = 1.2, k = 0.16,
                       shape_gamma = 6, scale_gamma = 1) {
  
  'Simulate from the Negative Binomial model'
  
  #INTIALISE VECTORS
  x = vector('numeric', num_days); x[1] = 2
  eta_vec = vector('numeric', num_days); 
  
  #INFECTIOUSNESS (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #DAYS OF THE EPIDEMIC
  for (t in 2:num_days) {
    
    #ETA (t-1)
    eta_vec[t-1] <- rgamma(1, shape = x[t-1]*k, scale = alphaX/k) #Draw eta from previous time step
    #INFECTIVITY
    infectivity = rev(prob_infect[1:t]) 
    #POISSON; OFFSPRINT DISTRIBUTION
    total_rate = sum(eta_vec*infectivity) #DOT PRODUCT
    x[t] = rpois(1, total_rate)
    
  }
  return(list(epidemic_data = x, eta_vec = eta_vec))
}


#**************************************
#* APPLY + SIMULATIONS 
#***************************************

#DATA I
seedX = 4
set.seed(seedX)
simX1 = SIMULATE_NU()
dataI = simX1$epidemic_data
plot.ts(dataI) #DATA II LOOKS GOOD; SEED = 7. seed 4 (data I)

#LIKELIHOOD
loglike1 = LOG_LIKELIHOOD_NU(dataI, c(1.2,0.16), simX1$eta_vec)
loglike1

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcIB = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcIB$time_elap = time_elap

#PLOT
dfIB = PLOT_MCMC_ETA_GRID(dataI, mcmcIB, seedX, simX1$eta_vec, loglike1)

#****************
#*DATA II
#****************
seedX = 7
set.seed(seedX)
simX2 = SIMULATE_NU()
dataII = simX2$epidemic_data
plot.ts(dataII) #DATA II LOOKS GOOD; SEED = 7. seed 4 (data I)

#LIKELIHOOD
loglike2 = LOG_LIKELIHOOD_NU(dataII, c(1.2,0.16), simX2$eta_vec)
loglike2

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcII = MCMC_ADAPTIVE_ETA(dataII, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcII$time_elap = time_elap

#PLOT 
dfII = PLOT_MCMC_ETA_GRID(dataII, mcmcII, seedX, simX2$eta_vec, loglike2)


