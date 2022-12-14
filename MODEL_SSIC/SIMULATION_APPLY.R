#********************************************************
#1. INDIVIDUAL R0 - SIMULATE DATA & SAMPLE USING MCMC                         
#********************************************************
library(MASS)
source("~/Github/epidemic_modelling/helper_functions.R")
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_individual_nu/v1"
ifelse(!dir.exists(file.path(OUTER_FOLDER)), dir.create(file.path(OUTER_FOLDER), recursive = TRUE), FALSE)

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

#RUN MCMC MULTIPLE TIMES
run_mcmc(num_iters = 10){
  
  vec_lme = vector(numeric, length = num_iters)
  for (i in 1:num){
    
    print(paste0('ITER ', i))
    #Run mcmc
    mcmcX = mcmc() #STORE MCMC?? 
    log_model_evidence = log_model_evidence(mcmcX$log_like_vec)
    print(paste0('log_model_evidence =  ', log_model_evidence))
    vec_lme[i] = log_model_evidence
    
  }
  
  return(vec_lme)
}

#**************************************
#* APPLY + SIMULATIONS 
#***************************************

#DATA I
seedX = 4
set.seed(seedX)
#simX1 = SIMULATE_NU()
#dataI = simX1$epidemic_data
plot.ts(dataI) #DATA II LOOKS GOOD; SEED = 7. seed 4 (data I)
plot.ts(simX1$eta_vec)

#LIKELIHOOD
loglike1 = LOG_LIKELIHOOD_ETA(dataI, c(1.2,0.16), simX1$eta_vec)
loglike1

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcI = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcI$time_elap = time_elap

#PLOT
dfI = PLOT_MCMC_ETA_GRID(dataI, mcmcI, seedX, simX1$eta_vec, loglike1)

#ETA plots
eta_start = 1; eta_step = 18
eta1 = mcmcI$eta_matrix
PLOT_ETA(dataI, mcmcI$eta_matrix, simX1$eta_vec, seedX) 

#Plot eta II
eta_start = eta_start + eta_step
eta_start = eta_start + 1
print(eta_start)
PLOT_ETA(dataI, mcmcI$eta_matrix, simX1$eta_vec, seedX, eta_start = eta_start)

#BAYES FACTOR
#bf1 = get_bayes_factor(mcmcI$log_like_vec, mcmcIB$log_like_vec)
#bf1

#****************
#*DATA II
#****************
seedX = 7
set.seed(seedX)
#simX2 = SIMULATE_NU()
#dataII = simX2$epidemic_data
#plot.ts(dataII, ylab = 'Daily Infection count',
#        main = 'Individual reproduction number model') #DATA II LOOKS GOOD; SEED = 7. seed 4 (data I)

#LIKELIHOOD
#loglike2 = LOG_LIKELIHOOD_ETA(dataII, c(1.2,0.16), simX2$eta_vec)
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

#ETA plots
eta_start = 1; eta_step = 18
eta2 = mcmcII$eta_matrix
PLOT_ETA(dataII, mcmcII$eta_matrix, simX2$eta_vec, seedX) 

#Plot eta II
eta_start = eta_start + eta_step
print(eta_start)
PLOT_ETA(dataII, mcmcII$eta_matrix, simX2$eta_vec, seedX, eta_start = eta_start)

#****************
#*DATA III (alpha = 0.9, k = 0.16)
#****************
seedX = 3
set.seed(seedX)
#simX3 = SIMULATE_NU(alphaX = 0.9)
#dataIII = simX3$epidemic_data
#plot.ts(dataIII, ylab = 'Daily Infection count',
#        main = 'Individual reproduction number model') #DATA II LOOKS GOOD; SEED = 7. seed 4 (data I)

#LIKELIHOOD
loglike3 = LOG_LIKELIHOOD_ETA(dataIII, c(0.9,0.16), simX3$eta_vec)
loglike3

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcIII = MCMC_ADAPTIVE_ETA(dataIII, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcIII$time_elap = time_elap

#PLOT 
dfIII = PLOT_MCMC_ETA_GRID(dataIII, mcmcIII, seedX, simX3$eta_vec, loglike3)


#****************
#*DATA IV
#****************
seedX = 4
set.seed(seedX)
#simX4 = SIMULATE_NU(k = 1.0)
#dataIV = simX4$epidemic_data
plot.ts(dataIV, ylab = 'Daily Infection count',
        main = 'Individual reproduction number model') #DATA II LOOKS GOOD; SEED = 7. seed 4 (data I)

#LIKELIHOOD
loglike4 = LOG_LIKELIHOOD_ETA(dataIV, c(1.2,1.0), simX4$eta_vec)
loglike4

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcIV = MCMC_ADAPTIVE_ETA(dataIV, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcIV$time_elap = time_elap

#PLOT 
dfIV = PLOT_MCMC_ETA_GRID(dataIV, mcmcIV, seedX, simX4$eta_vec, loglike4)


#****************
#*DATA I
#****************

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcIB = MCMC_ADAPTIVE_ETA(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcIB$time_elap = time_elap

#PLOT
dfIB = PLOT_MCMC_ETA_GRID(dataI, mcmcIB, seedX, simX1$eta_vec, loglike1)

#BAYES FACTOR (Different start points)
bf2 = get_bayes_factor(mcmcII$log_like_vec, mcmcIIB$log_like_vec)
bf2

#**********************************
#* DATA SIMULATIONS
#***********************************

seedX = 8; #seed 18, 23, 32 interesting
seedX = seedX + 1
set.seed(seedX)
simXXX = SIMULATE_NU()
dataXX = simXXX$epidemic_data
plot.ts(dataXX)

#****************
#*DATA III
#****************
seedX = 32
set.seed(seedX)
simX3 = SIMULATE_NU()
dataIII = simX3$epidemic_data
plot.ts(dataIII) #DATA II LOOKS GOOD; SEED = 7. seed 4 (data I)

#LIKELIHOOD
loglike3 = LOG_LIKELIHOOD_ETA(dataIII, c(1.2,0.16), simX3$eta_vec)
loglike3

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmcIII = MCMC_ADAPTIVE_ETA(dataIII, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmcIII$time_elap = time_elap

#PLOT 
dfIII = PLOT_MCMC_ETA_GRID(dataIII, mcmcIII, seedX, simX3$eta_vec, loglike3)

#BAYES FACTOR (Different start points)
#bf2 = get_bayes_factor(mcmcII$log_like_vec, mcmcIIB$log_like_vec)
#bf2

#Test
for(i in 1:3) { 
  nam <- paste("A", i, sep = "")
  assign(nam, rnorm(3)+5)
}

#BQUOTE
bquote(alpha ~ "; Mean of " ~ eta ~ "~ Ga(). Simulated " .(mcmc_specs$simulated$m1))

plot.ts(c(1,2,4,5,3), main = 
          bquote(bold(alpha ~ "; Mean of" ~ eta ~ "~ Ga(). Simulated: " ~ .(mcmc_specs$simulated$m1))))

plot.ts(c(1,2,4,5,3), main = 
bquote(bold(alpha ~ "Prior: exp("~.(priors_list$alpha_prior[1])~")")))

for(day in seq(eta_start, eta_start + 17, by = 1)){
  print(day)
}