#MODEL INDIVIDUAL NU
library(MASS)
source("~/Github/epidemic_modelling/helper_functions.R")
#****************************
#SIMULATION
#*****************************
SIMULATE_NU = function(num_days = 110, alphaX = 1.2, k = 0.16,
                       shape_gamma = 6, scale_gamma = 1) {
  
  'Simulate from the Negative Binomial model'
  
  #INTIALISE VECTORS
  x = vector('numeric', num_days); x[1] = 2
  
  #INFECTIOUSNESS (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  vec_sum_nu_t = vector('numeric', num_days); 
  
  #DAYS OF THE EPIDEMIC
  for (t in 2:num_days) {
    
    #NU: INDIVIDUAL R0; GAMMA(ALPHA, K)
    vec_sum_nu_t[t-1] <- rgamma(1, shape = x[t-1]*k, scale = alphaX/k)
    
    #OFFSPRING; POISSON(). Why t-1?? 
    infectivity = rev(prob_infect[1:t]) #x[1:(t-1)]*rev(prob_infect[1:(t-1)])
    total_rate = sum(vec_sum_nu_t*infectivity)
    #print(total_rate)
    x[t] = rpois(1, total_rate)
    
  }
  return(x)
}


#**********************************************
#LOG LIKELIHOOD
#**********************************************
LOG_LIKELIHOOD_NU <- function(x, nu_params, eta){ #eta - a vector of length x. eta[1] = infectivity of xt[1]
  
  #Params
  num_days = length(x)
  shape_gamma = 6; scale_gamma = 1
  alpha = nu_params[1]; k = nu_params[2]
  count_inf = 0
    
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - 
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  loglike = 0
  
  for (t in 2:num_days) {
    
    #OFFSPRING; POISSON()
    infectivity = rev(prob_infect[1:t-1]) #Current infectivity dependent on people already infected #rev(prob_infect[1:(t-1)]) 
    total_rate = sum(eta[1:(t-1)]*infectivity) #SHOULD ETA BE LENGTH NUM_DAYS --> WORK OUT
    #print(paste0('total rate = ', total_rate))
    
    #dgamma #EXPLAIN DGAMMA??
    eta_prob = dgamma(eta[t-1], shape = x[t-1]*k, scale = alpha/k, log = TRUE)
    if (!is.infinite(eta_prob)) loglike = loglike + eta_prob
    loglike = loglike + x[t]*log(total_rate) - total_rate - lfactorial(x[t]) #Need to include normalizing constant 
    #print(paste0('t:', t, 'loglike 2', loglike))
    #if (is.infinite(loglike))  count_inf = count_inf + 1 #print(paste0('t:', t, 'loglike poisson'))
    
  }
  
  #print(paste0('count_inf: ', count_inf ))
  
  return(loglike)
  
}

#********************************************************
#1. INDIVIDUAL R0 MCMC ADAPTIVE SHAPING                           
#********************************************************
#NOTE NO REFLECTION, NO TRANSFORMS, MORE INTELLIGENT ADAPTATION
MCMC_ADAPTIVE_MODEL_NU <- function(dataX,
                          mcmc_inputs = list(n_mcmc = 100000,
                                             mod_start_points = c(1.2, 0.16),
                                             dim = 2, alpha_star = 0.4, v0 = 100,  #priors_list = list(alpha_prior = c(1, 0), k_prior = c()),
                                             thinning_factor = 10),
                          FLAGS_LIST = list(ADAPTIVE = TRUE, THIN = TRUE)) {    
  
  #NOTE:
  #i - 1 = n (Simon's paper)
  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  
  #MCMC PARAMS + VECTORS
  time = length(dataX); n_mcmc = mcmc_inputs$n_mcmc;
  dim = mcmc_inputs$dim; count_accept = 0; count_accept_da = 0
  vec_min = rep(0, mcmc_inputs$dim) 
  
  #THINNING FACTOR
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1; mcmc_vec_size = n_mcmc
  }
  
  #MODEL PARAMS
  eta = dataX
  nu_params_matrix = matrix(NA, mcmc_vec_size, dim);   #Changed from 0 to NA (As should be overwriting all cases)
  nu_params_matrix[1,] <- mcmc_inputs$mod_start_points;
  nu_params = nu_params_matrix[1,] #2x1 #as.matrix
  eta_matrix = matrix(NA, mcmc_vec_size, num_days);  #eta_mean_vec <- vector('numeric', mcmc_vec_size);
  log_like_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec[1] <- LOG_LIKELIHOOD_NU(dataX, nu_params, eta);  log_like = log_like_vec[1]
  
  #ADAPTIVE SHAPING PARAMS + VECTORS
  lambda_vec <- vector('numeric', mcmc_vec_size); lambda_vec[1] <- 1
  c_star = (2.38^2)/dim; termX = mcmc_inputs$v0 + dim
  delta = 1/(mcmc_inputs$alpha_star*(1 - mcmc_inputs$alpha_star))
  print(paste0('delta = ', delta))
  sigma_i = diag(dim); lambda_i = 1
  
  #MCMC
  for(i in 2:n_mcmc) {
    #print(paste0('i = ', i))
    #PRINT PROGRESS
    if(i%%(n_mcmc/50) == 0) print(paste0('i = ', i))
    
    #SIGMA ITERATION NO.1
    if (i == 2){
      x_bar = 0.5*(nu_params + nu_params_matrix[1,])
      sigma_i = (1/(termX + 3))*(tcrossprod(nu_params_matrix[1,]) + tcrossprod(nu_params) -
                                   2*tcrossprod(x_bar) + (termX + 1)*sigma_i) #CHANGE TO USE FUNCTIONS
    }
    #print(paste0('sigma_i = ', sigma_i))
    #print(paste0('lambda_i = ', lambda_i))
    #print(paste0('c_star = ', c_star))
    #PROPOSAL
    
    nu_params_dash = c(nu_params + mvrnorm(1, mu = rep(0, dim), Sigma = lambda_i*c_star*sigma_i)) #Vectorise using c()
    
    #ONLY KEEP POSTIVE
    if (min(nu_params_dash - vec_min) < 0){ #REJECT IF C - 1 <0 I.E C < 1 #Model specific. #
      
      #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
      if (i%%thinning_factor == 0) {
        nu_params_matrix[i/thinning_factor,] = nu_params  # i%/%
        log_like_vec[i/thinning_factor] <- log_like
        lambda_vec[i/thinning_factor] <- lambda_i #Taking role of sigma, overall scaling constant. Sigma becomes estimate of the covariance matrix of the posterior
      }
      
      next
    }
    
    #LOG LIKE
    logl_new = LOG_LIKELIHOOD_NU(dataX, nu_params_dash, eta)
    #print(paste0('logl_new', logl_new))
    
    #ACCEPTANCE RATIO
    log_accept_ratio = logl_new - log_like #PRIORS?

    #METROPOLIS ACCEPTANCE STEP
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      nu_params <- nu_params_dash
      count_accept = count_accept + 1
      log_like = logl_new
    }
    
    #SIGMA - ADAPTIVE SHAPING
    if (i > 2) {
      #SIGMA - ADAPTIVE SHAPING
      xbar_prev = x_bar
      x_bar = (i-1)/i*xbar_prev + (1/i)*nu_params
      sigma_i = (1/(i + termX + 1))*( (i + termX)*sigma_i +tcrossprod(nu_params)
                                      + (i-1)*tcrossprod(xbar_prev)
                                      -i*tcrossprod(x_bar))
    }
    
    #LAMBDA - ADAPTIVE SCALING
    accept_prob = min(1, exp(log_accept_ratio))
    lambda_i =  lambda_i*exp(delta/i*(accept_prob - mcmc_inputs$alpha_star))
    #print(paste0('lambda_i = ', lambda_i))
    #print(paste0('accept_prob = ', accept_prob))
    
    #************************************
    #DATA AUGMENTATION
    #************************************
    for(t in 1:time){
      
      v = rep(0, length(eta)); v[t] = 1
      
      #METROPOLIS STEP 
      eta_dash = abs(eta + rnorm(1,0,1)*v) #normalise the t_th element of eta #or variance = x[t]
      
      #LOG LIKELIHOOD
      logl_new = LOG_LIKELIHOOD_NU(dataX, nu_params_dash, eta)
      log_accept_ratio = logl_new - log_like
      
      #METROPOLIS ACCEPTANCE STEP
      if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
        
        #ACCEPT
        eta <- eta_dash
        log_like <- logl_new
        count_accept_da = count_accept_da + 1
      }
    }
    
    #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0) {
      nu_params_matrix[i/thinning_factor,] = nu_params
      log_like_vec[i/thinning_factor] <- log_like
      lambda_vec[i/thinning_factor] <- lambda_i #Taking role of sigma, overall scaling constant. Sigma becomes estimate of the covariance matrix of the posterior
      eta_matrix[i/thinning_factor, ] <- eta 
    }
    
  } #END FOR LOOP
  
  #Final stats
  accept_rate = 100*count_accept/(n_mcmc-1)
  accept_rate_da = 100*count_accept_da/((n_mcmc-1)*time)

  #Return a, acceptance rate
  return(list(nu_params_matrix = nu_params_matrix, eta_matrix = eta_matrix,
              log_like_vec = log_like_vec, lambda_vec = lambda_vec,
              accept_rate = accept_rate, accept_rate_da = accept_rate_da))
} 

#********************************************************
#1. INDIVIDUAL R0 MCMC                            
#********************************************************
MCMC_MODEL_NU <- function(dataX,
                          mcmc_inputs = list(n_mcmc = 100,
                                             mod_start_points = c(1.2, 0.16),  #priors_list = list(alpha_prior = c(1, 0), k_prior = c()),
                                             thinning_factor = 3, dim = 2, seed_count = 1),
                          FLAGS_LIST = list(ADAPTIVE = TRUE, THIN = FALSE)) {    
  
  #NOTE:
  #i - 1 = n (Simon's paper)
  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  
  #MCMC PARAMS + VECTORS
  time = length(dataX); n_mcmc = mcmc_inputs$n_mcmc;
  count_accept = 0; count_accept_da = 0; dim = mcmc_inputs$dim
  vec_min = rep(0, dim) 
  
  #THINNING FACTOR
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1; mcmc_vec_size = n_mcmc
  }
  
  #MODEL PARAMS
  eta = dataX
  #SIGMA
  sigmaX = diag(2); 
  sigmaX[1,1] = 0.5*mcmc_inputs$mod_start_points[1];  sigmaX[2,2] = 0.5*mcmc_inputs$mod_start_points[2]
  #NU
  nu_params_matrix = matrix(NA, mcmc_vec_size, dim);   #Changed from 0 to NA (As should be overwriting all cases)
  nu_params_matrix[1,] <- mcmc_inputs$mod_start_points;
  nu_params = nu_params_matrix[1,] #2x1 #as.matrix
  eta_matrix = matrix(NA, mcmc_vec_size, num_days);  #eta_mean_vec <- vector('numeric', mcmc_vec_size);
  log_like_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec[1] <- LOG_LIKELIHOOD_NU(dataX, nu_params, eta);  log_like = log_like_vec[1]
  
  #MCMC
  for(i in 2:n_mcmc) {

    if(i%%(n_mcmc/50) == 0) print(paste0('i = ', i))
    
    #NU dash
    nu_params_dash = c(nu_params + mvrnorm(1, mu = rep(0, dim), Sigma = sigmaX)) #Vectorise using c()
    
    #ONLY KEEP POSTIVE
    if (min(nu_params_dash - mcmc_inputs$vec_min) < 0){ #REJECT IF C - 1 <0 I.E C < 1 #Model specific. #
      
      #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
      if (i%%thinning_factor == 0) {
        nu_params_matrix[i/thinning_factor,] = nu_params  # i%/%
        log_like_vec[i/thinning_factor] <- log_like
      }
      
      next
    }
    
    #LOG LIKE
    logl_new = LOG_LIKELIHOOD_NU(dataX, nu_params_dash, eta)
    
    #ACCEPTANCE RATIO
    log_accept_ratio = logl_new - log_like #PRIORS?
    
    #METROPOLIS ACCEPTANCE STEP
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      nu_params <- nu_params_dash
      count_accept = count_accept + 1
      log_like = logl_new
    }
    
    #************************************
    #DATA AUGMENTATION
    #************************************
    for(t in 1:time){
      
      v = rep(0, length(eta)); v[t] = 1
      
      #METROPOLIS STEP 
      eta_dash = abs(eta + rnorm(1,0,1)*v) #normalise the t_th element of eta #or variance = x[t]
      
      #LOG LIKELIHOOD
      logl_new = LOG_LIKELIHOOD_NU(dataX, nu_params_dash, eta)
      log_accept_ratio = logl_new - log_like
      
      #METROPOLIS ACCEPTANCE STEP
      if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
        
        #ACCEPT
        eta <- eta_dash
        log_like <- logl_new
        count_accept_da = count_accept_da + 1
      }
    }
    
    #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0) {
      nu_params_matrix[i/thinning_factor,] = nu_params
      log_like_vec[i/thinning_factor] <- log_like
      eta_matrix[i/thinning_factor, ] <- eta 
    }
    
    if (i == 20000) {
      saveRDS(nu_params_matrix, file = paste0(OUTER_FOLDER, 'nu_params_matrix_', seed_count, '.rds' ))
      saveRDS(eta_matrix, file = paste0(OUTER_FOLDER, 'eta_matrix_', seed_count, '.rds' ))
      saveRDS(log_like_vec, file = paste0(OUTER_FOLDER, 'log_like_vec_', seed_count, '.rds' ))
    }
  } #END FOR LOOP
  
  #Final stats
  accept_rate = 100*count_accept/(n_mcmc-1)
  accept_rate_da = 100*count_accept_da/((n_mcmc-1)*time)
  
  #Return a, acceptance rate
  return(list(nu_params_matrix = nu_params_matrix, 
              eta_matrix = eta_matrix, log_like_vec = log_like_vec,
              accept_rate = accept_rate, accept_rate_da = accept_rate_da))
} 

#**************************
#APPLY MODELS
#**************************

#APPLY I
mcmc_nu = MCMC_ADAPTIVE_MODEL_NU(canadaX)

#SAVE
iter = 'II'
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_individual_nu/"
ifelse(!dir.exists(file.path(OUTER_FOLDER)), dir.create(file.path(OUTER_FOLDER), recursive = TRUE), FALSE)
saveRDS(mcmc_nu, file = paste0(OUTER_FOLDER, 'mcmc_nu_', iter, '.rds' ))

#RESULTS
#OUTPUT
burn_in = 100
alpha_mcmc = mcmc_nu2$nu_params_matrix[,1]
plot.ts(alpha_mcmc)
#BURN-IN
alpha_mcmcB = alpha_mcmc[burn_in:length(alpha_mcmc)]
plot.ts(alpha_mcmcB)

#k
k_mcmc = mcmc_nu2$nu_params_matrix[,2]
plot.ts(k_mcmc)
#BURN-IN
k_mcmcB = k_mcmc[burn_in:length(k_mcmc)]
plot.ts(k_mcmcB)

#ETA
eta_matrix = mcmc_nu$eta_matrix

#LOGLIKE
log_like_vec =  mcmc_nu$log_like_vec
plot.ts(log_like_vec)

#PLOT
mcmc_nu$time_elap = "25 mins"
PLOT_NU_MCMC_GRID(canadaX, mcmc_nu)

#**********************
#APPLY II

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_nu2 = MCMC_MODEL_NU(canadaX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_nu2$time_elap = time_elap

#PLOT
PLOT_NU_MCMC_GRID(canadaX, mcmc_nu2)
