#NOTE NO REFLECTION, NO TRANSFORMS, MORE INTELLIGENT ADAPTATION
MCMC_ADAPTIVE_MODEL_NU <- function(dataX, OUTER_FOLDER, seed_count,
                                   mcmc_inputs = list(n_mcmc = 20000, #000,
                                                      mod_start_points = c(1.2, 0.16),
                                                      dim = 2, target_acceptance_rate = 0.4, v0 = 100,  #priors_list = list(alpha_prior = c(1, 0), k_prior = c()),
                                                      thinning_factor = 2),
                                   FLAGS_LIST = list(ADAPTIVE = TRUE, THIN = TRUE)) {    
  
  #NOTE:
  #i - 1 = n (Simon's paper)
  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  
  #MCMC PARAMS + VECTORS
  num_days = length(dataX); n_mcmc = mcmc_inputs$n_mcmc;
  dim = mcmc_inputs$dim; count_accept = 0; count_accept_da = 0
  vec_min = rep(0, mcmc_inputs$dim);
  
  #DIRECTORY - SAVING
  ifelse(!dir.exists(file.path(OUTER_FOLDER)), dir.create(file.path(OUTER_FOLDER), recursive = TRUE), FALSE)
  
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
  delta = 1/(mcmc_inputs$target_acceptance_rate*(1 - mcmc_inputs$target_acceptance_rate))
  x_bar = 0.5*(nu_params + nu_params_matrix[1,])
  sigma_i = diag(dim); lambda_i = 1
  sigma_i = (1/(termX + 3))*(tcrossprod(nu_params_matrix[1,]) + tcrossprod(nu_params) -
                               2*tcrossprod(x_bar) + (termX + 1)*sigma_i) #CHANGE TO USE FUNCTIONS
  
  #ETA (ONE DIMENSION - ROBINS MUNROE) 
  sigma_eta =  0.5*rep(1, num_days)
  sigma_eta_matrix = matrix(0, mcmc_vec_size, num_days); sigma_eta_matrix[1,] =  sigma_eta;
  delta_eta = 1/(mcmc_inputs$target_acceptance_rate*(1-mcmc_inputs$target_acceptance_rate))
  
  #MCMC
  for(i in 2:n_mcmc) {
    
    if(i%%(n_mcmc/100) == 0) print(paste0('i = ', i))
    
    #PROPOSAL
    nu_params_dash = c(nu_params + mvrnorm(1, mu = rep(0, dim), Sigma = lambda_i*c_star*sigma_i)) 
    
    #POSTIVE ONLY
    if (min(nu_params_dash - vec_min) >= 0){ 
    
    #LOG LIKELIHOOD
    logl_new = LOG_LIKELIHOOD_NU(dataX, nu_params_dash, eta)
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
    
    } else {
      
      accept_prob = 0
    }
    
    lambda_i =  lambda_i*exp(delta/i*(accept_prob - mcmc_inputs$target_acceptance_rate))
    
    #************************************
    #DATA AUGMENTATION
    #************************************
    for(t in 1:num_days){
      
      v = rep(0, length(eta)); v[t] = 1
      
      #METROPOLIS STEP 
      eta_dash = abs(eta + rnorm(1,0,sigma_eta[t])*v) #normalise the t_th element of eta #or variance = x[t]
      
      #LOG LIKELIHOOD
      logl_new = LOG_LIKELIHOOD_NU(dataX, nu_params, eta_dash)
      log_accept_ratio = logl_new - log_like
      
      #METROPOLIS ACCEPTANCE STEP
      if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
        
        #ACCEPT
        eta <- eta_dash
        log_like <- logl_new
        count_accept_da = count_accept_da + 1
      }
      
      accept_prob = min(1, exp(log_accept_ratio)) #Acceptance PROB = MIN(1, EXP(ACCPET_PROB))
      sigma_eta[t] =  sigma_eta[t]*exp(delta_eta/(1+i)*(accept_prob - mcmc_inputs$target_acceptance_rate))
    }
    
    #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0) {
      nu_params_matrix[i/thinning_factor,] = nu_params
      log_like_vec[i/thinning_factor] <- log_like
      lambda_vec[i/thinning_factor] <- lambda_i #Taking role of sigma, overall scaling constant. Sigma becomes estimate of the covariance matrix of the posterior
      eta_matrix[i/thinning_factor, ] <- eta 
      sigma_eta_matrix[i/thinning_factor, ] = sigma_eta
    }
    
  } #END FOR LOOP
  
  #SAVE
  saveRDS(nu_params_matrix, file = paste0(OUTER_FOLDER, 'nu_params_matrix_', seed_count, '.rds' ))
  saveRDS(eta_matrix, file = paste0(OUTER_FOLDER, 'eta_matrix_', seed_count, '.rds' ))
  saveRDS(log_like_vec, file = paste0(OUTER_FOLDER, 'log_like_vec_', seed_count, '.rds' ))
  
  #Final stats
  accept_rate = 100*count_accept/(n_mcmc-1)
  accept_rate_da = 100*count_accept_da/((n_mcmc-1)*num_days)
  
  #Return a, acceptance rate
  return(list(nu_params_matrix = nu_params_matrix, eta_matrix = eta_matrix,
              sigma_eta_matrix = sigma_eta_matrix,
              log_like_vec = log_like_vec, lambda_vec = lambda_vec, 
              accept_rate = accept_rate, accept_rate_da = accept_rate_da))
} 
