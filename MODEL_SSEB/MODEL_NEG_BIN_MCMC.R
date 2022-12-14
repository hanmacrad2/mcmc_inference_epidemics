#MODEL NEGATIVE BINOMIAL II

#****************************
#SIMULATION
SIMULATE_NU = function(num_days = 110, alphaX = 1.2, k = 0.16,
                               shape_gamma = 6, scale_gamma = 1) {
  
  'Simulate from the Negative Binomial model'
  
  #INTIALISE VECTORS
  x = vector('numeric', num_days); x[1] = 2
  
  #INFECTIOUSNESS (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  vec_nu_t = vector('numeric', num_days); 
  
  #DAYS OF THE EPIDEMIC
  for (t in 2:num_days) {
    
    #NU: INDIVIDUAL R0; GAMMA(ALPHA, K)
    vec_nu_t[t-1] <- rgamma(1, shape = x[t-1]*k, scale = alphaX/k)
    
    #OFFSPRING; POISSON()
    infectivity = rev(prob_infect[1:(t-1)]) #x[1:(t-1)]*rev(prob_infect[1:(t-1)])
    total_rate = sum(vec_nu_t*infectivity)
    #print(total_rate)
    x[t] = rpois(1, total_rate)
    
  }
  return(x)
}

#****************************
#LOG LIKELIHOOD

#*******************************************************************************
#LOG LIKELIHOOD
LOG_LIKELIHOOD_NU <- function(x, alphaX, k){
  
  #Params
  num_days = length(x)
  shape_gamma = 6; scale_gamma = 1
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - 
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  loglike = 0
  
  for (t in 2:num_days) {
    
    #NU: INDIVIDUAL R0; GAMMA(ALPHA, K)
    vec_nu_t[t-1] <- rgamma(1, shape = x[t-1]*k, scale = alphaX/k)
    
    #OFFSPRING; POISSON()
    infectivity = rev(prob_infect[1:(t-1)]) #x[1:(t-1)]*rev(prob_infect[1:(t-1)])
    total_rate = sum(vec_nu_t*infectivity)

    x[t] = rpois(1, total_rate)
    
    loglike = loglike + x[t]*log(total_rate) - total_rate - lfactorial(x[t]) #Need to include normalizing constant 
    
  }
  
  return(loglike)
  
}

#*******************************************************************************
# #DATA AUGMENTATION      (W/ DATA AUGMENTATION)
#************************************************************************

MODEL_NU_MCMC <- function(data,
                              mcmc_inputs = list(n_mcmc = 500000,
                                                 mod_start_points = list(m1 = 0.72, m2 = 0.0038), alpha_star = 0.4,
                                                 thinning_factor = 10),
                              priors_list = list(alpha_prior = c(1, 0), k_prior = c()),
                              FLAGS_LIST = list(ADAPTIVE = TRUE, THIN = TRUE)) {
  
  'Returns MCMC samples of MODEL_NU params 
  w/ acceptance rates.
  INCLUDES; ADAPTATION, DATA AUGMENTATION'
  
  'Priors
  p(a) = exp(rate) = rate*exp(-rate*x). log(r*exp(-r*x)) = log(r) - rx
      -> E.g exp(1) = 1*exp(-1*a) = exp(-a). log(exp(-a)) = - a
  p(b) = exp(1) or p(b) = g(shape, scale), for e.g g(3, 2)
  p(c) = exp(1) + 1 = 1 + exp(-c) = exp(c - 1)'
  
  
  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  time = length(data[[1]]); n_mcmc = mcmc_inputs$n_mcmc;
  print(paste0('num mcmc iters = ', n_mcmc))
  
  #THINNING FACTOR
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1
    mcmc_vec_size = n_mcmc
  }
  
  #INITIALISE MCMC VECTORS
  alpha_vec <- vector('numeric', mcmc_vec_size); k_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec <- vector('numeric', mcmc_vec_size);
  
  #INITIALISE MCMC[1]
  alpha_vec[1] <- mcmc_inputs$mod_start_points$m1; k_vec[1] <- mcmc_inputs$mod_start_points$m2
  log_like_vec[1] <- LOG_LIKELIHOOD_NU(data, alpha_vec[1], k_vec[1])
  
  #INITIALISE RUNNING PARAMS
  alpha = a_vec[1];  k = b_vec[1]; log_like = log_like_vec[1]
  
  #SIGMA
  sigma1 =  0.4*mcmc_inputs$mod_start_points$m1;  sigma2 = 0.4*mcmc_inputs$mod_start_points$m2
  
  #SIGMA; INITIALISE FOR ADAPTIVE MCMC
  if (FLAGS_LIST$ADAPTIVE){
    
    #SIGMA
    sigma1_vec <- vector('numeric', mcmc_vec_size); sigma2_vec <- vector('numeric', mcmc_vec_size)
    #SIGMA; INITIALISE FIRST ELEMENT
    sigma1_vec[1] =  sigma1; sigma2_vec[1] =  sigma2
    #SIGMA; List of sigma vectors for each iteration of the MCMC algorithm
    sigma = list(sigma1_vec = sigma1_vec, sigma2_vec = sigma2_vec)
    
    #Other adaptive parameters
    delta = 1/(mcmc_inputs$alpha_star*(1-mcmc_inputs$alpha_star))
    
  } else {
    
    #SIGMA; List of sigma vectors for each iteration of the MCMC algorithm
    sigma = list(sigma1 <- sigma1, sigma2 <- sigma2)
  }
  
  #INITIALISE: ACCEPTANCE COUNTS
  list_accept_counts = list(count_accept1 = 0, count_accept2 = 0)
  
  #DATA AUG OUTPUT
  #mat_count_da = matrix(0, mcmc_vec_size, time) #i x t
  #non_ss = matrix(0, mcmc_vec_size, time) #USE THINNING FACTOR
  #ss = matrix(0, mcmc_vec_size, time) #USE THINNING FACTOR
  
  #******************************
  #MCMC CHAIN
  #******************************
  for(i in 2:n_mcmc) {
    
    #****************************************************** 
    #MCMC PARAMS
    
    #FOR EACH S_T
    for(t in 1:time){
      
      #Copy of data 
      #data_dash = data
      
      #STOCHASTIC PROPOSAL for s
      if (runif(1) < 0.5) {
        param_dash = param + 1
      } else {
        st_dash = param - 1
      }
      
      #CHECK
      if (param_dash < 0) {
        param_dash = 0
      }
      
      #ACCEPTANCE PROBABILITY
      # data_dash[[2]][t] = st_dash #s_t = st_dash
      # data_dash[[1]][t] =  data[[1]][t] + data[[2]][t] - st_dash #n_t = x_t - s_t
      # 
      # if (data_dash[[1]][t] < 0){
      #   data_dash[[1]][t] = 0
      # }
      # 
      # #CRITERIA FOR S_T & N_T
      # if((data_dash[[2]][t] < 0) || (data_dash[[1]][t] < 0)){
      #   #Store
      #   non_ss[i/thinning_factor, t] = data[[1]][t]
      #   ss[i/thinning_factor, t] = data[[2]][t]
      #   next
      # }
      
      #LOG LIKELIHOOD
      logl_new = LOG_LIKELIHOOD_NU(data_dash, alpha, k)
      log_accept_ratio = logl_new - log_like
      
      #METROPOLIS ACCEPTANCE STEP
      if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
        
        #ACCEPT
        data <- data_dash
        log_like <- logl_new
        #mat_count_da[i, t] = mat_count_da[i, t] + 1
        list_accept_counts$count_accept1 = list_accept_counts$count_accept1 + 1
      }
      
      #Store
      # non_ss[i/thinning_factor, t] = data[[1]][t] #TAKE MEAN ACROSS MCMC DIMENSION (PLOT 0 > 50)
      # ss[i/thinning_factor, t] = data[[2]][t]
    }
  }
  
  #Loglikelihood Check (Passing - no error)
  #if (!(log_like && LOG_LIKE_SSI(data, a, b, c))) print('loglike doesnt exist')
  #else (log_like!=LOG_LIKE_SSI(data, a, b, c)) print(paste0('ERROR! logl diff = ', log_like - LOG_LIKE_SSI(data, a, b, c)))
  
  #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
  if (i%%thinning_factor == 0) {
    #print(paste0('i = ', i))
    i_thin = i/thinning_factor
    alpha_vec[i_thin] <- alpha; k_vec[i_thin] <- k
    log_like_vec[i_thin] <- log_like 
    sigma$sigma1[i_thin] = sigma1; sigma$sigma2[i_thin] = sigma2
  }

#Final stats
#accept_rate1 = 100*list_accept_counts$count_accept1/(n_mcmc-1)
#accept_rate2 = 100*list_accept_counts$count_accept2/(n_mcmc-1) #(list_accept_counts$count_accept2 + list_reject_counts$count_accept2)
accept_rate1 = 100*list_accept_counts$count_accept1/((n_mcmc-1)*time) #i x t

#Acceptance rates
list_accept_rates = list(accept_rate1 = accept_rate1,
                         accept_rate2 = accept_rate2)
print(list_accept_rates)

#Return a, acceptance rate
return(list(alpha_vec = alpha_vec, k_vec = k_vec,
            log_like_vec = log_like_vec, sigma = sigma,
            list_accept_rates = list_accept_rates))
}
