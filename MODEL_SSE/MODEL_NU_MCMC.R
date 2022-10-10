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
  
  vec_sum_nu_t = vector('numeric', num_days); 
  
  #DAYS OF THE EPIDEMIC
  for (t in 2:num_days) {
    
    #NU: INDIVIDUAL R0; GAMMA(ALPHA, K)
    vec_sum_nu_t[t-1] <- rgamma(1, shape = x[t-1]*k, scale = alphaX/k)
    
    #OFFSPRING; POISSON()
    infectivity = rev(prob_infect[1:(t-1)]) #x[1:(t-1)]*rev(prob_infect[1:(t-1)])
    total_rate = sum(vec_sum_nu_t*infectivity)
    #print(total_rate)
    x[t] = rpois(1, total_rate)
    
  }
  return(x)
}

#****************************
#LOG LIKELIHOOD

#*******************************************************************************
#LOG LIKELIHOOD
LOG_LIKELIHOOD_NU <- function(x, alphaX, k, eta){ #eta - a vector of length x. eta[1] = infectivity of x[1]
  
  #Params
  num_days = length(x)
  shape_gamma = 6; scale_gamma = 1
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - 
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  loglike = 0
  
  for (t in 2:num_days) {
    
    #NU: INDIVIDUAL R0; GAMMA(ALPHA, K)
    #vec_nu_t[t-1] <- rgamma(1, shape = x[t-1]*k, scale = alphaX/k)
    
    #OFFSPRING; POISSON()
    infectivity = rev(prob_infect[1:(t-1)]) 
    total_rate = sum(eta[1:(t-1)]*infectivity)
    
    #dgamma #EXPLAIN DGAMMA??
    loglike = loglike + dgamma(eta[t-1], shape = x[t-1]*k, scale = alphaX/k, log = TRUE)
    loglike = loglike + x[t]*log(total_rate) - total_rate - lfactorial(x[t]) #Need to include normalizing constant 
    
  }
  
  return(loglike)
  
}

LOG_LIKELIHOOD_NU(c(1,2), alphaX = 1.2, k = 100000, eta = c(1,2))

#*******************************************************************************
# #DATA AUGMENTATION      (W/ DATA AUGMENTATION)
#************************************************************************
MODEL_NU_MCMC <- function(data,
                              mcmc_inputs = list(n_mcmc = 1000,
                                                 mod_start_points = list(m1 = 0.72, m2 = 0.0038), alpha_star = 0.4,
                                                 thinning_factor = 10),
                              priors_list = list(alpha_prior = c(1, 0), k_prior = c()),
                              FLAGS_LIST = list(ADAPTIVE = TRUE, THIN = TRUE)) {
  
  'Returns MCMC samples of MODEL_NU params 
  w/ acceptance rates.
  INCLUDES; ADAPTATION, DATA AUGMENTATION'
  
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
  eta = x #Initialise as data 
  
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
  
  #******************************
  #MCMC CHAIN
  #******************************
  for(i in 2:n_mcmc) {
    
    #alpha (metropolis + k; jointly) multivariate normal for proposal
    #START PRETEND KNOW ALPHA + K 
    
    #****************************************************** 
    #MCMC PARAMS
    
    #FOR EACH S_T
    for(t in 1:time){
      
      v = rep(0, length(eta)); v[t] = 1
      #Metropolis 
      eta_dash = abs(eta + rnorm(1,0,1)*v) #normalise the t_th element of eta #or variance = x[t]
     
      #LOG LIKELIHOOD
      logl_new = LOG_LIKELIHOOD_NU(data_dash, alpha, k, eta_dash)
      log_accept_ratio = logl_new - log_like
      
      #METROPOLIS ACCEPTANCE STEP
      if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
        
        #ACCEPT
        eta <- eta_dash
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

#PARAMS CORRECT?
n = 10000
zz = rgamma(n, shape = x[t-1]*k, scale = alphaX/k)
