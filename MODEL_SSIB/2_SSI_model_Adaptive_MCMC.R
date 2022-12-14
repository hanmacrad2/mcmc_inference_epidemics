#****************************************************************
#1. SSI MODEL MCMC + DATA AUGMENTATION
#****************************************************************
#*
#*******************************
#1. SSI MODEL - LOG LIKELIHOOD
#*******************************
LOG_LIKE_SSI <- function(sim_data, aX, bX, cX){
  
  #Data
  n = sim_data[[1]]; s = sim_data[[2]]
  
  #Params
  num_days = length(n)
  shape_gamma = 6; scale_gamma = 1
  logl = 0
  
  #INFECTIOUSNESS  - Difference of 2 GAMMA distributions. Discretized 
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  for (t in 1:num_days) { #*1 or 2
    
    #INFECTIOUS PRESSURE - SUM OF ALL INDIVIDUALS INFECTIOUSNESS 
    lambda_t = sum((n[1:(t-1)] + cX*s[1:(t-1)])*rev(prob_infect[1:(t-1)]))
    
    #LOG-LIKELIHOOD 
    logl = logl - lambda_t*(aX + bX) + n[t]*(log(aX) + log(lambda_t)) + s[t]*(log(bX) + log(lambda_t))  + 2*log(1) - lfactorial(n[t]) - lfactorial(s[t])
  }
  
  logl
}

#CHECK ############
# LOG_LIKE_SSI <- function(sim_data, aX, bX, cX){
#   
#   #Data
#   n = sim_data[[1]]; s = sim_data[[2]]
#   
#   #Params
#   num_days = length(n)
#   shape_gamma = 6; scale_gamma = 1
#   logl = 0;
#   
#   #INFECTIOUSNESS  - Difference of 2 GAMMA distributions. Discretized 
#   prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
#     pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
#   
#   for (t in 1:num_days) { #*1 or 2
#     
#     #INFECTIOUS PRESSURE - SUM OF ALL INDIVIDUALS INFECTIOUSNESS 
#     lambda_t = sum((n[1:(t-1)] + cX*s[1:(t-1)])*rev(prob_infect[1:(t-1)]))
#     
#     #LOG-LIKELIHOOD 
#     log_like_add = - lambda_t*(aX + bX) + n[t]*(log(aX) + log(lambda_t)) + s[t]*(log(bX) + log(lambda_t))  + 2*log(1) - lfactorial(n[t]) - lfactorial(s[t])
#     
#     #LOG LIKE
#     logl = logl + log_like_add
#     
#     # if (is.na(log_like_add)){
#     #   # print(paste0('WARNING NaN,  log_like_add = ', log_like_add))
#     #   # print(paste0('lambda_t = ', lambda_t))
#     #   # print(paste0('n[t] = ', n[t]))
#     #   # print(paste0('s[t] = ', s[t]))
#     #   next
#     # } else if (is.infinite(log_like_add)) {
#     #   # print(paste0('WARNING Inf,  log_like_add = ', log_like_add))
#     #   # print(paste0('lambda_t = ', lambda_t))
#     #   # print(paste0('n[t] = ', n[t]))
#     #   # print(paste0('s[t] = ', s[t]))
#     #   
#     # } else {
#     #   logl = logl + log_like_add
#     # } 
#       #logl = logl - lambda_t*(aX + bX) + n[t]*(log(aX) + log(lambda_t)) + s[t]*(log(bX) + log(lambda_t))  + 2*log(1) - lfactorial(n[t]) - lfactorial(s[t])
#     
#     }
#     
#   #   if (is.na(logl)){
#   #     logl = 0
#   #     print(paste0('1. logl  = ', logl))
#   #     # print('WARNING: LOGLIKE IS NA')
#   #     # print(paste0('lambda_t = ', lambda_t))
#   #     # print(paste0('n[t] = ', n[t]))
#   #     # print(paste0('s[t] = ', s[t]))
#   #     # print('****************')
#   #     # break
#   #   }
#   # 
#   #print(paste0('2. logl  = ', logl))
#   
#   logl
# }

#************************************************************************
#1. SSI MCMC                              (W/ DATA AUGMENTATION OPTION)
#************************************************************************
MCMC_SSI_ADAPTIVE <- function(data,
                     mcmc_inputs = list(n_mcmc = n_mcmc,
                                        mod_start_points = mod_start_points, alpha_star = 0.4), #THINNING FACTOR, burn_in  
                     priors_list = list(a_prior = c(1, 0), b_prior = c(10, 2/100), b_prior_exp = c(1,0),  #c(0.1,0),#10, 1/100
                                        c_prior = c(10, 1), c_prior_exp = c(0.1,0)),
                     FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                       PRIOR = TRUE,
                                       B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE  )) {
                                       
  'Returns MCMC samples of SSI model parameters (a, b, c, r0 = a + b*c) 
  w/ acceptance rates.
  INCLUDES; DATA AUGMENTATION, B-C transform' 
  
  'Priors
  p(a) = exp(rate) = rate*exp(-rate*x). log(r*exp(-r*x)) = log(r) - rx
      -> E.g exp(1) = 1*exp(-1*a) = exp(-a). log(exp(-a)) = - a
  p(b) = exp(1) or p(b) = g(shape, scale), for e.g g(3, 2)
  p(c) = exp(1) + 1 = 1 + exp(-c) = exp(c - 1)'
  
  
  #**********************************************
  #INITIALISE PARAMS 
  #**********************************************
  time = length(data[[1]]); n_mcmc = mcmc_inputs$n_mcmc;
  delta = 1/(mcmc_inputs$alpha_star*(1-mcmc_inputs$alpha_star))
  
  #INITIALISE MCMC VECTORS
  a_vec <- vector('numeric', n_mcmc); b_vec <- vector('numeric', n_mcmc)
  c_vec <- vector('numeric', n_mcmc); r0_vec <- vector('numeric', n_mcmc)
  log_like_vec <- vector('numeric', n_mcmc);
  
  #INITIALISE MCMC[1]    
  a_vec[1] <- mcmc_inputs$mod_start_points$m1; b_vec[1] <- mcmc_inputs$mod_start_points$m2
  c_vec[1] <- mcmc_inputs$mod_start_points$m3; r0_vec[1] <- mcmc_inputs$mod_start_points$true_r0;
  log_like_vec[1] <- LOG_LIKE_SSI(data, a_vec[1], b_vec[1], c_vec[1])
  
  #INITIALISE RUNNING PARAMS
  a = mcmc_inputs$mod_start_points$m1; b =  mcmc_inputs$mod_start_points$m2; 
  c = mcmc_inputs$mod_start_points$m3; log_like = log_like_vec[1]
  
  #SIGMA INITIALISE FOR ADAPTIVE MCMC
  sigma = list(sigma1 <- vector('numeric', n_mcmc), sigma2 <- vector('numeric', n_mcmc), sigma3 <- vector('numeric', n_mcmc),
               sigma4 <- vector('numeric', n_mcmc), sigma5 <- vector('numeric', n_mcmc))
  
  #SIGMA INITIALISE FIRST ELEMENT 
  sigma$sigma1[1] =  0.4*mcmc_inputs$mod_start_points$m1;  sigma$sigma2[1] = 0.3*mcmc_inputs$mod_start_points$m2
  sigma$sigma3[1] = 0.5*mcmc_inputs$mod_start_points$m3; sigma$sigma4[1] =   0.85*mcmc_inputs$mod_start_points$m3
  sigma$sigma5[1] =   0.85*mcmc_inputs$mod_start_points$m3
  
  sigma1 = sigma$sigma1[1]; sigma2 = sigma$sigma2[1]; sigma3 = sigma$sigma3[1]; 
  sigma4 = sigma$sigma4[1]; sigma5 = sigma$sigma5[1]
  
  print(paste0('sigma1 = ', sigma1))
  
  #INITIALISE: ACCEPTANCE COUNTS 
  list_accept_counts = list(count_accept1 = 0, count_accept2 = 0, count_accept3 = 0,
                            count_accept4 = 0, count_accept5 = 0, count_accept6 = 0)
  
  #DATA AUG OUTPUT
  mat_count_da = matrix(0, n_mcmc, time) #i x t
  non_ss = matrix(0, n_mcmc, time) #USE THINNING FACTOR
  ss = matrix(0, n_mcmc, time) #USE THINNING FACTOR
  
  #******************************
  #MCMC CHAIN
  #******************************
  for(i in 2:n_mcmc) {

    #****************************************************** s
    #a
    a_dash <- a + rnorm(1, sd = sigma1) 
      
    if(a_dash < 0){
      a_dash = abs(a_dash)
    }
    
    #log a
    logl_new = LOG_LIKE_SSI(data, a_dash, b, c)
    log_accept_ratio = logl_new - log_like  #+ prior1 - prior
    #Priors
    if (FLAGS_LIST$PRIOR){
      log_accept_ratio = log_accept_ratio - a_dash + a #*Actually this is the Acceptance RATIO. ACCEPTANCE PROB = MIN(1, EXP(ACCPET_PROB))
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      a <- a_dash
      list_accept_counts$count_accept1 = list_accept_counts$count_accept1 + 1
      log_like = logl_new
    } 
    
    #Sigma (Adpative)
    accept_prob = min(1, exp(log_accept_ratio))  
    sigma1 =  sigma1*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    
    #************************************************************************ Only if (b > 0){ ?
    #b  
    b_dash <- b + rnorm(1, sd = sigma2) 
    if(b_dash < 0){
      b_dash = abs(b_dash)
    }
    
    #loglikelihood
    logl_new = LOG_LIKE_SSI(data, a, b_dash, c)
    log_accept_ratio = logl_new - log_like
    
    #Priors
    if (FLAGS_LIST$B_PRIOR_GAMMA){
      log_accept_ratio = log_accept_ratio +
        dgamma(b_dash, shape = priors_list$b_prior[1], scale = priors_list$b_prior[2], log = TRUE) -
        dgamma(b, shape = priors_list$b_prior[1], scale = priors_list$b_prior[2], log = TRUE)
    } else {
      log_accept_ratio = log_accept_ratio - b_dash + b 
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      b <- b_dash
      log_like = logl_new
      list_accept_counts$count_accept2 = list_accept_counts$count_accept2 + 1
    } 
    #Sigma (Adpative)
    accept_prob = min(1, exp(log_accept_ratio))  
    sigma2 =  sigma2*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    
    #************************************************************************
    #c
    c_dash <- c + rnorm(1, sd = sigma3) 
    if(c_dash < 1){
      c_dash = 2 - c_dash #Prior on c: > 1
    }
    #Acceptance Probability
    logl_new = LOG_LIKE_SSI(data, a, b, c_dash)
    log_accept_ratio = logl_new - log_like 
    
    #Priors
    if(FLAGS_LIST$C_PRIOR_GAMMA){
      log_accept_ratio = log_accept_ratio + dgamma(c_dash, shape = priors_list$c_prior[1], scale = priors_list$c_prior[1], log = TRUE) -
        dgamma(c, shape = priors_list$c_prior[1], scale = priors_list$c_prior[2], log = TRUE)
    } else {
      log_accept_ratio = log_accept_ratio - priors_list$c_prior_exp[1]*c_dash + priors_list$c_prior_exp[1]*c
      if (i == 3)print('exp prior on')
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      c <- c_dash
      log_like <- logl_new
      list_accept_counts$count_accept3 = list_accept_counts$count_accept3 + 1
    }
    
    #Sigma (Adpative)
    accept_prob = min(1, exp(log_accept_ratio))  
    sigma3 =  sigma3*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    
    #*****************************************************
    #B-C TRANSFORM
    if(FLAGS_LIST$BCA_TRANSFORM){
      
      c_dash <- c + rnorm(1, sd = sigma4)
      #Prior > 1 #*** TRY WITHOUT REFLECTION!
      if(c_dash < 1){
        c_dash = 2 - c_dash 
      }
      #New b
      b_transform = ((a + b*c) - a)/c_dash #b = (r0 - a)c
      
      if(b_transform >= 0){ #Only accept values of b > 0
        
        logl_new = LOG_LIKE_SSI(data, a, b_transform, c_dash)
        log_accept_ratio = logl_new - log_like
        
        #PRIORS
        #b prior
        if (FLAGS_LIST$B_PRIOR_GAMMA) {
          tot_b_prior = dgamma(b_transform, shape = priors_list$b_prior[1], scale = priors_list$b_prior[2], log = TRUE) -
            dgamma(b, shape = priors_list$b_prior[1], scale = priors_list$b_prior[2], log = TRUE)
        } else { 
          tot_b_prior = - b_transform + b #exp(1) piror
        }
        
        #c prior
        if (FLAGS_LIST$C_PRIOR_GAMMA) {
          tot_c_prior = dgamma(c_dash, shape = priors_list$c_prior[1], scale = priors_list$c_prior[2], log = TRUE) -
            dgamma(c, shape = priors_list$c_prior[1], scale = priors_list$c_prior[2], log = TRUE)
        } else { 
          tot_c_prior = - priors_list$c_prior_exp[1]*c_dash + priors_list$c_prior_exp[1]*c 
        }
        
        #LOG ACCEPT PROB
        log_accept_ratio = log_accept_ratio + tot_b_prior + tot_c_prior 
        
        #Metropolis Step
        if (!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
          b <- b_transform
          c <- c_dash
          log_like <- logl_new
          list_accept_counts$count_accept4 = list_accept_counts$count_accept4 + 1
        }
        
        #Sigma (Adpative)
        accept_prob = min(1, exp(log_accept_ratio))  
        sigma4 = sigma4*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
      }
    }
    
    #*****************************************************
    #A-C TRANSFORM
    if(FLAGS_LIST$BCA_TRANSFORM){
      
      c_dash <- c + rnorm(1, sd = sigma5)
      #Prior > 1
      if(c_dash < 1){
        c_dash = 2 - c_dash
      }
      #New a
      a_transform = (a + b*c) - b*c_dash #a = (r0 - b*c
      
      if(a_transform >= 0){ #Only accept values of b > 0
        
        logl_new = LOG_LIKE_SSI(data, a_transform, b, c_dash)
        log_accept_ratio = logl_new - log_like
        
        #PRIORS
        #c prior
        if (FLAGS_LIST$C_PRIOR_GAMMA) {
          tot_c_prior = dgamma(c_dash, shape = priors_list$c_prior[1], scale = priors_list$c_prior[2], log = TRUE) -
            dgamma(c, shape = priors_list$c_prior[1], scale = priors_list$c_prior[2], log = TRUE)
        } else { 
          tot_c_prior = - priors_list$c_prior_exp[1]*c_dash + priors_list$c_prior_exp[1]*c 
        }
        
        #LOG ACCEPT PROB
        log_accept_ratio = log_accept_ratio + - a_transform + a + tot_c_prior 
        
        #Metropolis Step
        if (!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
          a <- a_transform
          c <- c_dash
          log_like <- logl_new
          list_accept_counts$count_accept5 = list_accept_counts$count_accept5 + 1
        }
        
        #Sigma (Adpative)
        accept_prob = min(1, exp(log_accept_ratio))  
        sigma5 = sigma5*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
        }
      }
    
    #************************************
    #DATA AUGMENTATION 
    #************************************
    if (FLAGS_LIST$DATA_AUG){

      #FOR EACH S_T
      for(t in 1:time){
        
        #Copy of data (or update as necessary)
        data_dash = data
        
        #STOCHASTIC PROPOSAL for s
        if (runif(1) < 0.5) {
          st_dash = data[[2]][t] + 1
        } else {
          st_dash = data[[2]][t] - 1 
        }
        
        #CHECK
        if (st_dash < 0) {
          st_dash = 0
        }
        
        #ACCEPTANCE PROBABILITY
        data_dash[[2]][t] = st_dash #s_t = st_dash 
        data_dash[[1]][t] =  data[[1]][t] + data[[2]][t] - st_dash #n_t = x_t - s_t
        
        if (data_dash[[1]][t] < 0){
          data_dash[[1]][t] = 0
        }
        #CRITERIA FOR S_T & N_T  
        if((data_dash[[2]][t] < 0) || (data_dash[[1]][t] < 0)){
          #print(paste0(i, t, 'WARNING'))
          #Store
          non_ss[i, t] = data[[1]][t]
          ss[i, t] = data[[2]][t]
          # print(paste0('data_dash[[1]][t = ', data_dash[[1]][t]))
          # print(paste0('data_dash[[2]][t = ', data_dash[[2]][t]))
          next  
        } 
        #print(paste0('i after next = ', i))
        
        logl_new = LOG_LIKE_SSI(data_dash, a, b, c)
        log_accept_ratio = logl_new - log_like  
        
        #METROPOLIS ACCEPTANCE STEP
        if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
          
          #ACCEPT
          data <- data_dash
          log_like <- logl_new
          mat_count_da[i, t] = mat_count_da[i, t] + 1
          list_accept_counts$count_accept6 = list_accept_counts$count_accept6 + 1
        }
        
        #Store
        non_ss[i, t] = data[[1]][t] #TAKE MEAN ACROSS MCMC DIMENSION (PLOT 0 > 50)
        ss[i, t] = data[[2]][t]
      }
    }
    
    #Loglikelihood Check (Passing - no error)
    #if (!(log_like && LOG_LIKE_SSI(data, a, b, c))) print('loglike doesnt exist')
    #else (log_like!=LOG_LIKE_SSI(data, a, b, c)) print(paste0('ERROR! logl diff = ', log_like - LOG_LIKE_SSI(data, a, b, c)))
    
    #POPPULATE MODEL PARAMETERS W/ CURRENT VALUES
    #print(paste0('i populating = ', i))
    a_vec[i] <- a; b_vec[i] <- b 
    c_vec[i] <- c; r0_vec[i] <- a + b*c
    log_like_vec[i] <- log_like #PLOT!!
    sigma$sigma1[i] = sigma1; sigma$sigma2[i] = sigma2; sigma$sigma3[i] = sigma3
    sigma$sigma4[i] = sigma4; sigma$sigma5[i] = sigma5
    
  }
  
  #Final stats
  accept_rate1 = 100*list_accept_counts$count_accept1/(n_mcmc-1)
  accept_rate2 = 100*list_accept_counts$count_accept2/(n_mcmc-1) #(list_accept_counts$count_accept2 + list_reject_counts$count_accept2)
  accept_rate3 = 100*list_accept_counts$count_accept3/(n_mcmc-1) 
  accept_rate4 = 100*list_accept_counts$count_accept4/(n_mcmc-1)
  accept_rate5 = 100*list_accept_counts$count_accept5/(n_mcmc-1)
  accept_rate6 = 100*list_accept_counts$count_accept6/((n_mcmc-1)*time) #i x t
  
  #Acceptance rates 
  list_accept_rates = list(accept_rate1 = accept_rate1,
                           accept_rate2 = accept_rate2, accept_rate3 = accept_rate3,
                           accept_rate4 = accept_rate4, accept_rate5 = accept_rate5,
                           accept_rate6 = accept_rate6)
  print(list_accept_rates)
  #Return a, acceptance rate
  return(list(a_vec = a_vec, b_vec = b_vec, c_vec = c_vec, r0_vec = r0_vec,
              log_like_vec = log_like_vec, sigma = sigma,
              list_accept_rates = list_accept_rates, 
              data = data, mat_count_da = mat_count_da, #13, 14
              non_ss = non_ss, ss = ss, #15, 16
              non_ss_mean = colMeans(non_ss),
              ss_mean = colMeans(ss))) #16 
}

