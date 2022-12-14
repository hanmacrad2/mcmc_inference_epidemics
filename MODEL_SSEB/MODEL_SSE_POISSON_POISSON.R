#SSE MODEL - POISSON-POISSON COMPOUND

#PARAMETERS
alphaX = 0.8; betaX = 0.1; gammaX = 10

#1. LOG LIKELIHOOD
LOG_LIKE_SSE_POISSON <- function(x, alphaX, betaX, gammaX){
  
  #Params
  num_days = length(x)
  shape_gamma = 6; scale_gamma = 1
  
  #Infectiousness (Discrete gamma) i,e Prob less than x2 - prob less than x1; the area in between 
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  logl = 0
  
  for (t in 2:num_days) {
    
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)]))
    inner_sum_xt = 0
    
    for (yt in 0:x[t]){ #Sum for all values of yt
      
      #Log likelihood
      zt = x[t] - yt
      inner_sum_xt = inner_sum_xt + 
                        exp(-alphaX*lambda_t)*(1/factorial(yt))*(alphaX*lambda_t)^yt*
                        PROBABILITY_ZT(zt, lambda_t, alphaX, betaX, gammaX)
                      

      
    } 
    
    logl = logl + log(inner_sum_xt) 
    
  }
  
  return(logl)
  
}

loglike1 = LOG_LIKE_SSE_POISSON_COMPOUND(x, alphaX, betaX, gammaX)
loglike1

#2. PROBABILITY OF ZT
PROBABILITY_ZT <- function(zt, lambda_t, alphaX, betaX, gammaX, max_nt = 5) {
  
  'Probability of zt'
  
  #Initialise
  prob_zt = 0

  for (nt in 0:max_nt){
    
    prob_zt = prob_zt + dpois(nt, betaX*lambda_t)*dpois(zt, gammaX*nt)
    
  }
  
  return(prob_zt)
}

#Apply
prob_zt2 = PROB_ZT(x, alphaX, betaX, gammaX)
prob_zt2

kk = vector("numeric");  count = 1
for (t in days){
  kk[count] = t
  count = count + 1
}
kk

#************
#* PROB ZT
PLOT_PROB_ZT <- function(x, alphaX, betaX, gammaX, max_nt = 5, days = c(15,19,36,25,96,102)) {
  
  'Plot total probability of zt for certain days of epidemic'
  
  #INITIALISE
  num_days = length(x); shape_gamma = 6; scale_gamma = 1
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  #Plot
  plot.new(); par(mfrow = c(2,3))
  
  #SPECIFIC DAY
  for (t in days) {
    print(paste0('t: ', t))
    #t
    vec_prob_zt = vector("numeric"); vec_prob_yzt = vector("numeric")
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)]))
    total_prob = 0
    
    for (zt in 0:x[t]){
      prob_zt = 0; prob_yt = dpois((x[t] - zt), alphaX*lambda_t)
      for (nt in 0:max_nt){
        
        prob_zt = prob_zt + dpois(nt, betaX*lambda_t)*dpois(zt, gammaX*nt)
       
      }
      
    vec_prob_zt[zt+1] = prob_zt
    prob_yt = dpois((x[t] - zt), alphaX*lambda_t)
    total_prob = total_prob + prob_zt*prob_yt
    vec_prob_yzt[zt+1] = prob_zt*prob_yt #renormalize at end of loop
    }
    vec_prob_yzt =  vec_prob_yzt/total_prob
    
  #Plot
    title = paste0('P(zt|X_1:t, alpha, beta, gamma), Day ', t, '. Xt = ', x[t], ' yt = Xt - zt')
  plot(seq_along(vec_prob_yzt) - 1, vec_prob_yzt, col = 'blue',
       xlab = 'z_t', ylab = 'P(z_t)', #type = '*',
       main = title, lwd = 3)
  #print(paste0('Should be = 1. total prob = ', total_prob))
  print(vec_prob_yzt)
  print(sum(vec_prob_yzt))
  }
}

#PLOT
PLOT_PROB_ZT(canadaX, alphaX, betaX, gammaX)

#***********
#*PLOT PROBABILITY OF NT

#TOTAL PROB
PLOT_PROB_NT <- function(x, alphaX, betaX, gammaX, max_nt = 5) {
  
  'Probability of zt'
  plot.new()
  par(mfrow = c(2,3))
  
  for (zt in c(0,1,3,5,10,15)){
  
  #Initialise (max = 25, 97)
  vec_prob_zt = vector("numeric"); vec_nt =  vector("numeric");  vec_nt2 =  vector("numeric")
  num_days = length(x); prob_zt = 0; shape_gamma = 6; scale_gamma = 1
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  for (t in 2:num_days){
    #print(paste0('t = ', t))
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)]))
    
    for (nt in 0:max_nt){
      
      #print(paste0('lambda_t: ', round(lambda_t, 3),  ' betaX*lambda_t: ',  round(betaX*lambda_t, 3))) 
      prob_t = dpois(nt, betaX*lambda_t)*dpois(zt, gammaX*nt)
      vec_nt[nt] = prob_t 
      prob_zt = prob_zt + prob_t
      #print(paste0('nt: ', nt, ', prob_t: ', round(prob_t, 4), ', rate nt: ', round(betaX*lambda_t, 4),  ' rate zt: ', gammaX*nt))
    }
    
    if (t == 25){
      vec_nt2 = vec_nt
    }
    
    #Title
    if (zt == 0){
      title = paste0('P(Zt = ', zt, '), max SSE events K = 5. Day 97 (blue), 25(rd), 99(gr)')
    } else {
      title = paste0('P(Zt = ', zt, '), max SSE events K = 5')
    } 
    
    #Plot for nt
    if (t == 97){
      print(paste0('t = ', t))
      plot(seq_along(vec_nt), vec_nt, col = 'blue',
           xlab = 'n_t', ylab = 'P(n_t)', #type = '*',
              main = title, lwd = 3)
      points(seq_along(vec_nt2), vec_nt2, pch="*", col = 'red', lwd = 3)
    } else if (t == 99){
      points(seq_along(vec_nt), vec_nt, pch="*", col = 'green', lwd = 3)
    }
    #print(paste0('prob_t = ', prob_t))
    vec_prob_zt[t] = prob_t
  }
  
  #plot.ts(vec_prob_zt)
  
  }
  #return(prob_zt)
}

#Apply
PLOT_PROB_NT(canadaX, alphaX, betaX, gammaX)


#BRAINSTORM
SSE_DPOIS <- function(xt, alphaX, betaX, gammaX, max_ss = 10) {
  
  ss_lim = min(xt, max_ss)
  ss_poi_prob = 0
  
  for (z_t in 0:ss_lim){
    print(z_t)
    ss_poi_prob = ss_poi_prob + dpois(z_t, gammaX*(x_t - z_t), log = TRUE) #dpois or manual
    
  }
  
  return(ss_poi_prob) #Probability? :D See if it works :)
}

#Log Likelihood 
log_like_sse <- function(x, alphaX, betaX, gammaX){
  
  #Params
  num_days = length(x)
  shape_gamma = 6
  scale_gamma = 1
  
  #Infectiousness (Discrete gamma) i,e Prob less than x2 - prob less than x1; the area in between 
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  logl = 0
  
  for (t in 2:num_days) {
    
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)]))
    inner_sum_xt = 0
    
    for (y_t in 0:x[t]){ #Sum for all values of y_t
      
      #Log likelihood
      inner_sum_xt = (inner_sum_xt + exp(-alphaX*lambda_t)*(1/factorial(y_t))*(alphaX*lambda_t)^y_t*
                        (gamma((x[t] - y_t) + betaX*lambda_t))/(gamma(betaX*lambda_t)*
                                                                  factorial(x[t] - y_t))*(1/(gammaX +1))^(betaX*lambda_t)*
                        (gammaX/(gammaX + 1))^(x[t] - y_t))
      
    } 
    
    logl = logl + log(inner_sum_xt) 
    
  }
  
  logl
  
}


#**********************************
#MCMC
#***********************************

#************************************************************************
#1. SSE MCMC
#************************************************************************
SSE_POI_MCMC_ADAPTIVE <- function(epidemic_data,
                              mcmc_inputs = list(n_mcmc = 500000,
                                                 mod_start_points = list(m1 = 0.8, m2 = 0.1, m3 = 10), alpha_star = 0.4,
                                                 thinning_factor = 10), #10
                              priors_list = list(alpha_prior_exp = c(1, 0), beta_prior_ga = c(10, 2/100), beta_prior_exp = c(0.1,0),
                                                 gamma_prior_ga = c(10, 1), gamma_prior_exp = c(0.1,0)),
                              FLAGS_LIST = list(ADAPTIVE = TRUE, ABG_TRANSFORM = TRUE,
                                                PRIOR = TRUE, BETA_PRIOR_GA = FALSE, GAMMA_PRIOR_GA = FALSE,
                                                THIN = TRUE)) {
  
  'Returns MCMC samples of SSI model parameters (alpha, beta, gamma, r0 = alpha + beta*gamma)
  w/ acceptance rates.
  INCLUDES; ADAPTATION, beta-gamma & alpha-gamma transform'
  
  'Priors
  p(a) = exp(rate) = rate*exp(-rate*x). log(r*exp(-r*x)) = log(r) - rx
      -> E.g exp(1) = 1*exp(-1*a) = exp(-a). log(exp(-a)) = - a
  p(b) = exp(1) or p(b) = g(shape, scale), for e.g g(3, 2)
  p(c) = exp(1) + 1 = 1 + exp(-c) = exp(c - 1)'
  
  
  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  n_mcmc = mcmc_inputs$n_mcmc;
  print(paste0('num mcmc iters = ', n_mcmc))
  
  #THINNING FACTOR
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1; mcmc_vec_size = n_mcmc
  }
  
  #INITIALISE MCMC VECTORS
  alpha_vec <- vector('numeric', mcmc_vec_size); beta_vec <- vector('numeric', mcmc_vec_size)
  gamma_vec <- vector('numeric', mcmc_vec_size); r0_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec <- vector('numeric', mcmc_vec_size);
  
  #INITIALISE MCMC[1]
  alpha_vec[1] <- mcmc_inputs$mod_start_points$m1; beta_vec[1] <- mcmc_inputs$mod_start_points$m2
  gamma_vec[1] <- mcmc_inputs$mod_start_points$m3; r0_vec[1] <- alpha_vec[1] + beta_vec[1]*gamma_vec[1]
  log_like_vec[1] <- LOG_LIKE_SSE_POISSON(epidemic_data, alpha_vec[1], beta_vec[1], gamma_vec[1])
  
  #INITIALISE RUNNING PARAMS
  alpha = alpha_vec[1]; beta = beta_vec[1]; gamma = gamma_vec[1]; log_like = log_like_vec[1]
  
  #SIGMA
  sigma1 =  0.4*mcmc_inputs$mod_start_points$m1;  sigma2 = 0.3*mcmc_inputs$mod_start_points$m2
  sigma3 = 0.4*mcmc_inputs$mod_start_points$m3; sigma4 = 0.85*mcmc_inputs$mod_start_points$m3
  sigma5 = 0.85*mcmc_inputs$mod_start_points$m3
  
  #SIGMA; INITIALISE FOR ADAPTIVE MCMC
  if (FLAGS_LIST$ADAPTIVE){
    
    #SIGMA
    sigma1_vec <- vector('numeric', mcmc_vec_size); sigma2_vec <- vector('numeric', mcmc_vec_size)
    sigma3_vec <- vector('numeric', mcmc_vec_size); sigma4_vec <- vector('numeric', mcmc_vec_size)
    sigma5_vec <- vector('numeric', mcmc_vec_size);
    
    #SIGMA; INITIALISE FIRST ELEMENT
    sigma1_vec[1] =  sigma1; sigma2_vec[1] =  sigma2; sigma3_vec[1] =  sigma3
    sigma4_vec[1] =  sigma4; sigma5_vec[1] =  sigma5
    
    #SIGMA; List of sigma vectors for each iteration of the MCMC algorithm
    sigma = list(sigma1_vec = sigma1_vec, sigma2_vec = sigma2_vec, sigma3_vec = sigma3_vec,
                 sigma4_vec = sigma4_vec, sigma5_vec = sigma5_vec)
    
    #Other adaptive parameters
    delta = 1/(mcmc_inputs$alpha_star*(1-mcmc_inputs$alpha_star))
    
  } else {
    
    #SIGMA; List of sigma vectors for each iteration of the MCMC algorithm
    sigma = list(sigma1 = sigma1, sigma2 = sigma2,
                 sigma3 = sigma3, sigma4 = sigma4,
                 sigma5 = sigma5)
  }
  
  #INITIALISE: ACCEPTANCE COUNTS
  list_accept_counts = list(count_accept1 = 0, count_accept2 = 0, count_accept3 = 0,
                            count_accept4 = 0, count_accept5 = 0)
  
  #******************************
  #MCMC CHAIN
  #******************************
  for(i in 2:n_mcmc) {
    
    if (i%%100 == 0) {
      
      print(paste0('i = ', i))
    }
    
    #****************************************************** s
    #alpha
    alpha_dash <- alpha + rnorm(1, sd = sigma1)
    
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }
    
    #log a
    logl_new = LOG_LIKE_SSE_POISSON(epidemic_data, alpha_dash, beta, gamma)
    log_accept_ratio = logl_new - log_like  #+ prior1 - prior
    #Priors
    if (FLAGS_LIST$PRIOR){
      log_accept_ratio = log_accept_ratio - alpha_dash + alpha #*Actually this is the Acceptance RATIO. ACCEPTANCE PROB = MIN(1, EXP(ACCPET_PROB))
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      alpha <- alpha_dash
      list_accept_counts$count_accept1 = list_accept_counts$count_accept1 + 1
      log_like = logl_new
    }
    
    #Sigma (Adaptive)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma1 =  sigma1*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }
    
    #************************************************************************ Only if (b > 0)
    #beta
    beta_dash <- beta + rnorm(1, sd = sigma2)
    if(beta_dash < 0){
      beta_dash = abs(beta_dash)
    }
    
    #loglikelihood
    logl_new = LOG_LIKE_SSE_POISSON(epidemic_data, alpha, beta_dash, gamma)
    log_accept_ratio = logl_new - log_like
    
    #Priors
    if (FLAGS_LIST$BETA_PRIOR_GA){
      log_accept_ratio = log_accept_ratio +
        dgamma(beta_dash, shape = priors_list$beta_prior_ga[1], scale = priors_list$beta_prior_ga[2], log = TRUE) -
        dgamma(beta, shape = priors_list$beta_prior_ga[1], scale = priors_list$beta_prior_ga[2], log = TRUE)
    } else {
      log_accept_ratio = log_accept_ratio - beta_dash + beta
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      beta <- beta_dash
      log_like = logl_new
      list_accept_counts$count_accept2 = list_accept_counts$count_accept2 + 1
    }
    
    #Sigma (Adpative)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma2 =  sigma2*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }
    
    #************************************************************************
    #gamma
    gamma_dash <- gamma + rnorm(1, sd = sigma3)
    if(gamma_dash < 1){
      gamma_dash = 2 - gamma_dash #Prior on c: > 1
    }
    #Acceptance Probability
    logl_new = LOG_LIKE_SSE_POISSON(epidemic_data, alpha, beta, gamma_dash)
    log_accept_ratio = logl_new - log_like
    
    #Priors
    if(FLAGS_LIST$GAMMA_PRIOR_GA){
      log_accept_ratio = log_accept_ratio + dgamma(gamma_dash, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[1], log = TRUE) -
        dgamma(gamma, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[2], log = TRUE)
    } else {
      log_accept_ratio = log_accept_ratio - priors_list$gamma_prior_exp[1]*gamma_dash + priors_list$gamma_prior_exp[1]*gamma
      if (i == 3) print('exp prior on')
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      gamma <- gamma_dash
      log_like <- logl_new
      list_accept_counts$count_accept3 = list_accept_counts$count_accept3 + 1
    }
    
    #Sigma (Adpative)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma3 =  sigma3*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }
    
    #*****************************************************
    #Beta-Gamma TRANSFORM
    if(FLAGS_LIST$ABG_TRANSFORM){
      gamma_dash <- gamma + rnorm(1, sd = sigma4)
      
      #Prior > 1 #* TRY WITHOUT REFLECTION
      if(gamma_dash < 1){
        gamma_dash = 2 - gamma_dash
      }
      #New b
      beta_transform = ((alpha + beta*gamma) - alpha)/gamma_dash #beta = (r0 - a)c
      
      if( beta_transform >= 0){ #Only accept values of beta> 0
        
        logl_new = LOG_LIKE_SSE_POISSON(epidemic_data, alpha, beta_transform, gamma_dash)
        log_accept_ratio = logl_new - log_like
        
        #PRIORS
        #Beta prior
        if (FLAGS_LIST$BETA_PRIOR_GA) {
          tot_beta_prior = dgamma(beta_transform, shape = priors_list$beta_prior_ga[1], scale = priors_list$beta_prior_ga[2], log = TRUE) -
            dgamma(beta, shape = priors_list$beta_prior_ga[1], scale = priors_list$beta_prior_ga[2], log = TRUE)
        } else {
          tot_beta_prior = - beta_transform + beta #exp(1) prior
        }
        
        #gamma prior
        if (FLAGS_LIST$GAMMA_PRIOR_GA) {
          tot_gamma_prior = dgamma(gamma_dash, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[2], log = TRUE) -
            dgamma(gamma, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[2], log = TRUE)
        } else {
          tot_gamma_prior = - priors_list$gamma_prior_exp[1]*gamma_dash + priors_list$gamma_prior_exp[1]*gamma
        }
        
        #LOG ACCEPT PROB
        log_accept_ratio = log_accept_ratio + tot_beta_prior + tot_gamma_prior
        
        #Metropolis Step
        if (!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
          beta <- beta_transform
          gamma <- gamma_dash
          log_like <- logl_new
          list_accept_counts$count_accept4 = list_accept_counts$count_accept4 + 1
        }
        
        #Sigma (Adpative)
        if (FLAGS_LIST$ADAPTIVE){
          accept_prob = min(1, exp(log_accept_ratio))
          sigma4 = sigma4*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
        }
      }
    }
    
    #*****************************************************
    #Alpha-Gamma TRANSFORM
    if(FLAGS_LIST$ABG_TRANSFORM){
      
      gamma_dash <- gamma+ rnorm(1, sd = sigma5)
      #Prior > 1
      if(gamma_dash < 1){
        gamma_dash = 2 - gamma_dash
      }
      #New alpha
      alpha_transform = (alpha + beta*gamma) - beta*gamma_dash #alpha = (r0 - beta*gamma)
      
      if( alpha_transform >= 0){ #Only accept values of beta> 0
        
        logl_new = LOG_LIKE_SSE_POISSON(epidemic_data, alpha_transform, beta, gamma_dash)
        log_accept_ratio = logl_new - log_like
        
        #PRIORS
        #gamma prior
        if (FLAGS_LIST$GAMMA_PRIOR_GA) {
          tot_gamma_prior = dgamma(gamma_dash, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[2], log = TRUE) -
            dgamma(gamma, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[2], log = TRUE)
        } else {
          tot_gamma_prior = - priors_list$gamma_prior_exp[1]*gamma_dash + priors_list$gamma_prior_exp[1]*gamma
        }
        
        #LOG ACCEPT PROB
        log_accept_ratio = log_accept_ratio - alpha_transform + alpha + tot_gamma_prior
        
        #Metropolis Step
        if (!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
          alpha <- alpha_transform
          gamma <- gamma_dash
          log_like <- logl_new
          list_accept_counts$count_accept5 = list_accept_counts$count_accept5 + 1
        }
        
        #Sigma (Adpative)
        if (FLAGS_LIST$ADAPTIVE){
          accept_prob = min(1, exp(log_accept_ratio))
          sigma5 = sigma5*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
        }
      }
    }
    
    #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0) {
      #print(paste0('i = ', i))
      i_thin = i/thinning_factor
      alpha_vec[i_thin] <- alpha; beta_vec[i_thin] <- beta
      gamma_vec[i_thin] <- gamma; r0_vec[i_thin] <- alpha + beta*gamma
      log_like_vec[i_thin] <- log_like
      sigma$sigma1_vec[i_thin] = sigma1; sigma$sigma2_vec[i_thin] = sigma2; sigma$sigma3_vec[i_thin] = sigma3
      sigma$sigma4_vec[i_thin] = sigma4; sigma$sigma5_vec[i_thin] = sigma5
    }
  }
  
  #Final stats
  accept_rate1 = 100*list_accept_counts$count_accept1/(n_mcmc-1)
  accept_rate2 = 100*list_accept_counts$count_accept2/(n_mcmc-1) #(list_accept_counts$count_accept2 + list_reject_counts$count_accept2)
  accept_rate3 = 100*list_accept_counts$count_accept3/(n_mcmc-1)
  accept_rate4 = 100*list_accept_counts$count_accept4/(n_mcmc-1)
  accept_rate5 = 100*list_accept_counts$count_accept5/(n_mcmc-1)
  
  #Acceptance rates
  list_accept_rates = list(accept_rate1 = accept_rate1,
                           accept_rate2 = accept_rate2, accept_rate3 = accept_rate3,
                           accept_rate4 = accept_rate4, accept_rate5 = accept_rate5)
  print(list_accept_rates)
  
  #Return a, acceptance rate
  return(list(alpha_vec = alpha_vec, beta_vec = beta_vec, gamma_vec = gamma_vec, r0_vec = r0_vec,
              log_like_vec = log_like_vec, sigma = sigma,
              list_accept_rates = list_accept_rates))
}

#Apply

#MCMC SPECS
mcmc_specs = list(model_type = 'SSI', n_mcmc = 500000, 
                  mod_start_points = list(m1 = 0.8, m2 = 0.1, m3 = 10),
                  mod_par_names = c('alpha', 'beta', 'gamma'),
                  seed_count = 1, burn_in_pc = 0.05, thinning_factor = 10)

#2. RUN SSE MODEL
mcmc_sse_output = SSE_POI_MCMC_ADAPTIVE(canadaX)
saveRDS(mcmc_sse_output, file = 'mcmc_sse_output_poisson_comp_.rds')


#***********************
#* PLOT SAMPLING DISTRIBUTION OF ZT
#Params
alphaX = 0.8; betaX = 0.1; gammaX = 10

lambda_t =
  
apply_zt_range <- function(max_zt = 15){
  
  #Intialise
  vec_zt = vector("numeric")
  PROBABILITY_ZT <- function(zt, lambda_t, alphaX, betaX, gammaX, max_nt = 5)
    
  for (i in 1:max_zt) {
    vec_zt[i] = PROBABILITY_ZT(i, )
  }
  
}
