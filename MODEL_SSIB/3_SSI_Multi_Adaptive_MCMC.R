#****************************************************************
#1. SSI MODEL MCMC + MULTI-PARAMETER STEP
#****************************************************************
#library(MASS)

LOG_LIKE_SSI_MULTI <- function(data, x){

  #PARAMS
  aX =  x[1]; bX = x[2]; cX = x[3]
  #print(paste0('aX = ', aX))
  #print(paste0('bX = ', bX))
  #print(paste0('cX = ', cX))

  #data
  n = data[[1]]; s = data[[2]]

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
    logl_add = - lambda_t*(aX + bX) + n[t]*(log(aX) + log(lambda_t)) + s[t]*(log(bX) + log(lambda_t))  + 2*log(1) - lfactorial(n[t]) - lfactorial(s[t])
    logl = logl + logl_add

    # if (is.na(logl_add)) {
    #   print(paste0('WARNING NaN,  logl = ', logl))
    #   # print(paste0('lambda_t = ', lambda_t))
    #   # print(paste0('n[t] = ', n[t]))
    #   # print(paste0('s[t] = ', s[t]))
    #   next
    # } else if (is.infinite(logl_add)) {
    #   print(paste0('WARNING Inf,  log_like_add = ', log_like_add))
    #   # print(paste0('lambda_t = ', lambda_t))
    #   # print(paste0('n[t] = ', n[t]))
    #   # print(paste0('s[t] = ', s[t]))
    #   next
    # } else {
    #   logl = logl + logl_add
    # }
  }

  logl
  #print(logl)
}

#************************************************************************
#1. SSI MCMC ADAPTIVE SHAPING                            (W/ data AUGMENTATION OPTION)
#************************************************************************

MCMC_SSI_MULTI_ADADPT <- function(data,
                                  mcmc_inputs = list(n_mcmc = n_mcmc,
                                                     mod_start_points = mod_start_points,  burn_in_pc = 0.05,
                                                     dim = 3, alpha_star = 0.4, v0 = 100, vec_min = c(0,0,0),
                                                     thinning_factor = 10),
                                  priors_list = list(b_prior_exp = c(1,0), c_prior_exp = c(0.1,0)),
                                  FLAGS_LIST = list(DATA_AUG = TRUE, B_PRIOR_EXP = TRUE, C_PRIOR_EXP = TRUE,
                                                    THIN = TRUE)){

  #**********************************************
  #INITIALISE PARAMS
  #**********************************************

  #REPLACE: n = i-1 (Simon's paper)

  #MCMC PARAMS + VECTORS
  time = length(data[[1]]); n_mcmc = mcmc_inputs$n_mcmc;
  dim = mcmc_inputs$dim; count_accept = 0

  #THINNING FACTOR
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1
    mcmc_vec_size = n_mcmc
  }

  #SIZE NEEDS TO
  x_matrix = matrix(NA, mcmc_vec_size, dim);   #Changed from 0 to NA (As should be overwriting all cases )
  x_matrix[1,] <- mcmc_inputs$mod_start_points;
  x = x_matrix[1,] #3x1 #as.matrix
  log_like_vec <- vector('numeric', mcmc_vec_size);
  log_like_vec[1] <- LOG_LIKE_SSI_MULTI(data, x);  log_like = log_like_vec[1]
  lambda_vec <- vector('numeric', mcmc_vec_size); lambda_vec[1] <- 1

  #ADAPTIVE SHAPING PARAMS + VECTORS
  c_star = (2.38^2)/dim; termX = mcmc_inputs$v0 + dim
  delta = 1/(mcmc_inputs$alpha_star*(1 - mcmc_inputs$alpha_star))
  sigma_i = diag(3); lambda_i = 1

  #DATA AUG OUTPUT
  non_ss = matrix(0, mcmc_vec_size, time) #USE THINNING FACTOR
  ss = matrix(0, mcmc_vec_size, time) #USE THINNING FACTOR
  count_accept_da = 0

  for(i in 2:n_mcmc) {

    #PRINT PROGRESS
    if(i%%(n_mcmc/50) == 0) print(paste0('i = ', i))

    #SIGMA ITERATION NO.1
    if (i == 2){
      x_bar = 0.5*(x + x_matrix[1,])
      sigma_i = (1/(termX + 3))*(tcrossprod(x_matrix[1,]) + tcrossprod(x) -
                                   2*tcrossprod(x_bar) + (termX + 1)*sigma_i) #CHANGE TO USE FUNCTIONS
    }

    #PROPOSAL
    x_dash = c(x + mvrnorm(1, mu = rep(0, dim), Sigma = lambda_i*c_star*sigma_i)) #Vectorise using c()

    #ONLY KEEP POSTIVE
    if (min(x_dash - mcmc_inputs$vec_min) < 0){ #REJECT IF C - 1 <0 I.E C < 1 #Model specific. #

      #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
      if (i%%thinning_factor == 0) {
        x_matrix[i/thinning_factor,] = x  # i%/%
        log_like_vec[i/thinning_factor] <- log_like
        lambda_vec[i/thinning_factor] <- lambda_i #Taking role of sigma, overall scaling constant. Sigma becomes estimate of the covariance matrix of the posterior
      }

      next
    }

    #LOG LIKE
    logl_new = LOG_LIKE_SSI_MULTI(data, x_dash)

    #ACCEPTANCE RATIO
    log_accept_ratio = logl_new - log_like #- x_dash + x

    #PRIORS
    log_accept_ratio = log_accept_ratio - x_dash[1] + x[1] #PRIOR ON A
    if (FLAGS_LIST$B_PRIOR_EXP){
      log_accept_ratio = log_accept_ratio - x_dash[2] + x[2] #PRIOR ON B
    }
    if (FLAGS_LIST$C_PRIOR_EXP){
      log_accept_ratio = log_accept_ratio - priors_list$c_prior_exp[1]*x_dash[3] + priors_list$c_prior_exp[1]*x[3] #PRIOR ON C
    }

    #METROPOLIS ACCEPTANCE STEP
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      x <- x_dash
      count_accept = count_accept + 1
      log_like = logl_new
    }

    #SIGMA - ADAPTIVE SHAPING
    if ( i > 2) {
      #SIGMA - ADAPTIVE SHAPING
      xbar_prev = x_bar
      x_bar = (i-1)/i*xbar_prev + (1/i)*x
      sigma_i = (1/(i + termX + 1))*( (i + termX)*sigma_i +tcrossprod(x)
                                      + (i-1)*tcrossprod(xbar_prev)
                                      -i*tcrossprod(x_bar))
    }

    #LAMBDA - ADAPTIVE SCALING
    accept_prob = min(1, exp(log_accept_ratio))
    lambda_i =  lambda_i*exp(delta/i*(accept_prob - mcmc_inputs$alpha_star))

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
          #Store
          non_ss[i/thinning_factor, t] = data[[1]][t]
          ss[i/thinning_factor, t] = data[[2]][t]
          next
        }

        #LOG LIKELIHOOD FOR METROPOLIS STEP
        logl_new = LOG_LIKE_SSI_MULTI(data_dash, x)
        log_accept_ratio = logl_new - log_like

        #METROPOLIS ACCEPTANCE STEP
        if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {

          #ACCEPT
          data <- data_dash
          log_like <- logl_new
          #mat_count_da[i, t] = mat_count_da[i, t] + 1
          accept_count_da = count_accept_da + 1
        }

        #STORE NEW DATA POINT
        non_ss[i/thinning_factor, t] = data[[1]][t] #TAKE MEAN ACROSS MCMC DIMENSION (PLOT 0 > 50)
        ss[i/thinning_factor, t] = data[[2]][t]
      }
    }

    #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0) {
      x_matrix[i/thinning_factor,] = x  # i%/%
      log_like_vec[i/thinning_factor] <- log_like
      lambda_vec[i/thinning_factor] <- lambda_i #Taking role of sigma, overall scaling constant. Sigma becomes estimate of the covariance matrix of the posterior
    }

  } #END FOR LOOP

  #Final stats
  accept_rate = 100*count_accept/(n_mcmc-1)
  accept_rate_da = 100*count_accept_da/((n_mcmc-1)*time)

  #print(list_accept_rates)
  #Return a, acceptance rate
  return(list(x_matrix = x_matrix, log_like_vec = log_like_vec, lambda_vec = lambda_vec,
              accept_rate = accept_rate, accept_rate_da = accept_rate_da,
              non_ss = non_ss, ss = ss, #15, 16
              non_ss_mean = colMeans(non_ss),
              ss_mean = colMeans(ss))) #16
}

#Apply
#mcmc_multiX1 = MCMC_SSI_MULTI_ADADPT(sim_data_canadaX1, mcmc_inputs)

#NOTE
#NO REFLECTION, NO TRANSFORMS, MORE INTELLIGENT ADAPTATION
