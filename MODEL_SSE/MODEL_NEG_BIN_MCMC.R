#MODEL NEGATIVE BINOMIAL II

#****************************
#SIMULATION
SIMULATE_NEGBIN_SSI = function(num_days = 110, alphaX = 1.2,
                               shape_gamma = 6, scale_gamma = 1, k = 0.16) {
  
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
# #DATA AUGMENTATION
