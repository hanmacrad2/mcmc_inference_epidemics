#SIMULATE FROM v ~ Gamma(R0, k), Z ~ NegBin(R0, k)

#PACKAGES
library(simstudy)

#*********************************************************
#I. BASELINE SIMULATION
SIMULATE_BRANCHING_NEGBIN = function(num_days = 110, alphaX = 1.2, shape_gamma = 6, scale_gamma = 1, k = 0.16) {
  
  'Negative binomial simulation model'
  
  #Initialise params
  x = vector('numeric', num_days)
  x[1] = 2
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Total rate
    #lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)])) #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    #r0 = alphaX*lambda_t
    
    #For Each Previous day 
    tot_infectiousness = 0
    for (tt in 1:t){
      
      #Params
      lambda_t = sum(x[1:(tt-1)]*rev(prob_infect[1:(tt-1)])) #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
      r0 = alphaX*lambda_t
      #Get gamma params, shape and scale from mean 
      gamma_params <- gammaGetShapeRate(r0, k)
      c(gamma_params$shape, gamma_params$rate)
      
      #INDIVIDUAL R0, I.E NU
      vec_nu <- rgamma(x[tt], shape = gamma_params$shape, rate = gamma_params$rate)
      tot_infectiousness = tot_infectiousness + sum(vec_nu)
      
    }
    
    x[t] = rnbinom(1, mu = tot_infectiousness, size = k) 
    
    #IDEA
    #tot_rate = (vec_nu1 + vecnu_2 + vecnu_3 + vecnu4)*rev(prob_infect[1:(t-1)])
    }
  
  return(x)
}

#Apply
y = SIMULATE_BRANCHING_NEGBIN()
plot.ts(y)