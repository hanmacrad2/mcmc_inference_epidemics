#SIMULATE FROM v ~ Gamma(R0, k), Z ~ NegBin(R0, k)

#PACKAGES
library(simstudy)

#*********************************************************
#I. BASELINE SIMULATION
SIMULATE_BRANCHING_NEGBIN_SSI = function(num_days = 110, alphaX = 1.2, shape_gamma = 6, scale_gamma = 1, k = 0.16) {
  
  'Simulate from the Negative Binomial model'
  
  #INTIALISE VECTORS
  x = vector('numeric', num_days); x[1] = 2
  
  #INFECTIOUSNESS (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  vec_nu = vector('numeric', num_days); #sum of nus 
  #DAYS OF THE EPIDEMIC
  for (t in 2:num_days) {
    
    #Total rate
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)])) #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    r0 = alphaX*lambda_t
    
    #GET GAMMA PARAMS, SHAPE & SCALE FROM THE MEAN
    gamma_params <- gammaGetShapeRate(r0, k)
    c(gamma_params$shape, gamma_params$rate)
    
    #INDIVIDUAL R0, I.E NU
    vec_nu[t] <- sum(rgamma(x[tt], shape = gamma_params$shape, rate = gamma_params$rate))
    nu_t = sum(rgamma(x[tt], shape = gamma_params$shape, rate = gamma_params$rate))
    x[t] = rpois(1, mu = nu_t*lambda_t, size = k) #Poisson 
    
  }
  return(x)
}

x2 = SIMULATE_BRANCHING_NEGBIN_SSI()
plot.ts(x2)
    
    #tot_infectiousness = tot_infectiousness + sum(vec_nu)
    
#     #FOR EACH PREVIOUS DATY 
#     tot_infectiousness = 0
#    
#     for (tt in 1:t){
#       
#       #PARAMS
#       lambda_t = sum(x[1:(tt-1)]*rev(prob_infect[1:(tt-1)])) #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
#       r0 = alphaX*lambda_t
#       
#       #GET GAMMA PARAMS, SHAPE & SCALE FROM THE MEAN
#       gamma_params <- gammaGetShapeRate(r0, k)
#       c(gamma_params$shape, gamma_params$rate)
#       
#       #INDIVIDUAL R0, I.E NU
#       vec_nu <- rgamma(x[tt], shape = gamma_params$shape, rate = gamma_params$rate)
#       tot_infectiousness = tot_infectiousness + sum(vec_nu)
#       
#     }
#     
#     x[t] = rnbinom(1, mu = tot_infectiousness, size = k) #Poisson 
#     
#     #IDEA
#     #tot_rate = (vec_nu1 + vecnu_2 + vecnu_3 + vecnu4)*rev(prob_infect[1:(t-1)])
#     }
#   
#   return(x)
# }

#Apply
y = SIMULATE_BRANCHING_NEGBIN()
plot.ts(y)


#ATTEMPT TWO
SIMULATE_BRANCHING_NEGBIN_SSE = function(num_days = 110, alphaX = 1.2, shape_gamma = 6, scale_gamma = 1, k = 0.16) {
  
  'Simulate from the Negative Binomial model'
  
  #INTIALISE VECTORS
  x = vector('numeric', num_days); x[1] = 2
  
  #INFECTIOUSNESS (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #DAYS OF THE EPIDEMIC
  for (t in 2:num_days) {
    print(paste0('t = ', t))
    #Total rate
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)])) #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    r0 = alphaX*lambda_t
    print(r0)
    
    #GET GAMMA PARAMS, SHAPE & SCALE FROM THE MEAN
    #gamma_params <- gammaGetShapeRate(r0, k)
    #c(gamma_params$shape, gamma_params$rate)
    
    x[t] = rnbinom(1, mu = r0, size = k) 
    
  }
  
  return(x)
}

#Apply
y = SIMULATE_BRANCHING_NEGBIN()
plot.ts(y)
    
    #FOR EACH PREVIOUS DATY 
    #tot_infectiousness = 0
    #vec_nu = rgamma(sum(x))
    
   # for (tt in 1:t){ #x = [2, 0, 3] rgamma(sum(x))
      
      #PARAMS
      #lambda_t = sum(x[1:(tt-1)]*rev(prob_infect[1:(tt-1)])) #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
      #r0 = alphaX*lambda_t
      
      #GET GAMMA PARAMS, SHAPE & SCALE FROM THE MEAN
      #gamma_params <- gammaGetShapeRate(r0, k)
      #c(gamma_params$shape, gamma_params$rate)
      
      #INDIVIDUAL R0, I.E NU
      vec_nu <- rgamma(x[tt], shape = gamma_params$shape, rate = gamma_params$rate)
      
      tot_infectiousness = tot_infectiousness + sum(vec_nu)
      
#    }
    
    x[t] = rnbinom(1, mu = r0, size = k) 
    
    #IDEA
    #tot_rate = (vec_nu1 + vecnu_2 + vecnu_3 + vecnu4)*rev(prob_infect[1:(t-1)])
#   }
#   
#   return(x)
# }

#Apply
y2 = SIMULATE_BRANCHING_NEGBIN()
plot.ts(y2)

data_nb = y

#NOTE
#SSI
#If you want to ensure each individual has a unique R0 then;
# v ~ Gamma(R0, k) #If x1 = 3, v1 + v2 + v3 (For every day)
# Z ~ Poisson(sum(v)) 