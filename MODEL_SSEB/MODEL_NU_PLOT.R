#SIMULATE FROM v ~ Gamma(R0, k), Z ~ NegBin(R0, k)

#PACKAGES
library(simstudy)


#NOTE
#SSI
#If you want to ensure each individual has a unique R0 then;
# v ~ Gamma(R0, k) #If x1 = 3, v1 + v2 + v3 (For every day)
# Z ~ Poisson(sum(v)) 

#*********************************************************
#I. SIMULATE SSE NEGBIN
#*********************************************************

SIMULATE_NU_SSE = function(num_days = 110, shape_gamma = 6, scale_gamma = 1,
                                         alphaX = 1.2, k = 0.16) {
  
  'Simulate from the Negative Binomial model'
  
  #INTIALISE VECTORS
  x = vector('numeric', num_days); x[1] = 2
  
  #INFECTIOUSNESS (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #DAYS OF THE EPIDEMIC
  for (t in 2:num_days) {
    
    #Total rate
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)])) #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    r0 = alphaX*lambda_t
    
    x[t] = rnbinom(1, mu = r0, size = k) 
    
  }
  
  return(x)
}

#Apply
y = SIMULATE_NEGBIN_SSE()
plot.ts(y)

#*********************************************************
#II. SIMULATE SSI NEGBIN
#*********************************************************

SIMULATE_NU_SSI = function(num_days = 110, alphaX = 1.2,
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

#APPLY
y3 = SIMULATE_NU_SSI()
plot.ts(y3)

#*********************************************************
#SAMPLING DISTRIBUTION
#*********************************************************

PLOT_SAMP_DIST_NU <- function(seedX = 1, FLAG_SSI = TRUE,
                                   FREQ_SETTING = FALSE,
                                   num_days = 110,
                                   num_samps = 100000, alphaXX = 1.35,
                               range_k = c(0.1, 0.5, 1, 4)){
  
  #Plot
  set.seed(seedX); print(paste0('seed = ', seedX))
  plot.new(); par(mfrow = c(2, 2))
  count = 0;
  
  #Outer loop: dispersion
  for (kX in range_k) {
  vec_infections = vector("numeric", length = num_samps)
  
  for (i in 1:num_samps){
      
      #Simulate for a given alpha, k
      if (FLAG_SSI){
        y = SIMULATE_NU_SSI(alphaX = alphaXX, k = kX)  
      } else {
        y = SIMULATE_NU_SSE(alphaX = alphaXX, k = kX)  
      }
      
    #ADD
    #print(paste0('sum y = ', sum(y)))
    vec_infections[i] = sum(y)
  }
    #PLOT
   print(paste0('kX  = ', sum(kX)))
    if (count == 0){
      titleX = bquote("Sum of xt, T = " ~ .(num_days) ~ "xt ~ Poisson(" ~ nu ~ "*" ~ lambda[t] ~ ")," ~ alpha ~ .(alphaXX) ~       # ~ nu ~ "~ Ga(" ~ .(alphaXX)
                      ~ ", k = " ~ .(kX))
    } else {
      titleX = bquote(alpha ~ "=" ~ .(alphaXX)
                      ~ ", k = " ~ .(kX))
    }
    #Plot
    hist(vec_infections, main = titleX, freq = FREQ_SETTING, xlab = 'sum of infections', #ylab = 'Daily infection count', 
            cex.lab=1.5, cex.axis=1.5, cex.main= 2.0, cex.sub=1.5)
    #count = count + 1
  }
  
}

#APPLY
PLOT_SAMP_DIST_NEG_BIN()


#*********************************************************
#APPLY TO A RANGE OF PARAMS
#*********************************************************

PLOT_RANGE_NU_VALS <- function(seedX = 1, FLAG_SSI = TRUE,
                               range_rate = c(1.35, 1.35, 1.35, 1.35),
                               range_k = c(0.1, 0.5, 1, 4)){
  
  #Plot
  set.seed(seedX); print(paste0('seed = ', seedX))
  plot.new(); par(mfrow = c(4, 4))
  count = 0
  
  #Outer loop: dispersion
  for (kX in range_k){
    
    for (alphaXX in range_rate) {
      
      #Simulate for a given alpha, k
      if (FLAG_SSE){
        y = SIMULATE_NU_SSI(alphaX = alphaXX, k = kX)  
      } else {
        y = SIMULATE_NU_SSE(alphaX = alphaXX, k = kX)  
      }
      
      if (count == 0){
        titleX = bquote("Negbin(" ~ alpha ~ "*" ~ lambda[t] ~ ", k)," ~ alpha ~ "=" ~ .(alphaXX)
                        ~ ", k = " ~ .(kX))
      } else {
        titleX = bquote(alpha ~ "=" ~ .(alphaXX)
                        ~ ", k = " ~ .(kX))
      }
      #Plot
      plot.ts(y, main = titleX, ylab = 'Daily infection count', 
              cex.lab=1.5, cex.axis=1.5, cex.main= 2.0, cex.sub=1.5)
      count = count + 1
      
    }
  }
}

#APPLY
PLOT_RANGE_NU_VALS()


#OLD!
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
    #gamma_params = c(shape = k, scale = alphaX/k)
    vec_nu_t[t-1] <- rgamma(1, shape = x[t-1]*k, scale = alphaX/k)
    #gamma_params <- gammaGetShapeRate(alphaX, k) 
    #c(gamma_params$shape, gamma_params$rate)
    #vec_nu_t[t-1] <- sum(rgamma(x[t-1], shape = gamma_params$shape, rate = gamma_params$rate))
    #vec_nu_t[t-1] <- rgamma(1, shape = x[t-1]*gamma_params$shape, rate = gamma_params$rate)
    
    #OFFSPRING; POISSON()
    infectivity = rev(prob_infect[1:(t-1)]) #x[1:(t-1)]*rev(prob_infect[1:(t-1)])
    total_rate = sum(vec_nu_t*infectivity)
    print(total_rate)
    x[t] = rpois(1, total_rate)
    
  }
  return(x)
}

#NEGATIVE BINOMIAL
#PLOT
v1= rnbinom(10000, mu=2.4,size=0.16);
v2= rnbinom(10000,mu=1.2,size=0.16) + rnbinom(10000,mu=1.2,size=0.16);

c(var(v1),var(v2))


#PARAMS CORRECT?

