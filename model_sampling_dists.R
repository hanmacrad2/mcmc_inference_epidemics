#MODEL SAMPLING DISTRIBUTIONS

PLOT_SAMP_DIST_SSE <- function(seedX = 1, FLAG_SSE = TRUE,
                              FREQ_SETTING = FALSE,
                              num_days = 110,
                              num_samps = 50000,
                              range_alpha = c(1.0, 1.2, 1.4, 1.6)){
  
  #Plot
  set.seed(seedX); print(paste0('seed = ', seedX))
  plot.new(); par(mfrow = c(2, 2))
  count = 0;
  
  #Outer loop: dispersion
  for (alphaXX in range_alpha) {
    print(alphaXX)
    vec_infections = vector("numeric", length = num_samps)
    
    for (i in 1:num_samps){
      
      #Simulate for a given alpha, k
      if (FLAG_SSE){
        y = simulate_branching_ss(alphaX = alphaXX)
      } else {
        y = SIMULATE_NU_SSE(alphaX = alphaXX, k = kX)  
      }
      
      #ADD
      #print(paste0('sum y = ', sum(y)))
      vec_infections[i] = sum(y)
    }
    #PLOT
    titleX = bquote("Sum of xt, T = 110"  ~ alpha ~ .(alphaXX) ~ beta ~ "0.2" ~ gamma ~ "10")
    #titleX = bquote("Sum of xt, T = " ~ .(num_days) ~ ", xt ~ Compound Poisson()" ~ alpha ~ .(alphaXX) ~ beta ~ "0.2" ~ gamma ~ "10")

    hist(vec_infections, main = titleX, freq = FREQ_SETTING, xlab = 'sum of infections', #ylab = 'Daily infection count', 
         cex.lab=1.5, cex.axis=1.5, cex.main= 2.0, cex.sub=1.5)
    #count = count + 1
  }
  
}

#PLOT
PLOT_SAMP_DIST_SSE()

#plot
titleX = bquote("Sum of xt, T = 110"  ~ alpha ~ .(alphaXX) ~ beta ~ "0.2" ~ gamma ~ "10")
