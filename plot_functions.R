#*******************************************
#*
#* PLOT FUNCTIONS
#* 
#******************************************
library(coda)
#**************
#PLOT MCMC MEANS
PLOT_MCMC_MEANS <- function(mcmc_output, 
                            mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                                                    model_params_true = model_params_true, mod_par_names = c('a', 'b', 'c'),
                                                    sigma = sigma, model_typeX = 'SSI', typeX = 'Individuals',
                                                    seed_count = seed_count),
                            FLAGS_LIST = list(DATA_AUG = TRUE,
                                              PRIOR = TRUE, B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                              FLAG_SSI = FALSE, RJMCMC = FALSE,
                                              BURN_IN = TRUE,
                                              FLAG_ADAPTIVE = FALSE, MULTI_ALG = FALSE)){
  #PLOT
  plot.new()
  par(mfrow = c(2,3))
  
  #EXTRACT MCMC SAMPLES
  n_mcmc = mcmc_plot_inputs$n_mcmc;
  r0_start = mcmc_plot_inputs$mod_start_points$m1 + (mcmc_plot_inputs$mod_start_points$m2*mcmc_plot_inputs$mod_start_points$m3)
  log_like_mcmc = mcmc_output$log_like_vec; log_like_mcmc = unlist(log_like_mcmc)
  
  #R0
  if (!FLAGS_LIST$MULTI_ALG ) {
    m1_mcmc = mcmc_output[1]; m1_mcmc = unlist(m1_mcmc); m2_mcmc = mcmc_output[2]; m2_mcmc = unlist(m2_mcmc); 
    m3_mcmc = mcmc_output[3]; m3_mcmc = unlist(m3_mcmc); 
    r0_mcmc = mcmc_output[4]; r0_mcmc = unlist(r0_mcmc)
  } else {
    m1_mcmc = mcmc_output$x_matrix[,1]
    m2_mcmc = mcmc_output$x_matrix[,2]
    m3_mcmc = mcmc_output$x_matrix[,3] 
    r0_mcmc = m1_mcmc + m2_mcmc*m3_mcmc
  }
  
  #BURN IN
  if (FLAGS_LIST$BURN_IN){
    burn_in = 0.05*n_mcmc
    m1_mcmc = m1_mcmc[burn_in:n_mcmc]
    m2_mcmc = m2_mcmc[burn_in:n_mcmc]
    m3_mcmc = m3_mcmc[burn_in:n_mcmc] 
    r0_mcmc = r0_mcmc[burn_in:n_mcmc]
    log_like_mcmc = log_like_mcmc[burn_in:n_mcmc]
  } else {
    burn_in = 0
  }
  
  #CUMULATIVE MEANS + PARAM SAMPLE LIMITS 
  #m1
  m1_mean = cumsum(m1_mcmc)/seq_along(m1_mcmc)
  m1_lim =  max(mcmc_plot_inputs$mod_start_points$m1[[1]], max(m1_mcmc))
  #m2
  m2_mean = cumsum(m2_mcmc)/seq_along(m2_mcmc)
  m2_lim = max(mcmc_plot_inputs$mod_start_points$m2[[1]], max(m2_mcmc))
  #m3
  m3_mean = cumsum(m3_mcmc)/seq_along(m3_mcmc)
  m3_lim =  max(mcmc_plot_inputs$mod_start_points$m3[[1]], max(m3_mcmc))
  
  #log likelihood
  ll_mean = cumsum(log_like_mcmc)/seq_along(log_like_mcmc)
  
  #PLOTS
  #m1 mean
  plot(seq_along(m1_mean), m1_mean,
       ylim=c(min(m1_mean)-0.05*min(m1_mean), max(m1_mean)+0.05*max(m1_mean)),
       xlab = 'Time', ylab =  mcmc_plot_inputs$mod_par_names[1],
       main = paste(mcmc_plot_inputs$mod_par_names[1], "MCMC mean, Start:", mcmc_plot_inputs$mod_start_points$m1),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #m2 mean
  plot(seq_along(m2_mean), m2_mean,
       ylim=c(min(m2_mean)-0.05*min(m2_mean), max(m2_mean)+0.05*max(m2_mean)),
       xlab = 'Time', ylab = 'm2',
       main = paste(mcmc_plot_inputs$mod_par_names[2], "MCMC mean, Start:", mcmc_plot_inputs$mod_start_points$m2),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       lwd = 1)
  
  #m3 Mean
  plot(seq_along(m3_mean), m3_mean,
       ylim=c(min(m3_mean)-0.05*min(m3_mean), max(m3_mean)+0.05*max(m3_mean)),
       xlab = 'Time', ylab = 'm3', 
       main = paste(mcmc_plot_inputs$mod_par_names[3], "MCMC mean, Start:", mcmc_plot_inputs$mod_start_points$m3),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       lwd = 1)
  
  #loglike mean
  plot(seq_along(ll_mean), ll_mean, #ylim=c(min(ll_mean)+0.05*min(ll_mean), max(ll_mean)+0.05*max(ll_mean)),
       xlab = 'Time', ylab = 'log likelihood', 
       main = "Log Likelihood Mean",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       lwd = 1)
  
  #R0
  r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
  plot(seq_along(r0_mean), r0_mean,
       ylim=c(min(r0_mean)-0.05*min(r0_mean), max(r0_mean)+0.05*max(r0_mean)),
       xlab = 'Time', ylab = 'R0', main = "R0 MCMC Mean",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       lwd = 1)
  
}

#**************
#PLOT ADAPTIVE SIGMA
PLOT_SIGMA_ADADPTIVE <- function(mcmc_output, mcmc_inputs){
  
  #Burn_in 
  burn_in = mcmc_inputs$burn_in_pc*mcmc_inputs$n_mcmc
  #Plot
  par(mfrow = c(2, 2))
  plot.ts(mcmc_output$sigma$sigma1[burn_in:n_mcmc], ylab = 'sigma a', main = paste0('sigma a, burn in = ', burn_in))
  plot.ts(mcmc_output$sigma$sigma2[burn_in:n_mcmc], ylab = 'sigma b', main = paste0('sigma b, burn in = ', burn_in))
  plot.ts(mcmc_output$sigma$sigma3[burn_in:n_mcmc], ylab = 'sigma c', main = paste0('sigma c, burn in = ', burn_in))
  
}

#*******************************************
#*
#* GRID PLOT SUPER-SPREADING MODELS:
#* 
#* FUNCTION TO PLOT 4x4 DASHBOARD OF MCMC RESULTS FOR SUPER SPREADING MODELS (*MCMC GRID REAL DATA SSI)
#* 
#******************************************

PLOT_MCMC_SSI_GRID_REAL_DATA <- function(sim_data, mcmc_output, 
                                         mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                                                                 mod_par_names = c('a', 'b', 'c'),
                                                                 sigma = sigma, model_typeX = 'Real', typeX = 'Events',
                                                                 seed_count = seed_count, thinning_factor = 10),
                                         priors_list = list(a_prior = c(1, 0), b_prior = c(10, 2/100), b_prior_exp = c(0.1,0), #10, 1/100
                                                            c_prior = c(10, 1), c_prior_exp = c(0.1,0)),
                                         FLAGS_LIST = list(DATA_AUG = TRUE,
                                                           PRIOR = TRUE, B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                                           FLAG_SSI = FALSE, RJMCMC = FALSE,
                                                           BURN_IN = TRUE,
                                                           FLAG_ADAPTIVE = FALSE, MULTI_ALG = FALSE, THIN = FALSE)) { 
  #PLOT
  plot.new()
  par(mfrow=c(5,5))
  
  #DATA (STARTING DATA)
  if (FLAGS_LIST$FLAG_SSI){
    non_ss_start = sim_data[[1]]; ss_start = sim_data[[2]]
    sim_data = non_ss_start + ss_start
  } 
  
  #EXTRACT MCMC SAMPLES
  r0_start = mcmc_plot_inputs$mod_start_points$m1 + (mcmc_plot_inputs$mod_start_points$m2*mcmc_plot_inputs$mod_start_points$m3)
  log_like_mcmc = mcmc_output$log_like_vec; log_like_mcmc = unlist(log_like_mcmc)
  
  if (FLAGS_LIST$MULTI_ALG){
    m1_mcmc = mcmc_output$x_matrix[,1]; m1_mcmc = unlist(m1_mcmc); m1_mcmc = m1_mcmc[!is.na(m1_mcmc)]
    m2_mcmc = mcmc_output$x_matrix[,2]; m2_mcmc = unlist(m2_mcmc); m2_mcmc = m2_mcmc[!is.na(m2_mcmc)]
    m3_mcmc = mcmc_output$x_matrix[,3]; m3_mcmc = unlist(m3_mcmc); m3_mcmc = m3_mcmc[!is.na(m3_mcmc)]
    r0_mcmc = mcmc_output$x_matrix[,1] + mcmc_output$x_matrix[,2]*mcmc_output$x_matrix[,3]
    r0_mcmc = unlist(r0_mcmc); r0_mcmc = r0_mcmc[!is.na(r0_mcmc)]
    
  } else {
    m1_mcmc = mcmc_output[1]; m1_mcmc = unlist(m1_mcmc); m1_mcmc = m1_mcmc[!is.na(m1_mcmc)]
    m2_mcmc = mcmc_output[2]; m2_mcmc = unlist(m2_mcmc);  m2_mcmc = m1_mcmc[!is.na(m1_mcmc)]
    m3_mcmc = mcmc_output[3]; m3_mcmc = unlist(m3_mcmc);  m2_mcmc = m1_mcmc[!is.na(m1_mcmc)]
    r0_mcmc = mcmc_output[4]; r0_mcmc = unlist(r0_mcmc); r0_mcmc = r0_mcmc[!is.na(r0_mcmc)] 
  }
  
  #THINNING FACTOR
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_plot_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1
    mcmc_vec_size = n_mcmc
  }
  
  #BURN IN
  if (FLAGS_LIST$BURN_IN){
    burn_in = 0.05*mcmc_vec_size
    m1_mcmc = m1_mcmc[burn_in:mcmc_vec_size]
    m2_mcmc = m2_mcmc[burn_in:mcmc_vec_size]
    m3_mcmc = m3_mcmc[burn_in:mcmc_vec_size] 
    r0_mcmc = r0_mcmc[burn_in:mcmc_vec_size]
    log_like_mcmc = log_like_mcmc[burn_in:mcmc_vec_size]
  } else {
    burn_in = 0
  }
  
  #CUMULATIVE MEANS + PARAM SAMPLE LIMITS  
  #m1
  m1_mean = cumsum(m1_mcmc)/seq_along(m1_mcmc)
  m1_lim =  max(mcmc_plot_inputs$mod_start_points$m1[[1]], max(m1_mcmc, na.rm = TRUE))
  print(paste0('a lim =', m1_lim))
  #m2
  m2_mean = cumsum(m2_mcmc)/seq_along(m2_mcmc)
  m2_lim = max(mcmc_plot_inputs$mod_start_points$m2[[1]], max(m2_mcmc, na.rm = TRUE))
  print(paste0('b lim =', m2_lim))
  #m3
  m3_mean = cumsum(m3_mcmc)/seq_along(m3_mcmc)
  m3_lim =  max(mcmc_plot_inputs$mod_start_points$m3[[1]], max(m3_mcmc, na.rm = TRUE))
  print(paste0('b lim =', m3_lim))
  #r0
  r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
  r0_lim = max(r0_start, max(r0_mcmc, na.rm = TRUE))
  print(paste0('r0 lim =', r0_lim))
  
  #PRIORS
  #m1
  m1_prior =  paste0('exp(', priors_list$a_prior[1], ')')
  #m2
  if (FLAGS_LIST$B_PRIOR_GAMMA) {
    m2_prior = paste0('Ga(', priors_list$b_prior[1], ', ', priors_list$b_prior[2], ')')
  } else {
    m2_prior = paste0('exp(', priors_list$b_prior_exp[1], ')')
  }
  #m3
  if (FLAGS_LIST$C_PRIOR_GAMMA) {
    m3_prior = paste0('1 + Ga(',   priors_list$c_prior[1], ', ',  priors_list$c_prior[2], ')')
  } else {
    m3_prior = paste0('1 + exp(',   priors_list$c_prior_exp[1], ')')
  }
  
  #******************************************************************
  #* PLOTS *
  #******************************************************************
  
  #************************
  #ROW 1: DATA INFECTIONS
  #************************
  
  #i. TOTAL INFECTIONS
  inf_tite = paste0(mcmc_plot_inputs$seed_count, ' ', mcmc_plot_inputs$model_typeX, " Data")
  plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = inf_tite, 
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #**********************
  #ii DATA - NON SS 
  #**********************
  plot.ts(non_ss_start, ylab = 'count', main = 'Non SS, Start (blk), Final (aq), Mean (or)')
  lines(mcmc_output$non_ss[mcmc_vec_size, ], col = 'aquamarine', lwd = 2) #Final
  lines(colMeans(mcmc_output$non_ss), col = '#FF3300', lwd = 2) #mean
  
  #iii. MEAN T1, T2
  plot.ts(colMeans(mcmc_output$non_ss[((mcmc_vec_size/2)+1):mcmc_vec_size, ]), ylab = 'count', 
          main = 'Non SS - Mean T1 (prp), Mean T2 (or)', col = 'orange')
  lines(colMeans(mcmc_output$non_ss[1:(mcmc_vec_size/2), ]), col = 'blueviolet', lwd = 2)
  
  #**********************
  #v, iv. DATA SS
  #**********************
  plot.ts(ss_start, ylab = 'count', main = 'SS, Start (blk), Final (bl), Mean (or)')
  lines(mcmc_output$ss[mcmc_vec_size, ], col = 'aquamarine', lwd = 2) #Final
  lines(colMeans(mcmc_output$ss), col = '#FF3300', lwd = 2) #mean
  
  #v. MEAN T1, T2
  plot.ts(colMeans(mcmc_output$ss[((mcmc_vec_size/2)+1):mcmc_vec_size, ]), ylab = 'count', 
          main = 'SS - Mean T1 (prp), Mean T2 (or)', col = 'orange')
  lines(colMeans(mcmc_output$ss[1:(mcmc_vec_size/2), ]), col = 'blueviolet', lwd = 2)
  
  #************************
  #ROW 2: MCMC TRACE PLOTS 
  #************************
  
  #****
  #a
  if (!FLAGS_LIST$FLAG_ADAPTIVE){
    plot.ts(m1_mcmc, ylab = mcmc_plot_inputs$mod_par_names[1], #ylim=c(0, m1_lim),
            main = paste(mcmc_plot_inputs$mod_par_names[1], "MCMC",
                         "Start: ", mcmc_plot_inputs$mod_start_points$m1),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  } else {
    plot.ts(m1_mcmc, ylab = paste0(mcmc_plot_inputs$mod_par_names[1], ",sigma"), #ylim=c(0, m1_lim),
            main = paste(mcmc_plot_inputs$mod_par_names[1], "MCMC",
                         "Start: ", mcmc_plot_inputs$mod_start_points$m1, ', Sigma (red)'),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    lines(mcmc_output$sigma$sigma1, col = 'red')
  }
  
  #***************
  #b
  if (!FLAGS_LIST$FLAG_ADAPTIVE){
    plot.ts(m2_mcmc, ylab = mcmc_plot_inputs$mod_par_names[3], #ylim=c(0, m2_lim), 
            main = paste(mcmc_plot_inputs$mod_par_names[2], "MCMC", 
                         "Start: ", mcmc_plot_inputs$mod_start_points$m2),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  } else {
    plot.ts(m2_mcmc, ylab = paste0(mcmc_plot_inputs$mod_par_names[2], ",sigma"), #ylim=c(0, m2_lim), 
            main = paste(mcmc_plot_inputs$mod_par_names[2], "MCMC", 
                         "Start: ", mcmc_plot_inputs$mod_start_points$m2, ', Sigma (blue)'),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    lines(mcmc_output$sigma$sigma2, col = 'blue')
  }
  
  #***************
  #c
  if (!FLAGS_LIST$FLAG_ADAPTIVE){
    plot.ts(m3_mcmc,  ylab = mcmc_plot_inputs$mod_par_names[3], #ylim=c(0, m3_lim),
            main = paste(mcmc_plot_inputs$mod_par_names[3], "MCMC",
                         "Start: ", mcmc_plot_inputs$mod_start_points$m3),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  } else {
    plot.ts(m3_mcmc,  ylab =  paste0(mcmc_plot_inputs$mod_par_names[3], ",sigma"), #ylim=c(0, m3_lim),
            main = paste(mcmc_plot_inputs$mod_par_names[3], "MCMC",
                         "Start: ", mcmc_plot_inputs$mod_start_points$m3, ', Sigma (green)'),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    lines(mcmc_output$sigma$sigma3, col = 'green')
  }
  
  #***************
  #r0
  plot.ts(r0_mcmc, ylab = "R0", #ylim=c(0, r0_lim),
          main = paste("R0 MCMC, burn-in for params =  ", burn_in),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #***************
  #LOG LIKELIHOOD
  plot.ts(log_like_mcmc, ylab = "log likelihood",
          main = paste("Log Likelihood"),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  # if (FLAGS_LIST$FLAG_ADAPTIVE){
  #   plot.ts(log_like_mcmc, ylab = "log likelihood",
  #           main = paste("Log Likelihood"),
  #           cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  # } else {
  #   plot.new()
  # }
  
  #**********************************************************
  #ROW 3:  HISTOGRAMS OF PARARMS (a, b, c, r0, loglike)
  #************************************************************
  
  #***********
  #HIST m1 
  hist(m1_mcmc, freq = FALSE, breaks = 100,
       xlab = mcmc_plot_inputs$mod_par_names[1], #ylab = 'Density', 
       main = paste(mcmc_plot_inputs$mod_par_names[1],
                    " prior:", m1_prior),
       xlim=c(0, m1_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #PRIOR PLOT
  if (FLAGS_LIST$PRIOR) {
    xseq = seq(0, 1.5, length.out = 500)
    lines(xseq, dexp(xseq, priors_list$a_prior[1]),
          type = 'l', lwd = 2, col = 'red')
  } else {
    #m2_prior = paste0('exp(', priors_list$b_prior[1], ')')
  }
  
  #***********
  #HIST m2 
  hist(m2_mcmc, freq = FALSE, breaks = 100,
       xlab = mcmc_plot_inputs$mod_par_names[2], #ylab = 'Density', 
       main = paste(mcmc_plot_inputs$mod_par_names[2],
                    " prior:", m2_prior),
       xlim=c(0, m2_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = mcmc_plot_inputs$mod_start_points$m2[[1]], col = 'blue', lwd = 2)
  
  #PRIOR PLOT
  if (FLAGS_LIST$B_PRIOR_GAMMA) {
    xseq = seq(0, 0.3, length.out = 500)
    lines(xseq, dgamma(xseq, shape =  priors_list$b_prior[1], scale =  priors_list$b_prior[2]),
          type = 'l', lwd = 2, col = 'blue')
  } else {
    xseq = seq(0, 10, length.out = 5000)
    lines(xseq, dexp(xseq, priors_list$b_prior_exp[1]),
          type = 'l', lwd = 2, col = 'blue')
  }
  
  #***********
  #Hist m3 
  hist(m3_mcmc, freq = FALSE, breaks = 100,
       xlab = mcmc_plot_inputs$mod_par_names[3], #ylab = 'Density', 
       main = paste(mcmc_plot_inputs$mod_par_names[3],
                    " prior:", m3_prior),
       xlim=c(0, m3_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = mcmc_plot_inputs$mod_start_points$m3[[1]], col = 'green', lwd = 2)#
  
  #PRIOR PLOT
  if (FLAGS_LIST$C_PRIOR_GAMMA) {
    xseq = seq(0, 35, length.out = 500)
    lines(xseq, dgamma(xseq, shape =  priors_list$c_prior[1], scale =  priors_list$c_prior[2]),
          type = 'l', lwd = 2, col = 'green')
  } else {
    xseq = seq(0, 50, length.out = 5000)
    lines(xseq, dexp(xseq, priors_list$c_prior_exp[1]),
          type = 'l', lwd = 2, col = 'green')
  }
  
  #***********
  #HIST r0 
  hist(r0_mcmc, freq = FALSE, breaks = 100,
       xlab = 'R0 total', #ylab = 'Density', 
       main = 'R0 total', #xlim=c(0, r0_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #***********
  #HIST log_like_vec 
  hist(log_like_mcmc, freq = FALSE, breaks = 100,
       xlab = 'Log likelihood', #ylab = 'Density', 
       main = 'Log likelihood',
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  # if (FLAGS_LIST$FLAG_ADAPTIVE){
  # hist(log_like_mcmc, freq = FALSE, breaks = 100,
  #      xlab = 'Log likelihood', #ylab = 'Density', 
  #      main = 'Log likelihood',
  #      cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  # } else {
  #   plot.new()
  # }
  
  #************************************************
  #ROW 4: CUMULATIVE MEAN PLOTS
  #************************************************
  
  #m1 mean
  plot(seq_along(m1_mean), m1_mean,
       ylim=c(0, m1_lim),
       xlab = 'Time', ylab =  mcmc_plot_inputs$mod_par_names[1],
       main = paste(mcmc_plot_inputs$mod_par_names[1], "MCMC mean, Start:", mcmc_plot_inputs$mod_start_points$m1),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #m2 mean
  plot(seq_along(m2_mean), m2_mean,
       ylim=c(0, m2_lim),
       xlab = 'Time', ylab = 'm2',
       main = paste(mcmc_plot_inputs$mod_par_names[2], "MCMC mean, Start:", mcmc_plot_inputs$mod_start_points$m2),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       lwd = 1)
  
  #m3 Mean
  plot(seq_along(m3_mean), m3_mean,
       xlab = 'Time', ylab = 'm3', 
       main = paste(mcmc_plot_inputs$mod_par_names[3], "MCMC mean, Start:", mcmc_plot_inputs$mod_start_points$m3),
       ylim=c(0, m3_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       lwd = 1)
  
  #r0 mean
  plot(seq_along(r0_mean), r0_mean, #ylim=c(0, r0_lim),
       xlab = 'Time', ylab = 'R0', main = "R0 MCMC Mean",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       lwd = 1)
  
  #loglike mean
  plot(seq_along(cumsum(log_like_mcmc)/seq_along(log_like_mcmc)), cumsum(log_like_mcmc)/seq_along(log_like_mcmc),
       xlab = 'Time', ylab = 'log likelihood', 
       main = "Log Likelihood Mean",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       lwd = 1)
  # if (FLAGS_LIST$FLAG_ADAPTIVE){
  # plot(seq_along(cumsum(log_like_mcmc)/seq_along(log_like_mcmc)), cumsum(log_like_mcmc)/seq_along(log_like_mcmc),
  #      xlab = 'Time', ylab = 'log likelihood', 
  #      main = "Log Likelihood Mean",
  #      cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
  #      lwd = 1)
  # } else {
  #   plot.new()
  # }
  #*****************
  #ROW 5: JOINT DISTRIBUTIONS/MARGINALS
  #********************
  
  #a vs r0
  plot(m1_mcmc, r0_mcmc,
       xlab = mcmc_plot_inputs$mod_par_names[1], ylab = 'R0',
       main = paste0(mcmc_plot_inputs$mod_par_names[1], ' vs R0'),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       cex = 0.5)
  
  #b*c vs r0
  plot(m2_mcmc*m3_mcmc, r0_mcmc,
       xlab = 'b*c', ylab = 'R0',
       main = 'b*c vs R0',
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       cex = 0.5)
  
  #m1 vs m2
  plot(m1_mcmc, m2_mcmc,
       xlab = mcmc_plot_inputs$mod_par_names[1], ylab = mcmc_plot_inputs$mod_par_names[2],
       main = paste(mcmc_plot_inputs$mod_par_names[1], 'vs', mcmc_plot_inputs$mod_par_names[2]),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       cex = 0.5)
  
  #m1 vs m3
  plot(m1_mcmc, m3_mcmc,
       xlab = mcmc_plot_inputs$mod_par_names[1], ylab = mcmc_plot_inputs$mod_par_names[3], main = paste(mcmc_plot_inputs$mod_par_names[1], 'vs', mcmc_plot_inputs$mod_par_names[3]),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       cex = 0.5)
  
  #v. m2 vs m3
  plot(m2_mcmc, m3_mcmc,
       xlab = mcmc_plot_inputs$mod_par_names[2], ylab = mcmc_plot_inputs$mod_par_names[3], main = paste(mcmc_plot_inputs$mod_par_names[2], 'vs', mcmc_plot_inputs$mod_par_names[3]),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       cex = 0.5)
  
  #********************
  #v. DATAFRAME
  #********************
  
  #FINAL MEAN Stats
  data_10_pc = 0.5*mcmc_vec_size #50%
  m1_mean_tail = round(mean(m1_mcmc[mcmc_vec_size - data_10_pc:mcmc_vec_size]), 2) 
  m2_mean_tail = round(mean(m2_mcmc[mcmc_vec_size - data_10_pc:mcmc_vec_size]), 2)
  m3_mean_tail = round(mean(m3_mcmc[mcmc_vec_size - data_10_pc:mcmc_vec_size]), 2)
  m4_mean_tail = round(mean(r0_mcmc[mcmc_vec_size - data_10_pc:mcmc_vec_size]), 2)
  
  if (!FLAGS_LIST$MULTI_ALG){
    df_results <- data.frame(
      rep = mcmc_plot_inputs$seed_count,
      mcmc_vec_size = mcmc_vec_size,
      m1 = mcmc_plot_inputs$mod_start_points$m1[[1]],
      m1_mc = m1_mean_tail,
      m2 = mcmc_plot_inputs$mod_start_points$m2[[1]],
      m2_mc = m2_mean_tail,
      m3 = mcmc_plot_inputs$mod_start_points$m3[[1]],
      m3_mc = m3_mean_tail,
      R0 = r0_start, 
      R0_mc = m4_mean_tail,
      accept_rate_m1 = round(mcmc_output$list_accept_rates$accept_rate1, 2),
      a_rte_m2 = round(mcmc_output$list_accept_rates$accept_rate2, 2),
      a_rte_m3 = round(mcmc_output$list_accept_rates$accept_rate3, 2),
      a_rte_m2_m3 = round(mcmc_output$list_accept_rates$accept_rate4, 2),
      a_rte_m1_m3 = round(mcmc_output$list_accept_rates$accept_rate5, 2),
      a_rte_d_aug = round(mcmc_output$list_accept_rates$accept_rate6, 2),
      a_es = effectiveSize(as.mcmc(m1_mcmc))[[1]],
      b_es = effectiveSize(as.mcmc(m2_mcmc))[[1]],
      c_es = effectiveSize(as.mcmc(m3_mcmc))[[1]],
      d_es = effectiveSize(as.mcmc(r0_mcmc))[[1]],
      time_elap = format(mcmc_output$time_elap, format = "%H:%M:%S")[1])
  } else {
    
    df_results <- data.frame(
      rep = mcmc_plot_inputs$seed_count,
      mcmc_vec_size = mcmc_vec_size,
      m1 = mcmc_plot_inputs$mod_start_points$m1[[1]],
      m1_mc = m1_mean_tail,
      m2 = mcmc_plot_inputs$mod_start_points$m2[[1]],
      m2_mc = m2_mean_tail,
      m3 = mcmc_plot_inputs$mod_start_points$m3[[1]],
      m3_mc = m3_mean_tail,
      R0 = r0_start, 
      R0_mc = m4_mean_tail,
      accept_rate_x = round(mcmc_output$accept_rate, 2),
      a_rte_d_aug = round(mcmc_output$accept_rate_da, 2),
      a_es = effectiveSize(as.mcmc(m1_mcmc))[[1]],
      b_es = effectiveSize(as.mcmc(m2_mcmc))[[1]],
      c_es = effectiveSize(as.mcmc(m3_mcmc))[[1]],
      d_es = effectiveSize(as.mcmc(r0_mcmc))[[1]],
      time_elap = format(mcmc_output$time_elap, format = "%H:%M:%S")[1])
    
  }
  
  print(df_results)
  return(df_results)
  
}

#TEST
x1 <- c(1, 4, 3, NA, 7)
