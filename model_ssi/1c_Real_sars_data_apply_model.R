#****************************************************************
#APPLY SSI MODEL MCMC + DATA AUGMENTATION
#****************************************************************

#SETUP 
library(coda)
setwd("~/GitHub/epidemic_modelling") 
source("epidemic_functions.R") 
source("plot_functions.R") 
source("helper_functions.R") 
source("model_ssi/1_SSI_model_w_data_aug.R")
source("model_ssi/2_SSI_model_Adaptive_MCMC.R")
source("model_sse/SSE_model.R")

#*********
#MCMC PARAMS  
n_mcmc = 100000 #100000 

#MCMC STARTING POINTS (RANDOM GUESS)
mod_start_points = list(m1 = 2, m2 = 0.05, m3 = 15, true_r0 = 1) #m3: 15

#SP2
mod_start_points = list(m1 = 0.5, m2 = 0.02, m3 = 22, true_r0 = 1.2) #m3: 15
#model_params_true = list(m1 = , m2 = '', m3 = '', true_r0 = '')

#SIGMA
sigma = list(sigma1 = 0.4*mod_start_points$m1, sigma2 = 0.3*mod_start_points$m2 , #Acc rate too big -> Make sigma bigger. 
             sigma3 = 0.5*mod_start_points$m3, sigma4 =  0.85*mod_start_points$m3) #Acc rate too small -> make sigma smaller

mcmc_inputs = list(n_mcmc = n_mcmc, sigma = sigma, 
                   mod_start_points = mod_start_points)

#****************************************************************
#REAL SARS DATASETS - I. CHINA
#****************************************************************

#DATA
typeX = 'Chinese'; seed_count = 1;
sim_data_china = list(china, china)
sim_data_china[[2]] = rep(1, length(sim_data_china[[1]]))

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_china = MCMC_SSI(sim_data_china, mcmc_inputs = mcmc_inputs,
                        FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                          PRIOR = TRUE, JOINT = TRUE,
                                          B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_china$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count)

df_results_da_china = PLOT_MCMC_GRID_REAL_DATA(sim_data_china, mcmc_ssi_china,
                                mcmc_plot_inputs = mcmc_plot_inputs,
                                FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                  PRIOR = TRUE,  JOINT = TRUE,
                                                  B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE,  RJMCMC = FALSE,
                                                  FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = FALSE))


#****************************************************************
# IV. CANADA
#****************************************************************

#DATA
typeX = 'Canadian'; seed_count = 1;
sim_data_canada = list(canada, canada)
sim_data_canada[[2]] = rep(1, length(sim_data_canada[[1]]))

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_canada = MCMC_SSI(sim_data_canada, mcmc_inputs = mcmc_inputs,
                           FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                             PRIOR = TRUE, JOINT = TRUE,
                                             B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_canada$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count)

df_results_da_canada = PLOT_MCMC_GRID_REAL_DATA(sim_data_canada, mcmc_ssi_canada,
                                                mcmc_plot_inputs = mcmc_plot_inputs,
                                                FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                                  PRIOR = TRUE,  JOINT = TRUE,
                                                                  B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE,  RJMCMC = FALSE,
                                                                  FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = FALSE))

  
#****************************************************************
#REAL SARS DATASETS - II. HONG KONG
#****************************************************************

#DATA
typeX = 'Hong Kong'; seed_count = 2 #seed_count + 1;
sim_data_hk = list(hk, hk)
sim_data_hk[[2]] = rep(1, length(sim_data_hk[[1]]))

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_hk = MCMC_SSI(sim_data_hk, mcmc_inputs = mcmc_inputs,
                       FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                         PRIOR = TRUE, JOINT = TRUE,
                                         B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_hk$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count)

df_results_da_hk = PLOT_MCMC_GRID_REAL_DATA(sim_data_hk, mcmc_ssi_hk,
                                  mcmc_plot_inputs = mcmc_plot_inputs,
                                  FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                    PRIOR = TRUE,  JOINT = TRUE,
                                                    B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE,  RJMCMC = FALSE,
                                                    FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = FALSE))

#****************************************************************
#REAL SARS DATASETS - III. TAIWAN
#****************************************************************

#DATA
typeX = 'Taiwan'; seed_count = seed_count + 1;
sim_data_tw = list(tw, tw)
sim_data_tw[[2]] = rep(1, length(sim_data_tw[[1]]))

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_tw = MCMC_SSI(sim_data_tw, mcmc_inputs = mcmc_inputs,
                       FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                         PRIOR = TRUE, JOINT = TRUE,
                                         B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_tw$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma,model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count)

df_results_da_tw = PLOT_MCMC_GRID_REAL_DATA(sim_data_tw, mcmc_ssi_tw,
                                  mcmc_plot_inputs = mcmc_plot_inputs,
                                  FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                    PRIOR = TRUE,  JOINT = TRUE,
                                                    B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE,  RJMCMC = FALSE,
                                                    FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = FALSE))

#****************************************************************
#REAL SARS DATASETS - V. SINGAPORE
#****************************************************************

#DATA
typeX = 'Singapore'; seed_count = 1;
sim_data_singapore = list(singapore, singapore)
sim_data_singapore[[2]] = rep(1, length(sim_data_singapore[[1]]))

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_singapore = MCMC_SSI(sim_data_singapore, mcmc_inputs = mcmc_inputs,
                              FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                PRIOR = TRUE, JOINT = TRUE,
                                                B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_singapore$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count)

df_results_da_singapore = PLOT_MCMC_GRID_REAL_DATA(sim_data_singapore, mcmc_ssi_singapore,
                                                   mcmc_plot_inputs = mcmc_plot_inputs,
                                                   FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                                     PRIOR = TRUE,  JOINT = TRUE,
                                                                     B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE,  RJMCMC = FALSE,
                                                                     FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = FALSE))

#****************************************************************
#*
#1. SSE MODEL
#*
#****************************************************************


#DATA
typeX = 'Canadian'; seed_count = seed_count = seed_count + 1;

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_sse_canada = MCMC_SSI(sim_data_canada, mcmc_inputs = mcmc_inputs,
                           FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                             PRIOR = TRUE, JOINT = TRUE,
                                             B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_canada$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count)

df_results_da_canada = PLOT_MCMC_GRID_REAL_DATA(sim_data_canada, mcmc_ssi_canada,
                                                mcmc_plot_inputs = mcmc_plot_inputs,
                                                FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                                  PRIOR = TRUE,  JOINT = TRUE,
                                                                  B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE,  RJMCMC = FALSE,
                                                                  FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = FALSE))

#****************************************************************
#*
#1.ADAPTIVE SSI: TAIWAN, CANADA
#*
#****************************************************************
n_mcmc = 100000 #5 #100000 #burn_in #5 #burn_in #100000
mcmc_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                   alpha_star = 0.4)
burn_in = 500

#****************************************************************
# IV. CANADA
#****************************************************************

#DATA
typeX = 'Canadian'; seed_count = 1;
canada_bool = canada > 1
canada_ss = as.integer(canada_bool)
sim_data_canada = list(canada, canada_ss)

#sim_data_canada = list(canada, canada)
#sim_data_canada[[2]] =  pmin(1, canada) #rep(1, length(sim_data_canada[[1]])) #

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_adp_canada3 = MCMC_SSI_ADAPTIVE(sim_data_canada, mcmc_inputs = mcmc_inputs,
                                         FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                           PRIOR = TRUE, JOINT = TRUE,
                                                           B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_adp_canada3$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count)

df_results_da_canada3 = PLOT_MCMC_GRID_REAL_DATA(sim_data_canada, mcmc_ssi_adp_canada3,
                                                 mcmc_plot_inputs = mcmc_plot_inputs,
                                                 FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                                                   PRIOR = TRUE, JOINT = TRUE,
                                                                   B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                                                   FLAG_SSI = TRUE, RJMCMC = FALSE,
                                                                   FLAG_ADAPTIVE = TRUE))

#SIGMA
par(mfrow = c(2, 3))
plot.ts(mcmc_ssi_adp_canada3$sigma$sigma1[burn_in:n_mcmc], ylab = 'sigma a', main = paste0('sigma a, burn in = ', burn_in))
plot.ts(mcmc_ssi_adp_canada3$sigma$sigma2[burn_in:n_mcmc], ylab = 'sigma b', main = paste0('sigma b, burn in = ', burn_in))
plot.ts(mcmc_ssi_adp_canada3$sigma$sigma3[burn_in:n_mcmc], ylab = 'sigma c', main = paste0('sigma c, burn in = ', burn_in))
plot.ts(mcmc_ssi_adp_canada3$sigma$sigma4[burn_in:n_mcmc], ylab = 'sigma bc', main = paste0('sigma bc, burn in = ', burn_in))
plot.ts(mcmc_ssi_adp_canada3$sigma$sigma5[burn_in:n_mcmc], ylab = 'sigma ac', main = paste0('sigma ac, burn in = ', burn_in))

#****************************************************************
#REAL SARS DATASETS - III. TAIWAN
#****************************************************************

#DATA
typeX = 'Taiwanese'; seed_count = seed_count + 1;
tw_bool = tw > 1
tw_ss = as.integer(tw_bool)
sim_data_tw = list(tw, tw_bool)

#sim_data_tw = list(tw, tw)
#sim_data_tw[[2]] = rep(1, length(sim_data_tw[[1]])) #pmin(1, tw) 

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_adp_tw3 = MCMC_SSI_ADAPTIVE(sim_data_tw, mcmc_inputs = mcmc_inputs,
                       FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                         PRIOR = TRUE, JOINT = TRUE,
                                         B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_adp_tw3$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count)

df_results_da_tw3 = PLOT_MCMC_GRID_REAL_DATA(sim_data_tw, mcmc_ssi_adp_tw3,
                                            mcmc_plot_inputs = mcmc_plot_inputs,
                                            FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                                              PRIOR = TRUE, JOINT = TRUE,
                                                              B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE,
                                                              FLAG_SSI = TRUE, RJMCMC = FALSE,
                                                              FLAG_ADAPTIVE = TRUE))

#


#Check
#LOG-LIKELIHOOD 
# nt = 0; lambda_t = 0 #0.5
# st = 0
# 
# nt = 1; lambda_t = 0 #0.5
# st = 1
# 
# logl = - lambda_t*(aX + bX) + nt*(log(aX) + log(lambda_t)) + st*(log(bX) + log(lambda_t))  + 2*log(1) - lfactorial(nt)
# logl


#****************************************************************
# PART 2: COMPARE CHAINS - GELMAN RUBIN MCMC DIAGNOSTIC :D 
#****************************************************************
#a param
#lista = mcmc.list(as.mcmc(mcmc_ssi_da2$a_vec), as.mcmc(mcmc_ssi_da3$a_vec))
gelman_a = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_adp_canadaX1$a_vec), as.mcmc(mcmc_ssi_adp_canadaX2$a_vec)))
gelman_a

#b param
gelman_b = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_adp_canadaX1$b_vec), as.mcmc(mcmc_ssi_adp_canadaX2$b_vec)))
gelman_b

#c param
gelman_c = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_adp_canadaX1$c_vec), as.mcmc(mcmc_ssi_adp_canadaX2$c_vec)))
gelman_c

#r0 vec
gelman_r0 = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_adp_canadaX1$r0_vec), as.mcmc(mcmc_ssi_adp_canadaX2$r0_vec)))
gelman_r0

#Taiwan
#a param
#lista = mcmc.list(as.mcmc(mcmc_ssi_da2$a_vec), as.mcmc(mcmc_ssi_da3$a_vec))
gelman_a = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_adp_tw2$a_vec), as.mcmc(mcmc_ssi_adp_tw3$a_vec)))
gelman_a

#b param
gelman_b = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_adp_tw2$b_vec), as.mcmc(mcmc_ssi_adp_tw3$b_vec)))
gelman_b

#c param
gelman_c = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_adp_tw2$c_vec), as.mcmc(mcmc_ssi_adp_tw3$c_vec)))
gelman_c

#r0 vec
gelman_r0 = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_adp_tw2$r0_vec), as.mcmc(mcmc_ssi_adp_tw3$r0_vec)))
gelman_r0