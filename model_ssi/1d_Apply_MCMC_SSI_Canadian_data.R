#****************************************************************
#* 
#IV. #DATA CANADIAN SOURCE II (Xavier's)
#*
#****************************************************************
library(coda)
setwd("~/GitHub/epidemic_modelling") 
source("epidemic_functions.R") 
source("plot_functions.R") 
source("helper_functions.R") 
source("model_ssi/1_SSI_model_w_data_aug.R")
source("model_ssi/2_SSI_model_Adaptive_MCMC.R")

#****************************************************************
# DATA 
#****************************************************************
typeX = 'Canadian XD'; seed_count = 1;

canadaX = c(1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 2, 3, 0, 2, 2,
            5, 7, 9, 7, 3, 4, 1, 4, 5, 7, 7, 7, 7, 3, 3, 5, 3, 5, 7, 4, 4, 2,
            3, 1, 1, 1, 0, 0, 2, 1, 3, 2, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1,
            1, 1, 0, 1, 0, 0, 0, 1, 0, 2, 0, 3, 2, 2, 1, 2, 3, 4, 5, 5, 4, 6,
            6, 4, 8, 5, 6, 7, 5, 9, 1, 2, 3, 1, 1, 2, 0, 0, 0, 2, 0, 0, 0, 1)

#************************
#DATA I Extreme - No SS
canada_ss = rep(0, length(canadaX))
sim_data_canadaX1 = list(canadaX, canada_ss)

#************************
#DATA II Extreme SS; [1 0]
canada_bool = canadaX > 1
canada_ss = as.integer(canada_bool)
canada_ns = canadaX - canada_ss
sim_data_canadaX2 = list(canada_ns, canada_ss)

#****************************************************************
# MCMC 
#****************************************************************
n_mcmc = 100000 #5 
mod_start_points = list(m1 = 0.5, m2 = 0.02, m3 = 22, true_r0 = 1.2) #list(m1 = 2, m2 = 0.05, m3 = 15, true_r0 = 1) #m3: 15
mcmc_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                   alpha_star = 0.4)
burn_in = 500

#****************************************************************
# DATA I (Starting point: No SS)
#****************************************************************

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_adp_canadaX1 = MCMC_SSI_ADAPTIVE(sim_data_canadaX1, mcmc_inputs = mcmc_inputs,
                                          FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                            PRIOR = TRUE, JOINT = TRUE,
                                                            B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_adp_canadaX1$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count)

df_results_da_canadaX1 = PLOT_MCMC_SSI_GRID_REAL_DATA(sim_data_canadaX1, mcmc_ssi_adp_canadaX1,
                                                      mcmc_plot_inputs = mcmc_plot_inputs,
                                                      FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                                                        PRIOR = TRUE, JOINT = TRUE,
                                                                        B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                                                        FLAG_SSI = TRUE, RJMCMC = FALSE,
                                                                        FLAG_ADAPTIVE = TRUE))

#SIGMA
par(mfrow = c(2, 3))
plot.ts(mcmc_ssi_adp_canadaX1$sigma$sigma1[burn_in:n_mcmc], ylab = 'sigma a', main = paste0('sigma a, trim val = ', burn_in))
plot.ts(mcmc_ssi_adp_canadaX1$sigma$sigma2[burn_in:n_mcmc], ylab = 'sigma b', main = paste0('sigma b, trim val = ', burn_in))
plot.ts(mcmc_ssi_adp_canadaX1$sigma$sigma3[burn_in:n_mcmc], ylab = 'sigma c', main = paste0('sigma c, trim val = ', burn_in))
plot.ts(mcmc_ssi_adp_canadaX1$sigma$sigma4[burn_in:n_mcmc], ylab = 'sigma bc', main = paste0('sigma bc, trim val = ', burn_in))
plot.ts(mcmc_ssi_adp_canadaX1$sigma$sigma5[burn_in:n_mcmc], ylab = 'sigma ac', main = paste0('sigma ac, trim val = ', burn_in))

#****************************************************************
# IV. CANADA - XAVIERS (Starting point: No SS) + DATA AUG HALF WAY
#****************************************************************
#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_adp_canadaX1B = MCMC_SSI_ADAPTIVE(sim_data_canadaX1, mcmc_inputs = mcmc_inputs,
                                           FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                             PRIOR = TRUE, JOINT = TRUE,
                                                             B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_adp_canadaX1B$time_elap = time_elap

df_results_da_canadaX1B = PLOT_MCMC_SSI_GRID_REAL_DATA(sim_data_canadaX1, mcmc_ssi_adp_canadaX1B,
                                                       mcmc_plot_inputs = mcmc_plot_inputs,
                                                       FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                                                         PRIOR = TRUE, JOINT = TRUE,
                                                                         B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                                                         FLAG_SSI = TRUE, RJMCMC = FALSE,
                                                                         FLAG_ADAPTIVE = TRUE))

#SIGMA
par(mfrow = c(2, 3))
plot.ts(mcmc_ssi_adp_canadaX1B$sigma$sigma1[burn_in:n_mcmc], ylab = 'sigma a', main = paste0('sigma a, trim val = ', burn_in))
plot.ts(mcmc_ssi_adp_canadaX1B$sigma$sigma2[burn_in:n_mcmc], ylab = 'sigma b', main = paste0('sigma b, trim val = ', burn_in))
plot.ts(mcmc_ssi_adp_canadaX1B$sigma$sigma3[burn_in:n_mcmc], ylab = 'sigma c', main = paste0('sigma c, trim val = ', burn_in))
plot.ts(mcmc_ssi_adp_canadaX1B$sigma$sigma4[burn_in:n_mcmc], ylab = 'sigma bc', main = paste0('sigma bc, trim val = ', burn_in))
plot.ts(mcmc_ssi_adp_canadaX1B$sigma$sigma5[burn_in:n_mcmc], ylab = 'sigma ac', main = paste0('sigma ac, trim val = ', burn_in))


#****************************************************************
# IV. CANADA - XAVIERS (Starting point: No SS) + NO ADAPTIVE
#****************************************************************

sigma = list(sigma1 = 0.4*mod_start_points$m1, sigma2 = 0.3*mod_start_points$m2 , #Acc rate too big -> Make sigma bigger. 
             sigma3 = 0.5*mod_start_points$m3, sigma4 =  0.85*mod_start_points$m3) #Acc rate too small -> make sigma smaller

mcmc_inputs = list(n_mcmc = n_mcmc, sigma = sigma, 
                   mod_start_points = mod_start_points)
#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_adp_canadaX1C = MCMC_SSI(sim_data_canadaX1, mcmc_inputs = mcmc_inputs,
                                  FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                    PRIOR = TRUE, JOINT = TRUE,
                                                    B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_adp_canadaX1C$time_elap = time_elap

df_results_da_canadaX1C = PLOT_MCMC_SSI_GRID_REAL_DATA(sim_data_canadaX1, mcmc_ssi_adp_canadaX1C,
                                                       mcmc_plot_inputs = mcmc_plot_inputs,
                                                       FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                                                         PRIOR = TRUE, JOINT = TRUE,
                                                                         B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                                                         FLAG_SSI = TRUE, RJMCMC = FALSE,
                                                                         FLAG_ADAPTIVE = FALSE))

#****************************************************************
# IV. CANADA - XAVIERS (Starting point: SS; [1, 0]) + NO ADAPTIVE
#****************************************************************

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_adp_canadaX2C = MCMC_SSI(sim_data_canadaX2, mcmc_inputs = mcmc_inputs,
                                  FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                    PRIOR = TRUE, JOINT = TRUE,
                                                    B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_adp_canadaX2C$time_elap = time_elap

df_results_da_canadaX2C = PLOT_MCMC_SSI_GRID_REAL_DATA(sim_data_canadaX2, mcmc_ssi_adp_canadaX2C,
                                                       mcmc_plot_inputs = mcmc_plot_inputs,
                                                       FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                                                         PRIOR = TRUE, JOINT = TRUE,
                                                                         B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                                                         FLAG_SSI = TRUE, RJMCMC = FALSE,
                                                                         FLAG_ADAPTIVE = FALSE))



#****************************************************************
# IV. CANADA - XAVIERS (Starting point: Many SS)
#****************************************************************

#DATA
typeX = 'Canadian XD'; seed_count = seed_count + 1;

#Extreme 2: SS many
canada_bool = canadaX > 1
canada_ss = as.integer(canada_bool)
canada_ns = canadaX - canada_ss
sim_data_canadaX2 = list(canada_ns, canada_ss)


#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_adp_canadaX2 = MCMC_SSI_ADAPTIVE(sim_data_canadaX2, mcmc_inputs = mcmc_inputs,
                                          FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                            PRIOR = TRUE, JOINT = TRUE,
                                                            B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_adp_canadaX2$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count)

df_results_da_canadaX2 = PLOT_MCMC_SSI_GRID_REAL_DATA(sim_data_canadaX2, mcmc_ssi_adp_canadaX2,
                                                      mcmc_plot_inputs = mcmc_plot_inputs,
                                                      FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                                                        PRIOR = TRUE, JOINT = TRUE,
                                                                        B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                                                        FLAG_SSI = TRUE, RJMCMC = FALSE,
                                                                        FLAG_ADAPTIVE = TRUE))



#****************************************************************
# IV. CANADA - XAVIERS (Starting point: Many SS)
#****************************************************************

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_adp_canadaX2 = MCMC_SSI_ADAPTIVE(sim_data_canadaX2, mcmc_inputs = mcmc_inputs,
                                          FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                            PRIOR = TRUE, JOINT = TRUE,
                                                            B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_adp_canadaX2$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count)

df_results_da_canadaX2 = PLOT_MCMC_SSI_GRID_REAL_DATA(sim_data_canadaX2, mcmc_ssi_adp_canadaX2,
                                                      mcmc_plot_inputs = mcmc_plot_inputs,
                                                      FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                                                        PRIOR = TRUE, JOINT = TRUE,
                                                                        B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                                                        FLAG_SSI = TRUE, RJMCMC = FALSE,
                                                                        FLAG_ADAPTIVE = TRUE))

#SIGMA
par(mfrow = c(2, 3))
plot.ts(mcmc_ssi_adp_canadaX2$sigma$sigma1[burn_in:n_mcmc], ylab = 'sigma a', main = paste0('sigma a, trim val = ', burn_in))
plot.ts(mcmc_ssi_adp_canadaX2$sigma$sigma2[burn_in:n_mcmc], ylab = 'sigma b', main = paste0('sigma b, trim val = ', burn_in))
plot.ts(mcmc_ssi_adp_canadaX2$sigma$sigma3[burn_in:n_mcmc], ylab = 'sigma c', main = paste0('sigma c, trim val = ', burn_in))
plot.ts(mcmc_ssi_adp_canadaX2$sigma$sigma4[burn_in:n_mcmc], ylab = 'sigma bc', main = paste0('sigma bc, trim val = ', burn_in))
plot.ts(mcmc_ssi_adp_canadaX2$sigma$sigma5[burn_in:n_mcmc], ylab = 'sigma ac', main = paste0('sigma ac, trim val = ', burn_in))

#****************************************************************
#*
#1.# ADAPTIVE AND NO DATA AUG, CANADA X
#*
#****************************************************************
n_mcmc = 100000 #5 #100000 #burn_in #5 #burn_in #100000
mcmc_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                   alpha_star = 0.4)
#burn_in = 500

#****************************************************************
# N0 1 START NO SS
#****************************************************************

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_adp_canadaX4 = MCMC_SSI_ADAPTIVE(sim_data_canadaX1, mcmc_inputs = mcmc_inputs,
                                          FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                            PRIOR = TRUE, JOINT = TRUE,
                                                            B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_adp_canadaX4$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count)

df_results_da_canadaX4 = PLOT_MCMC_SSI_GRID_REAL_DATA(sim_data_canadaX1, mcmc_ssi_adp_canadaX4,
                                                      mcmc_plot_inputs = mcmc_plot_inputs,
                                                      FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                                                        PRIOR = TRUE, JOINT = TRUE,
                                                                        B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                                                        FLAG_SSI = TRUE, RJMCMC = FALSE,
                                                                        FLAG_ADAPTIVE = TRUE))

#NO 2
#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_adp_canadaX5 = MCMC_SSI_ADAPTIVE(sim_data_canadaX2, mcmc_inputs = mcmc_inputs,
                                          FLAGS_LIST = list(DATA_AUG = FALSE, BCA_TRANSFORM = TRUE,
                                                            PRIOR = TRUE,
                                                            B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_adp_canadaX5$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count)

df_results_da_canadaX5 = PLOT_MCMC_SSI_GRID_REAL_DATA(sim_data_canadaX2, mcmc_ssi_adp_canadaX5,
                                                      mcmc_plot_inputs = mcmc_plot_inputs,
                                                      FLAGS_LIST = list(DATA_AUG = FALSE,
                                                                        PRIOR = TRUE, B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                                                        FLAG_SSI = TRUE, RJMCMC = FALSE,
                                                                        BURN_IN = TRUE, FLAG_ADAPTIVE = TRUE))

#****************************************************************
#*
#1.# ADAPTIVE AND NO TRANSFORMS
#*
#****************************************************************
n_mcmc = 500000 
mcmc_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                   alpha_star = 0.4, burn_in_pc = 0.05)

#****************************************************************
# N0 1 START NO SS
#****************************************************************

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_adp_canadaX8 = MCMC_SSI_ADAPTIVE(sim_data_canadaX1, mcmc_inputs = mcmc_inputs,
                                          FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = FALSE,
                                                            PRIOR = TRUE, JOINT = TRUE,
                                                            B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_adp_canadaX8$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count)

df_results_da_canadaX8 = PLOT_MCMC_SSI_GRID_REAL_DATA(sim_data_canadaX1, mcmc_ssi_adp_canadaX8,
                                                      mcmc_plot_inputs = mcmc_plot_inputs,
                                                      FLAGS_LIST = list(DATA_AUG = FALSE,
                                                                        PRIOR = TRUE, B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                                                        FLAG_SSI = TRUE, RJMCMC = FALSE,
                                                                        BURN_IN = TRUE, FLAG_ADAPTIVE = TRUE))
#Plots
PLOT_MCMC_MEANS(mcmc_ssi_adp_canadaX8,
                mcmc_plot_inputs = mcmc_plot_inputs)

PLOT_SIGMA_ADADPTIVE(mcmc_ssi_adp_canadaX8, mcmc_inputs)

#****************************************************************
# No. II START SS; 1, 0
#****************************************************************
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_adp_canadaX9 = MCMC_SSI_ADAPTIVE(sim_data_canadaX2, mcmc_inputs = mcmc_inputs,
                                          FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = FALSE,
                                                            PRIOR = TRUE, JOINT = TRUE,
                                                            B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_adp_canadaX7$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count)

df_results_da_canadaX9 = PLOT_MCMC_SSI_GRID_REAL_DATA(sim_data_canadaX2, mcmc_ssi_adp_canadaX9,
                                                      mcmc_plot_inputs = mcmc_plot_inputs,
                                                      FLAGS_LIST = list(DATA_AUG = FALSE,
                                                                        PRIOR = TRUE, B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                                                        FLAG_SSI = TRUE, RJMCMC = FALSE,
                                                                        BURN_IN = TRUE, FLAG_ADAPTIVE = TRUE))
#Plots
PLOT_MCMC_MEANS(mcmc_ssi_adp_canadaX9,
                mcmc_plot_inputs = mcmc_plot_inputs)
#Plot
PLOT_SIGMA_ADADPTIVE(mcmc_ssi_adp_canadaX9, mcmc_inputs)

#Plots
PLOT_MCMC_MEANS(mcmc_ssi_adp_canadaX8,
                mcmc_plot_inputs = mcmc_plot_inputs)

#GELMAN RUBIN

#MAKE FUNCTION;!!!!!!!!!!!
#****************************************************************
# PART 2: COMPARE CHAINS - GELMAN RUBIN MCMC DIAGNOSTIC :D 
#****************************************************************
get_gelman_rubin(mcmc_ssi_adp_canadaX6, mcmc_ssi_adp_canadaX7)

#a param
#lista = mcmc.list(as.mcmc(mcmc_ssi_da2$a_vec), as.mcmc(mcmc_ssi_da3$a_vec))
gelman_a = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_adp_canadaX6$a_vec), as.mcmc(mcmc_ssi_adp_canadaX7$a_vec)))
gelman_a

#b param
gelman_b = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_adp_canadaX6$b_vec), as.mcmc(mcmc_ssi_adp_canadaX7$b_vec)))
gelman_b

#c param
gelman_c = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_adp_canadaX6$c_vec), as.mcmc(mcmc_ssi_adp_canadaX7$c_vec)))
gelman_c

#r0 vec
gelman_r0 = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_adp_canadaX6$r0_vec), as.mcmc(mcmc_ssi_adp_canadaX7$r0_vec)))
gelman_r0

#Print
gelman_a
gelman_b
gelman_c
gelman_r0

