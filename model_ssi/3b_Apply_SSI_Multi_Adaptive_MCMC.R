#****************************************************************
#*
# APPLY ADAPTIVE SHAPING
#*
#****************************************************************
library(coda)
library(MASS)
setwd("~/GitHub/epidemic_modelling")
source("epidemic_functions.R")
source("plot_functions.R")
source("helper_functions.R")


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
# MCMC INPUTS
#****************************************************************
n_mcmc = 500000 #100000 #0 #0 #500000
mod_start_points = list(m1 = 0.5, m2 = 0.02, m3 = 20)
mcmc_inputs = list(n_mcmc = n_mcmc, mod_start_points = c(0.5, 0.02, 20), burn_in_pc = 0.05,
                   dim = 3, alpha_star = 0.4, v0 = 100, vec_min = c(0,0,1), thinning_factor = 10)

#****************************************************************
# N0 1 START NO SS
#****************************************************************

#MOD START POINTS
n_mcmc = 500000 #100000 #0 #0 #500000
#SAME STARTING POINTS AS FINAL POINTS IN LAST MCMC RUN
mod_start_points = list(m1 = 0.72, m2 = 0.0038, m3 = 22)
mcmc_inputs = list(n_mcmc = n_mcmc, mod_start_points = c(0.5, 0.02, 20), burn_in_pc = 0.05,
                   dim = 3, alpha_star = 0.4, v0 = 100, vec_min = c(0,0,1), thinning_factor = 10)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_multiX3 = MCMC_SSI_MULTI_ADADPT(sim_data_canadaX1, mcmc_inputs) #, FLAGS_LIST = list(DATA_AUG = TRUE))
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_multiX3$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count, thinning_factor = 10)

df_results_M3 = PLOT_MCMC_SSI_GRID_REAL_DATA(sim_data_canadaX1, mcmc_multiX3,
                                             mcmc_plot_inputs = mcmc_plot_inputs,
                                             FLAGS_LIST = list(DATA_AUG = FALSE,
                                                               PRIOR = TRUE, B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                                               FLAG_SSI = TRUE, RJMCMC = FALSE,
                                                               BURN_IN = TRUE, FLAG_ADAPTIVE = FALSE,
                                                               MULTI_ALG = TRUE, THIN = TRUE))
#PLOTS ADDITIONAL
PLOT_MCMC_MEANS(mcmc_multiX3,
                mcmc_plot_inputs = mcmc_plot_inputs,
                FLAGS_LIST = list(DATA_AUG = TRUE,
                                  PRIOR = TRUE, B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                  FLAG_SSI = FALSE, RJMCMC = FALSE,
                                  BURN_IN = TRUE,
                                  FLAG_ADAPTIVE = FALSE, MULTI_ALG = TRUE))

# PLOT_SIGMA_ADADPTIVE(mcmc_ssi_adp_canadaX8, mcmc_inputs)

#****************************************************************
# N0 2 SS; 1, 0
#****************************************************************
#mcmc_multiX1 = MCMC_SSI_MULTI_ADADPT(sim_data_canadaX1, mcmc_inputs, FLAGS_LIST = list(DATA_AUG = TRUE))

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_multiX4 = MCMC_SSI_MULTI_ADADPT(sim_data_canadaX2, mcmc_inputs)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_multiX4$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = typeX, TYPEX = '',
                        seed_count = seed_count)

#*No prior
df_results_M4 = PLOT_MCMC_SSI_GRID_REAL_DATA(sim_data_canadaX2, mcmc_multiX4,
                                             mcmc_plot_inputs = mcmc_plot_inputs,
                                             FLAGS_LIST = list(DATA_AUG = FALSE,
                                                               PRIOR = TRUE, B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                                               FLAG_SSI = TRUE, RJMCMC = FALSE,
                                                               BURN_IN = TRUE, FLAG_ADAPTIVE = FALSE,
                                                               MULTI_ALG = TRUE))

#PLOTS ADDITIONAL
PLOT_MCMC_MEANS(mcmc_multiX4,
                mcmc_plot_inputs = mcmc_plot_inputs,
                FLAGS_LIST = list(DATA_AUG = TRUE,
                                  PRIOR = TRUE, B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                  FLAG_SSI = FALSE, RJMCMC = FALSE,
                                  BURN_IN = TRUE,
                                  FLAG_ADAPTIVE = FALSE, MULTI_ALG = TRUE))

#PLOT LAMBDA
#PLOT_SIGMA_ADADPTIVE(mcmc_ssi_adp_canadaX8, mcmc_inputs)

#****************************************************************
# PART 2: COMPARE CHAINS - GELMAN RUBIN MCMC DIAGNOSTIC :D
#****************************************************************
#get_gelman_rubin(mcmc_ssi_adp_canadaX6, mcmc_ssi_adp_canadaX7)

#a param
#lista = mcmc.list(as.mcmc(mcmc_ssi_da2$a_vec), as.mcmc(mcmc_ssi_da3$a_vec))
gelman_a = gelman.diag(mcmc.list(as.mcmc(mcmc_multiX3$x_matrix[,1]), as.mcmc(mcmc_multiX4$x_matrix[,1])))

#b param
gelman_b = gelman.diag(mcmc.list(as.mcmc(mcmc_multiX3$x_matrix[,2]), as.mcmc(mcmc_multiX4$x_matrix[,2])))

#c param
gelman_c = gelman.diag(mcmc.list(as.mcmc(mcmc_multiX3$x_matrix[,3]), as.mcmc(mcmc_multiX4$x_matrix[,3])))

#r0 vec
#gelman_r0 = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_adp_canadaX6$r0_vec), as.mcmc(mcmc_ssi_adp_canadaX7$r0_vec)))

#Print
gelman_a
gelman_b
gelman_c
#gelman_r0
