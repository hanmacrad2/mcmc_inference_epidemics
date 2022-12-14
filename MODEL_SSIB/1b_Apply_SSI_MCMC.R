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

#*********
#DATA SIMULATION PARAMS
num_days = 50
shape_g = 6; scale_g = 1 #Infectious pressure (lambda) - gamma params

#SSI specific (*TO DO: DESIGN OF EXPERIMENTS FOR PARAM COMBINATIONS)
model_typeX = 'SSI'; 
aX = 0.8; bX = 0.1; cX = 10 
true_r0 = aX + bX*cX
model_params_true = list(m1 = aX, m2 = bX, m3 = cX, true_r0 = true_r0)

#*********
#MCMC PARAMS  
n_mcmc = 500 #100000 #100 #10000 #100000 #10000 #1000

#MCMC STARTING POINTS (RANDOM GUESS)
mod_start_points = list(m1 = 2, m2 = 0.05, m3 = 15, true_r0 = 1) #Random guess before finding paper. Try m3 = 8

#SIGMA *?
sigma = list(sigma1 = 0.4*mod_start_points$m1, sigma2 = 0.3*mod_start_points$m2 , #Acc rate too big -> Make sigma bigger. 
             sigma3 = 0.5*mod_start_points$m3, sigma4 =  0.85*mod_start_points$m3) #Acc rate too small -> make sigma smaller

mcmc_inputs = list(n_mcmc = n_mcmc, sigma = sigma, 
                   mod_start_points = mod_start_points)
seed_count = 3
mcmc_plot_inputs = list(n_mcmc = n_mcmc,
                        mod_start_points = mod_start_points, model_params_true = model_params_true,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = 'SSI', TYPEX = 'Individuals',
                        seed_count = seed_count)

# #SIGMA *?
# sigma_a = 0.4*aX; sigma_b = 1.0*bX #0.1 #SHOULD SIGMA BE DEFINED MORE RIGOROUS
# sigma_c = 0.85*cX; sigma_bc = 1.5*cX
# sigma = list(sigma_a = sigma_a, sigma_b = sigma_b, #Acc rate too big -> Make sigma bigger. 
#              sigma_c = sigma_c, sigma_bc = sigma_bc) #Acc rate too small -> make sigma smaller


#****************************************************************
#DATASET - GENERATED USING SIMULATION FUNCTIONS
#****************************************************************
seed_count = 3 #seed_count = seed_count + 1 #print(paste0('i mcmc = ', i))
set.seed(seed_count)
sim_data = simulation_super_spreaders(num_days, shape_g, scale_g, aX, bX, cX)

#PLOTS
par(mfrow=c(2,1))
non_ss = sim_data[[1]]
plot.ts(non_ss, ylab = 'Daily Infections count', main = 'Non Super-Spreaders' )
ss = sim_data[[2]]
plot.ts(ss, ylab = 'Daily Infections count', main = 'Super-Spreaders')

#Total
sim_dataX = non_ss + ss
plot.ts(sim_dataX, ylab = 'Daily Infections count', main = paste0('Total count - Super-Spreaders Model, R0 = ', true_r0))

#****************************************************************
#DATASET - SUPER-SPREADING EVENTS 20/80 RULE
#****************************************************************

#Params
alphaX = 0.8 #Without ss event, ~r0.
betaX = 0.2; gammaX = 4
true_r0 = alphaX + betaX*gammaX
model_params_true_se = list(m1 = alphaX, m2 = betaX, m3 = gammaX, true_r0 = true_r0)
#Epidemic data (SSE data also good)
par(mfrow = c(2,1))
sim_data_sse2 = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX,
                                 FLAG_seperate = TRUE)
non_sse2 = sim_data_sse2[[1]]
plot.ts(non_sse2, ylab = 'Daily Infections count', main = 'Non SSE' )
sse2 = sim_data_sse2[[2]]
plot.ts(sse2, ylab = 'Daily Infections count', main = 'SSE')

plot.ts(sim_data_sse, ylab = 'Daily infection count',
        main = paste('Super-Spreading Events Model - Daily Infections. R0 = ', true_r0)) #,
                     #expression(alpha), ':', alphaX, 'b:', betaX, ':', gammaX))

#****************************************************************
# I. APPLY MCMC SSI MODEL   
#***************************************************************

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi = MCMC_SSI(sim_data, mcmc_inputs = mcmc_inputs)

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi$time_elap = time_elap

#PLOT RESULTs
df_results_da = PLOT_MCMC_GRID(sim_data, mcmc_ssi,
                               mcmc_plot_inputs = mcmc_plot_inputs,
                               FLAGS_LIST = list(DATA_AUG = FALSE, BCA_TRANSFORM = TRUE,
                                                 PRIOR = TRUE,  JOINT = TRUE,
                                                 B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,  RJMCMC = FALSE))

#****************************************************************
# II. APPLY MCMC SSI MODEL + DATA AUGMENTATION  
#***************************************************************

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_da = MCMC_SSI(sim_data, mcmc_inputs = mcmc_inputs)

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_da$time_elap = time_elap

#PLOT RESULTS
mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points, model_params_true = model_params_true,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = 'SSI', TYPEX = 'Individuals',
                        seed_count = seed_count)

df_results_da = PLOT_MCMC_GRID(sim_data, mcmc_ssi_da,
                               mcmc_plot_inputs = mcmc_plot_inputs,
                               FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                 PRIOR = TRUE,  JOINT = TRUE,
                                                 B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,  RJMCMC = FALSE,
                                                 FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = FALSE))

#****************************************************************
# II. APPLY MCMC SSI MODEL + DATA AUGMENTATION + exp(0.1) prior on c
#***************************************************************

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_da3 = MCMC_SSI(sim_data, mcmc_inputs = mcmc_inputs,
                        FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                           PRIOR = TRUE, JOINT = TRUE,
                           B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_da3$time_elap = time_elap

#PLOT RESULTS
mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points, model_params_true = model_params_true,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = 'SSI', TYPEX = 'Individuals',
                        seed_count = seed_count)

df_results_da = PLOT_MCMC_GRID(sim_data, mcmc_ssi_da3,
                               mcmc_plot_inputs = mcmc_plot_inputs,
                               FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                 PRIOR = TRUE,  JOINT = TRUE,
                                                 B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = FALSE,  RJMCMC = FALSE,
                                                 FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = FALSE))

#****************************************************************
# II. APPLY MCMC SSI MODEL + DATA AUGMENTATION + exp(1) prior on b, 1 + exp(0.1) prior on c 
#***************************************************************

#DATA
sim_data4 = sim_data
sim_data4[[1]] = sim_data4[[1]] + sim_data4[[2]]
sim_data4[[2]] = rep(1, length(sim_data3[[1]]))
#sim_data4[[1]] = sim_data4[[1]] - sim_data4[[2]] #sim_data4[[1]] = abs()

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_da4 = MCMC_SSI(sim_data4, mcmc_inputs = mcmc_inputs,
                        FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                          PRIOR = TRUE, JOINT = TRUE,
                                          B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_da4$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        model_params_true = model_params_true,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = 'SSI', TYPEX = 'Individuals',
                        seed_count = seed_count)

df_results_da4 = PLOT_MCMC_GRID(sim_data, mcmc_ssi_da4,
                               mcmc_plot_inputs = mcmc_plot_inputs,
                               FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                 PRIOR = TRUE,  JOINT = TRUE,
                                                 B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE,  RJMCMC = FALSE,
                                                 FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = FALSE))

#****************************************************************
# III. APPLY MCMC SSI MODEL + NON-SS EXTREME CASE
#***************************************************************

#DATA PREP
FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                  PRIOR = TRUE, JOINT = TRUE,
                  B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE)

sim_data2 = sim_data
FLAG_NS_DATA_AUG = TRUE
if(FLAG_NS_DATA_AUG){ 
  sim_data2[[1]] = sim_data2[[1]] + sim_data2[[2]]
  sim_data2[[2]] = rep(0, length(sim_data2[[2]]))
  print(sim_data2[1])
  print(sim_data2[2])
}

#**********
#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))

mcmc_ssi_da2 = MCMC_SSI(sim_data2, mcmc_inputs = mcmc_inputs,
                           FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                             PRIOR = TRUE, JOINT = TRUE,
                                             B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE))
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_da2$time_elap = time_elap

#************
#PLOT RESULTS
mcmc_plot_inputs = list(n_mcmc = n_mcmc, model_params = model_params, mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = 'SSI', TYPEX = 'Individuals',
                        seed_count = seed_count, x0 = 1)

df_results2 = PLOT_MCMC_GRID(sim_data2, mcmc_ssi_da2,
                             mcmc_plot_inputs = mcmc_plot_inputs,
               FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                 PRIOR = TRUE,  JOINT = TRUE,
                                 B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,  RJMCMC = FALSE,
                                 FLAG_NS_DATA_AUG = TRUE, FLAG_SS_DATA_AUG = FALSE))

#****************************************************************
#*
# IV. APPLY MCMC SSI MODEL + SS EXTREME CASE
#*
#***************************************************************
#DATA PREP
FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                  PRIOR = TRUE, JOINT = TRUE,
                  B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                  FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = TRUE)

sim_data3 = sim_data
if(FLAGS_LIST$FLAG_SS_DATA_AUG){ 
  sim_data3[[2]] = sim_data3[[1]] + sim_data3[[2]]
  sim_data3[[1]] = rep(0, length(sim_data3[[1]]))
  print(sim_data3[1])
  print(sim_data3[2])
}
#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))

#REGUALR SIM DATA?? WRONG!
mcmc_ssi_da3 = MCMC_SSI(sim_data, mcmc_inputs = mcmc_inputs,
                           FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                             PRIOR = TRUE, JOINT = TRUE,
                                             B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,  RJMCMC = FALSE,
                                             FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = TRUE))

end_time = Sys.time()
mcmc_ssi_da3$time_elap = time_elap

#PLOT RESULTS
df_results3b = PLOT_MCMC_GRID(sim_data3, mcmc_ssi_da3,
                             mcmc_plot_inputs = mcmc_plot_inputs,
                             FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                               PRIOR = TRUE,  JOINT = TRUE,
                                               B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,  RJMCMC = FALSE,
                                               FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = TRUE))

#****************************************************************
#*
# IV. APPLY MCMC SSI MODEL + SSE DATA
#*
#***************************************************************
#DATA PREP

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))

mcmc_ssi_sse2 = MCMC_SSI(sim_data_sse2, mcmc_inputs = mcmc_inputs,
                        FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                          PRIOR = TRUE, JOINT = TRUE,
                                          B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,  RJMCMC = FALSE,
                                          FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = FALSE))

end_time = Sys.time()
mcmc_ssi_sse2$time_elap = time_elap

#PLOT RESULTS
model_params_true = list(m1 = 0.8, m2 = 0.2, m3 = 4, true_r0 = true_r0)
mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        model_params_true = model_params_true, mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = 'SSI', TYPEX = 'Individuals',
                        seed_count = seed_count)
df_results5 = PLOT_MCMC_GRID(sim_data_sse2, mcmc_ssi_sse2,
                              mcmc_plot_inputs = mcmc_plot_inputs,
                              FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                PRIOR = TRUE,  JOINT = TRUE,
                                                B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,  RJMCMC = FALSE,
                                                FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = FALSE))


#****************************************************************
# II. APPLY MCMC SSI MODEL + DATA AUGMENTATION  
#***************************************************************

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_da = MCMC_SSI(sim_data, mcmc_inputs = mcmc_inputs)

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_da$time_elap = time_elap

#PLOT RESULTS
mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points, model_params_true = model_params_true,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = 'SSI', TYPEX = 'Individuals',
                        seed_count = seed_count)

df_results_da = PLOT_MCMC_GRID(sim_data, mcmc_ssi_da,
                               mcmc_plot_inputs = mcmc_plot_inputs,
                               FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                 PRIOR = TRUE,  JOINT = TRUE,
                                                 B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,  RJMCMC = FALSE,
                                                 FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = FALSE))

#****************************************************************
# II. APPLY ADAPTIVE MCMC SSI MODEL + DATA AUGMENTATION + exp(1) prior on b, 1 + exp(0.1) prior on c 
#***************************************************************

#DATA
sim_data6 = sim_data
sim_data6[[1]] = sim_data6[[1]] + sim_data6[[2]]
sim_data6[[2]] = rep(1, length(sim_data6[[1]]))
#sim_data4[[1]] = sim_data4[[1]] - sim_data4[[2]] #sim_data4[[1]] = abs()

mcmc_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                   alpha_star = 0.4)
#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_da6 = MCMC_SSI_ADAPTIVE(sim_data6, mcmc_inputs = mcmc_inputs,
                        FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                          PRIOR = TRUE, JOINT = TRUE,
                                          B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_da6$time_elap = time_elap

mcmc_plot_inputs = list(n_mcmc = n_mcmc, mod_start_points = mod_start_points,
                        model_params_true = model_params_true,
                        mod_par_names = c('a', 'b', 'c'),
                        sigma = sigma, model_typeX = 'SSI', TYPEX = 'Individuals',
                        seed_count = seed_count)

df_results_da6 = PLOT_MCMC_GRID(sim_data, mcmc_ssi_da6,
                                mcmc_plot_inputs = mcmc_plot_inputs,
                                FLAGS_LIST = list(DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                  PRIOR = TRUE,  JOINT = TRUE,
                                                  B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE,  RJMCMC = FALSE,
                                                  FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = FALSE))



#****************************************************************
# PART 2: COMPARE CHAINS - GELMAN RUBIC MCMC DIAGNOSTIC :D 
#****************************************************************
#a param
#lista = mcmc.list(as.mcmc(mcmc_ssi_da2$a_vec), as.mcmc(mcmc_ssi_da3$a_vec))
gelman_a = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_da2$a_vec), as.mcmc(mcmc_ssi_da3$a_vec)))
gelman_a

#b param
gelman_b = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_da2$b_vec), as.mcmc(mcmc_ssi_da3$b_vec)))
gelman_b

#c param
gelman_c = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_da2$c_vec), as.mcmc(mcmc_ssi_da3$c_vec)))
gelman_c

#r0 vec
gelman_r0 = gelman.diag(mcmc.list(as.mcmc(mcmc_ssi_da2$r0_vec), as.mcmc(mcmc_ssi_da3$r0_vec)))
gelman_r0
