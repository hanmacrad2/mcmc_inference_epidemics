#LIBRARIES
library(EpiEstim)
source("~/GitHub/SuperSpreadingEpidemicsMCMC/R/EPI_FUNCTIONS.R")

#Data

#Flu
data("Flu2009")
Flu2009$incidence[,2]
Flu2009$si_distr

res <- estimate_R(incid = Flu2009$incidence[, 2],
                  method = "non_parametric_si",
                  config = make_config(list(si_distr = Flu2009$si_distr)))

plot(res)

res2 <- estimate_R(incid = Flu2009$incidence[, 2],
                  method = "non_parametric_si",
                  config = make_config(list(si_distr = Flu2009$si_distr)))

plot(res2)

#lambda_vec = get_lambda(canadaX)
#lambda_vec = c(0, lambda_vec)

#SARS CANADA (Our model)
infectivity_vec = dgamma(0:15, shape = shape_gamma, scale = scale_gamma)
infectivity_vec <- infectivity_vec/sum(infectivity_vec)

#a = 1 - sum(dgamma(0:10000, shape = shape_gamma, scale = scale_gamma))
#infectivitiy_vec = c(infectivitiy_vec, a)
plot.ts(infectivitiy_vec)

#ESTIMATE R
result_II = estimate_R(canadaX, config = make_config(list(si_distr = infectivity_vec)))

plot(result_II)

sum(lambda_vec)
