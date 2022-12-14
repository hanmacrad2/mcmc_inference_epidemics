#DATA
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_individual_nu/"

#PLOT
plot.new()
par(mfrow=c(3,3))

#ETA
eta8 = 'eta_matrix_8.rds'
eta_df = readRDS(paste0(OUTER_FOLDER, '/', eta8))

#Function

plot_eta <- function(eta_df, day) {
  
  #DAY
  eta_dfX = eta_df[, day]
  
  plot.ts(eta_dfX,  ylab = paste0('Eta day', day), #ylim = m3_lim,
          main = paste0('Eta day', day),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
}

#Plot
plot_eta(eta_df, 1)
plot_eta(eta_df, 6)
plot_eta(eta_df, 12)
plot_eta(eta_df, 16)
plot_eta(eta_df, 19)
plot_eta(eta_df, 23)
plot_eta(eta_df, 24)
plot_eta(eta_df, 26)
plot_eta(eta_df, 49)
plot_eta(eta_df, 50)

plot_eta(eta_df, 80)
plot_eta(eta_df, 79)
plot_eta(eta_df, 81)
