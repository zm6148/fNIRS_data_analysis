rm(list = ls())#_# remove all objects__
library(plyr)
library(dplyr)
library(rmcorr)
library(lme4)
library(ggplot2)
library(reshape2)

# R_STG
full_data <- read.csv(file="C:/Users/mz86/Desktop/Fnirs data and analysis/fnirs_R_analysis_antje_using_d_prime experiment new 1/active vs passive - R_STG_HbT.csv",header = TRUE,sep=",")  
#rmcorr(full_data$SID, full_data$dprime, full_data$RAC, full_data)
R_STG <- rmcorr(participant = SID, measure1 = dprime, measure2 = R_STG, dataset = full_data)
library(RColorBrewer)
blueset <- brewer.pal(8, 'Spectral')
pal <- colorRampPalette(blueset)
plot(R_STG, overall = FALSE, palette = pal, overall.col = 'black', lty = 2, xlab = "d prime", ylab = "R_STG",ylim=c(-0.2,1.5))

# L_STG
full_data2 <- read.csv(file="C:/Users/mz86/Desktop/Fnirs data and analysis/fnirs_R_analysis_antje_using_d_prime experiment new 1/active vs passive - L_STG_HbT.csv",header = TRUE,sep=",")  
my.rmc <- rmcorr(participant = SID, measure1 = dprime, measure2 = L_STG, dataset = full_data2)
L_STG <-my.rmc
library(RColorBrewer)
blueset <- brewer.pal(8, 'Spectral')
pal <- colorRampPalette(blueset)
plot(L_STG, overall = FALSE, palette = pal, overall.col = 'black', lty = 2, xlab = "d prime", ylab = "L_STG", ylim=c(-0.2,1.5))

# R_cIFS
full_data3 <- read.csv(file="C:/Users/mz86/Desktop/Fnirs data and analysis/fnirs_R_analysis_antje_using_d_prime experiment new 1/active vs passive - R_cIFS_HbT.csv",header = TRUE,sep=",")  
my.rmc <- rmcorr(participant = SID, measure1 = dprime, measure2 = R_cIFS, dataset = full_data3)
R_cIFS <-my.rmc
library(RColorBrewer)
blueset <- brewer.pal(8, 'Spectral')
pal <- colorRampPalette(blueset)
plot(R_cIFS, overall = FALSE, palette = pal, overall.col = 'black', lty = 2, xlab = "d prime", ylab = "R_cIFS")

#L_cIFS
full_data4 <- read.csv(file="C:/Users/mz86/Desktop/Fnirs data and analysis/fnirs_R_analysis_antje_using_d_prime experiment new 1/active vs passive - L_cIFS_HbT.csv",header = TRUE,sep=",")  
my.rmc <- rmcorr(participant = SID, measure1 = dprime, measure2 = L_cIFS, dataset = full_data4)
L_cIFS <-my.rmc
library(RColorBrewer)
blueset <- brewer.pal(8, 'Spectral')
pal <- colorRampPalette(blueset)
plot(L_cIFS, overall = FALSE, palette = pal, overall.col = 'black', lty = 2, xlab = "d prime", ylab = "L_cIFS")


print(R_STG)
print(L_STG)
print(R_cIFS)
print(L_cIFS)



