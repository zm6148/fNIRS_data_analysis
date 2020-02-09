# Load eye-tracking data for one subject,
# mung it, plot it
library(ggplot2)
library(reshape2)
library(dplyr)
library(lme4)

# First remove all objects
rm(list=ls())

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Load the files
setwd("C:\\Users\\Matt\\Documents\\ListenLab\\Tools\\GCA")
load("demo_data.RData")
source("code_poly.R")
source("get_pvals.R")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Summarise for plotting
#
# First for individuals
data_sum_indiv <- demo_data %>%
  ungroup() %>%
  group_by(Participant, Hearing, Group, Context, Time) %>%
  summarise(pupil = mean(pupil, na.rm = TRUE),
            pupil_lpf = mean(pupil_lpf, na.rm = TRUE)) 

# Then whole groups, with variance
data_sum_group <- data_sum_indiv %>%
  ungroup() %>%
  group_by(Group, Hearing, Context, Time) %>%
  summarise(num_listeners = n(),
            pupil.se = sd(pupil, na.rm = TRUE)/sqrt(n()),
            pupil_lpf.se = sd(pupil_lpf, na.rm = TRUE)/sqrt(n()),
    pupil = round(mean(pupil, na.rm = TRUE),4),
            pupil_lpf = round(mean(pupil_lpf, na.rm = TRUE),4))

# Basic plot
ggplot(data_sum_group)+
  aes(x = Time, y = pupil_lpf, color = Context)+
  geom_ribbon(aes(ymin = pupil_lpf - pupil_lpf.se,
                  ymax = pupil_lpf + pupil_lpf.se,
                  fill = Context),
              alpha = 0.5, color = NA)+
  geom_line()+
  facet_grid(. ~ Group)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# establish analysis windows
window.part1.start <- -1000
window.part1.end <- 1000

window.part2.start <- 1000
window.part2.end <- 3000
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Add polynomials to windowed data frames

df_poly_win1 <- data_sum_indiv %>%
  filter(Time >= window.part1.start,
         Time <= window.part1.end) %>%
  code_poly(df = ., 
            predictor = "Time",
            poly.order = 3, 
            orthogonal = TRUE, 
            draw.poly = TRUE)

df_poly_win2 <- data_sum_indiv %>%
  filter(Time >= window.part2.start,
         Time <= window.part2.end) %>%
  code_poly(df = ., 
            predictor = "Time",
            poly.order = 3, 
            orthogonal = TRUE, 
            draw.poly = TRUE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Fit the first window
gca_model_win1 <- lmer(pupil_lpf ~ 
                         # Default effects
                         poly1 + poly2 + poly3 + 
                         Group + Context +
                         # two-way interactions
                         Group:Context +
                         poly1:Context + 
                         poly2:Context + 
                         poly3:Context +
                         poly1:Group + 
                         poly2:Group + 
                         poly3:Group +
                         # three-way interactions
                         poly1:Context:Group +
                         poly2:Context:Group +
                         poly3:Context:Group +
                         # Random effects
                         #-----------------------------------#
                         # these are the actual model effects
                         # (1+
                         #    poly1 + poly2 + 
                         #    poly3 + Context + Group +
                         #    Context:Group|Participant),
                         #-----------------------------------#
                         # this is a simplified random effect
                         # just to run quickly in this demo.
                         (1|Participant),
                         #-----------------------------------#
                       data=df_poly_win1)

model_win1_summary <- summary(gca_model_win1)
model_win1_summary

# Add p values
coefs_1 <- get_pvals(gca_model_win1)
coefs_1

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Verify model output on the plot
# add fitted valued to data frame
df_poly_win1$predicted <- fitted(gca_model_win1)

# summarize first
df_poly_win1_sum <- df_poly_win1 %>%
  ungroup() %>%
  group_by(Group, Hearing, Context, Time) %>%
  summarise(num_listeners = n(),
            pupil.se = sd(pupil, na.rm = TRUE)/sqrt(n()),
            pupil_lpf.se = sd(pupil_lpf, na.rm = TRUE)/sqrt(n()),
            pupil = round(mean(pupil, na.rm = TRUE),4),
            pupil_lpf = round(mean(pupil_lpf, na.rm = TRUE),4),
            predicted = mean(predicted))

# Plot model prediction on top of the real data
ggplot(df_poly_win1_sum)+
  aes(x = Time, y = pupil_lpf, color = Context)+
  geom_ribbon(aes(ymin = pupil_lpf - pupil_lpf.se,
                  ymax = pupil_lpf + pupil_lpf.se,
                  fill = Context),
              alpha = 0.5, color = NA)+
  # actual data
  geom_line()+
  # model prediction
  # first white line
  geom_line(aes(y = predicted, group = Context), color = "white", size = 2.2)+
  # dashed line
  geom_line(aes(y = predicted, group = Context), size = 1, linetype = "dashed")+
  facet_grid(. ~ Group)
