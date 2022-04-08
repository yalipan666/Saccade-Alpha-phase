# R scripts for paired-Ttest 

## set path
setwd("Z:/Saccade_AlphaPhase/Analyse_codes") 

## get needed packages
# effet size -- Cohen's d for t-test
install.packages("lsr")
library(lsr)

######################## first fixation ######################
## load in data
my_data <- read.csv("first_fixation.csv")

## Convert voriables as factors
my_data$wrd_frq <- factor(my_data$wrd_frq, 
                          levels = c(1, 2),
                          labels = c("low", "high"))
my_data$wrd_loc <- factor(my_data$wrd_loc, 
                          levels = c(1, 2),
                          labels = c("pre", "targ"))
my_data$subject <- factor(my_data$subject)
# pre-target: t(37)=0.054, p=0.958, d = 0.028
data_pre <- subset(my_data, wrd_loc == "pre")
t_pre <- t.test(duration ~ wrd_frq, data = data_pre, paired = TRUE, alternative = "two.sided")
t_pre
# cohen's d
pre_low <- subset(data_pre, wrd_frq == "low")
pre_high <- subset(data_pre, wrd_frq == "high")
cohensD(pre_low$duration, pre_high$duration, method = "paired")

# target: t(38)=6.939, p=2.975e-08, d = 1.111
data_targ <- subset(my_data, wrd_loc == "targ")
t_targ <- t.test(duration ~ wrd_frq, data = data_targ, paired = TRUE, alternative = "two.sided")
t_targ
# cohen's d
targ_low <- subset(data_targ, wrd_frq == "low")
targ_high <- subset(data_targ, wrd_frq == "high")
cohensD(targ_low$duration, targ_high$duration, method = "paired")


######################## power ######################
## load in data
my_data <- read.csv("power.csv")

## Convert voriables as factors
my_data$wrd_frq <- factor(my_data$wrd_frq, 
                          levels = c(1, 2),
                          labels = c("low", "high"))
my_data$subject <- factor(my_data$subject)
# power: t(37)=-1.189, p=0.242, d = 0.193
t_power <- t.test(power ~ wrd_frq, data = my_data, paired = TRUE, alternative = "two.sided")
t_power
# cohen's d
pre_low <- subset(my_data, wrd_frq == "low")
pre_high <- subset(my_data, wrd_frq == "high")
cohensD(pre_low$power, pre_high$power, method = "paired")
