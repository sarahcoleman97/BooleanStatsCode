setwd("~/Library/CloudStorage/Box-Box/BooleanStats/BooleanStatsCode/exampleAnalysis")

library(lme4) 
library(multcomp) 
library(dplyr) 
library(broom) 
library(msm)

# Read in data
data <- read.csv('data.csv')
data <- data |> 
  mutate(group = group |> as.factor(),
         on = on |> as.factor())

# Fit model
mod <- lmer(Results ~ on + (1|group) - 1, data = data) 

# Check for singular fit
Singular <- isSingular(mod)
if (Singular){ # refit without random effects
  print('Model is singular, refitting without random effects')
  mod <- lm(Results ~ on, data = working.df)
}

# Estimate beta (ON - OFF)
glht_beta <- glht(mod,linfct=c('on1-on0==0')) # equiv to on1/on0 == 1
est_beta <- glht_beta |> summary()

# Estimate gamma (ON/OFF)
mean_gam <- fixef(mod) # fixed effects of model
vcov_gam <- vcov(mod) # vcov of model
est_gamma <- (mean_gam['on1']/mean_gam['on0']) |> unname() # estimate gamma
se_gamma <- deltamethod(~x2/x1, mean_gam, vcov_gam) # estimate gamma se

# Statistical tests involving gamma require an approximation of the associated 
# degrees of freedom. We refer the interested reader to articles such as 
# https://doi.org/10.1007/s11222-014-9488-7