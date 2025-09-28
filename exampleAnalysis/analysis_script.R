setwd("~/Library/CloudStorage/Box-Box/BooleanStats/BooleanStatsCode/exampleAnalysis")

library(lme4) 
library(multcomp) 
library(dplyr) 
library(broom) 

# Read in data, make sure input group is treated as factor
data <- read.csv('data.csv') |> 
  mutate(group = group |> as.factor())

# Fit model
working.lmer <- lmer(Results ~ on + (1|group), data = data) 

# Check for singular fit
Singular <- isSingular(working.lmer)
if (Singular){ # refit without random effects
  print('Model is singular, refitting without random effects')
  working.lmer <- lm(Results ~ on, data = working.df)
}

# Hypothesis test
data.test <- glht(working.lmer,linfct=c('on==0')) |>
  tidy(conf.int = TRUE)

# Output test results
print(data.test)
