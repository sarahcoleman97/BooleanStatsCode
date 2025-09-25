setwd("~/Library/CloudStorage/Box-Box/BooleanStats/BooleanStatsCode/exampleAnalysis")

library(lme4) 
library(multcomp) 
library(dplyr) 
library(broom) 

# Read in data, make sure input group is treated as factor
data <- read.csv('data.csv') |> 
  mutate(group = group |> as.factor())

# Fit model
data.lmer <- lmer(Results ~ on + (1|group), data = data) 

# Hypothesis test
data.test <- glht(data.lmer,linfct=c('on==0')) |>
  tidy(conf.int = TRUE)

# Output test results
print(data.test)
