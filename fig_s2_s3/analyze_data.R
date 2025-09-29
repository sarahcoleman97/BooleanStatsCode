setwd("~/Library/CloudStorage/Box-Box/BooleanStats/BooleanStatsCode/figs3")

set.seed(2024)
library(lme4)
library(dplyr)
library(broom)
library(multcomp)

source("simulate_data.R")

resim_data <- TRUE
reanal_data <- TRUE

# Define simulation conditions
sample.size <- c(3,5,30)
new.cond <- read.csv("simulation_conditions.csv",row.names=1,header=TRUE)
conditions <- new.cond |> row.names()
conditions_logics <- c(conditions[1:8], "OR") # specific to this case

# Define logics
and.logic <- c(0, 0, 0, 1)
or.logic <- c(0, 1, 1, 1)
nand.logic <- c(1, 1, 1, 0)
nor.logic <- c(1, 0, 0, 0)
xnor.logic <- c(1, 0, 0, 1)
xor.logic <- c(0, 1, 1, 0)
imply.logic <- c(1,0,1,1)
nimply.logic <- c(0,1,0,0)

# Order logics
logic_list <- list('AND' = and.logic,
                   'OR' = or.logic,
                   'NAND' = nand.logic,
                   'NOR' = nor.logic,
                   'XNOR' = xnor.logic,
                   'XOR' = xor.logic,
                   'IMPLY' = imply.logic,
                   'NIMPLY' = nimply.logic)

# Simulating data
if (resim_data){
  for (j in 1:length(sample.size)){
    n <- sample.size[j]
    for (i in 1:length(conditions)){
      condition <- conditions[[i]]
      logic <- logic_list[[conditions_logics[[i]]]]
      results <- simulate_data(new.cond[i,], 
                               n, 
                               logic, 
                               paste0('sim_data/',
                                      condition),
                               write.file = TRUE)
    }
  }
}

# Analyzing data
if (reanal_data){
  results_df <- lapply(conditions, function(condition){
    nested_df <- lapply(sample.size, function(n){
      fileName <- paste0('sim_data/',
                         condition,
                         '_sample_size_',
                         n,'_simulated_data.csv')
      working.df <- read.csv(fileName, header=TRUE, row.names = 1)
      working.lmer <- lmer(Results ~ on + (1|group), data = working.df)
      Singular <- isSingular(working.lmer)
      if (Singular){ # refit without random effects
        working.lmer <- lm(Results ~ on, data = working.df)
      }
      working.test <- glht(working.lmer,linfct=c('on==0')) |>
        tidy(conf.int = TRUE) |> 
        mutate(n = n,
               isSingular = Singular)
    }) |> bind_rows() |> 
      mutate(Condition = condition)
  }) |> bind_rows() |> 
    dplyr::select(Condition, n, isSingular, estimate, conf.low, conf.high, adj.p.value)
  write.csv(results_df, 'figss3stats.csv', row.names = F)
}

