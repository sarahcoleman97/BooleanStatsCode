setwd("~/Library/CloudStorage/Box-Box/BooleanStats/BooleanStatsCode/fig_3")

set.seed(2024)
library(lme4)
library(dplyr)
library(broom)
library(multcomp)

resim_data <- FALSE
reanal_data <- FALSE

source("simulate_data.R")

# Read in condition information
new.cond <- read.csv('k_means_centroids.csv',header=TRUE,row.names=1)
groups <- row.names(new.cond)
sample.size <- 3
or.logic <- c(0, 1, 1, 1)

# Simulate data
if (resim_data){
  for (i in 1:length(groups)){ 
   working.row.name = paste0('sim_data/group_',groups[i])
   or.results.i <- simulate_data(new.cond[i,], 
                                 sample.size, 
                                 or.logic, 
                                 working.row.name, 
                                 write.file=TRUE)
  }
}

# Analyze data 
if (reanal_data){
  results_df <- lapply(groups, function(group){
    fileName <- paste0('sim_data/group_',
                       group,
                       '_sample_size_',
                       sample.size,'_simulated_data.csv')
    working.df <- read.csv(fileName, header=TRUE, row.names = 1)
    working.lmer <- lmer(Results ~ on + (1|group), data = working.df)
    Singular <- isSingular(working.lmer)
    if (Singular){ # refit without random effects
      working.lmer <- lm(Results ~ on, data = working.df)
    }
    working.test <- glht(working.lmer,linfct=c('on==0')) |>
      tidy(conf.int = TRUE) |> 
      mutate(Group = group,
             isSingular = Singular)
  }) |> bind_rows() |> 
    dplyr::select(Group, isSingular, estimate, conf.low, conf.high, adj.p.value)
  write.csv(results_df, 'fig3stats.csv', row.names = F)
}

# Formatting simulated data for Prism
n <- 3
group_mapping <- c('-/-', '+/-', '-/+','+/+')
sample_idx <- rep(1:n, times = length(group_mapping))

results_df <- lapply(groups, function(group){
  
  # Define filename
  fileName <- paste0('sim_data/group_',
                     group,
                     '_sample_size_',
                     sample.size,'_simulated_data.csv')
  
  # Read in dataframe
  df <- read.csv(fileName, header=TRUE, row.names = 1) |> 
    mutate(Group = group_mapping[group],
           Replicate = sample_idx) |> 
    dplyr::select(Results, Group, Replicate) |> 
    pivot_wider(names_from = Group, values_from = Results) |>
    dplyr::select(-Replicate)
  
  write.csv(df, 
            paste0('sim_data/group_', group, '_formatted_for_prism.csv'),
            row.names = F)
}) 


