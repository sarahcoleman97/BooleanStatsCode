setwd("~/Library/CloudStorage/Box-Box/BooleanStats/BooleanStatsCode/fig_s6_s7_s8_s9")

set.seed(2024) 

library(lme4)
library(effects)
library(broom)
library(dplyr)
library(tidyr)
library(multcomp)

# Define simulation groups
sim.means <- read.csv("k_means_centroids.csv",row.names=1,header=TRUE)
groups <- sim.means |> row.names()
samples <- 2:30
variance.values <- 1:10 * 0.1
n.mc <- 1e4

# Make a matrix of all the condtions I want to test
colNames <- c('Group', 'Variance', 'Sample.Size')
combs <- expand.grid(groups, variance.values, samples,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) |> 
  rename_with(~colNames) |> 
  as.matrix()

# Turn into list and get sim means for each row, and logic
entryNames <- c('Group', 'Variance', 'Sample.Size')
comb_list_pt1 <- split(combs, row(combs))

# Define cores for parallel stuff 
n_cores <- parallel::detectCores() - 2

group_list <- parallel::mclapply(comb_list_pt1, function(row){
  
  # Make each row a list, name it
  row <- row |> as.list() 
  names(row) <- entryNames 
  
  # Make sure classes are correct
  row[['Variance']] <- row[['Variance']] |> as.numeric()
  row[['Sample.Size']] <- row[['Sample.Size']] |> as.numeric()
  
  # Add additional components
  Means <- sim.means[row[['Group']],] # get the sim.means
  row_final <- list('Row' = row, 'Means' = Means)
}, mc.cores = n_cores)

# Simulate data (doing mclapply inside a for list due to memory issues)

pb <- txtProgressBar(min = 0, max = length(group_list), style = 3)
prog <- 0

# Final loop (takes a while)
for (ind_list in group_list){ # doing this as a for loop to reduce overwhelm
  
  # Get groups
  Group <- ind_list[['Row']][['Group']]
  variance <- ind_list[['Row']][['Variance']] |> as.numeric()
  sample.size <- ind_list[['Row']][['Sample.Size']] |> as.integer()
  working.params.df <- ind_list[['Means']]
  logic <- c(0, 1, 1, 1)
  
  # Making factors
  group <- rep(c(1:4), each = sample.size) |> as.factor()
  on <- rep(logic, each = sample.size)
  
  # Simulate data
  df_list <- parallel::mclapply(1:n.mc, function(mc){ # for each mc
  
    # Simulating
    var <- variance * (working.params.df |> unlist())
    Results <- lapply(names(var), function(name){
      sim <- rnorm(sample.size, working.params.df[[name]], var[[name]])
    }) |> unlist()

    # Construct df
    working.df <- data.frame(Results, on, group) 
  }) 
  
  # Stats
  result_df <- parallel::mclapply(1:n.mc, function(mc){
    
    # Get df
    df <- df_list[[mc]]
    
    # Fit model, check for singularity
    working.lmer <- lmer(Results ~ on + (1|group), data = df)
    Singular <- isSingular(working.lmer)
    
    if (Singular){ # refit without random effects
      working.lmer <- lm(Results ~ on, data = df)
      beta.estimate <- coef(working.lmer)[['on']]
    } else {
      beta.estimate <- coef(working.lmer)[['group']][['on']][1]
    }
    
    # Do test
    answer <- glht(working.lmer,linfct=c('on==0')) |> tidy()
    p.val <- (answer |> pull(adj.p.value))[1]
    
    result <- data.frame('MonteCarlo' = mc,
                         'Beta' = beta.estimate, 
                         'p.value' = p.val,
                         'singular' = Singular)
  }, mc.cores = n_cores) |> bind_rows()
  
  write.csv(result_df, 
            paste0('beta_montecarlo/',Group,
                   '_var_',variance,
                   '_n_',sample.size,
                   '.csv'),
            row.names = F)
  
  prog <- prog + 1
  setTxtProgressBar(pb, prog)
}

# Power analysis figures (easy for plotting)
var_for_out <- rep(variance.values, each = length(samples))
sam_for_out <- rep(samples, times = length(variance.values))

p_val_thresh <- 0.05

group_df_list <- vector(mode = 'list', length = length(groups))
names(group_df_list) <- groups

# go through and read all the files
for (group in groups){
  files <- paste0('beta_montecarlo/', group,
                  '_var_',var_for_out,'_n_',sam_for_out, '.csv')
  
  df <- parallel::mclapply(1:length(files), function(idx2){
    file <- files[[idx2]]
    var <- var_for_out[[idx2]]
    sam <- sam_for_out[[idx2]]
    dat <- read.csv(file) |> 
      summarise(Power = sum(p.value < p_val_thresh)/n.mc) |> 
      mutate(Variance = var,
             Sample.Size = sam)
  }) |> bind_rows() |> 
    mutate(Condition = group)
  
  group_df_list[[group]] <- df
} 

# make a final list formatted for prism
final_df_list <- lapply(group_df_list, function(df){
  cond <- df |> pull(Condition) |> unique()
  df_pivot <- df |>
    pivot_wider(names_from = Variance,
                values_from = Power) |> 
    dplyr::select(Sample.Size, all_of(variance.values |> as.character()))
  write.csv(df_pivot, 
            paste0('power_analysis/',cond,'_power_for_plotting.csv'),
            row.names = FALSE)
})
