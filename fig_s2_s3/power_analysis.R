setwd("~/Library/CloudStorage/Box-Box/BooleanStats/BooleanStatsCode/fig_s2_s3")

set.seed(2024) 

library(lme4)
library(effects)
library(broom)
library(dplyr)
library(multcomp)

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

# Define simulation conditions
sim.means <- read.csv("simulation_conditions.csv",row.names=1,header=TRUE)
conditions <- sim.means |> row.names()
samples <- 2:30
variance.values <- 1:10 * 0.1
n.mc <- 1e4

# Make a matrix of all the condtions I want to test
colNames <- c('Condition', 'Variance', 'Sample.Size')
combs <- expand.grid(conditions, variance.values, samples,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) |> 
  rename_with(~colNames) |> 
  mutate(Logic = if_else(Condition == 'OR1','OR',Condition)) |> 
  as.matrix()

# Turn into list and get sim means for each row, and logic
entryNames <- c('Condition', 'Variance', 'Sample.Size','Logic')
comb_list_pt1 <- split(combs, row(combs))

# Define cores for parallel stuff 
n_cores <- parallel::detectCores() - 2

colNames <- c('Variance', 'Sample.Size')
combs <- expand.grid(variance.values, samples,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) |> 
  rename_with(~colNames) |> 
  mutate(Condition = 'NAND', Logic = 'NAND') |> 
  dplyr::select(Condition, Variance, Sample.Size, Logic) |> 
  as.matrix()
comb_list_pt1 <- split(combs, row(combs))

condition_list <- parallel::mclapply(comb_list_pt1, function(row){
  
  # Make each row a list, name it
  row <- row |> as.list() 
  names(row) <- entryNames 
  
  # Make sure classes are correct
  row[['Variance']] <- row[['Variance']] |> as.numeric()
  row[['Sample.Size']] <- row[['Sample.Size']] |> as.numeric()
  
  # Add additional components
  Means <- sim.means[row[['Condition']],] # get the sim.means
  Logic <- logic_list[[row[['Logic']]]]
  row_final <- list('Row' = row, 'Means' = Means, 'Logic' = Logic)
}, mc.cores = n_cores)

# Simulate data (doing mclapply inside a for list due to memory issues)

pb <- txtProgressBar(min = 0, max = length(condition_list), style = 3)
prog <- 0

# Final loop (takes a while)
for (ind_list in condition_list){ # doing this as a for loop to reduce overwhelm
  
  # Get conditions
  condition <- ind_list[['Row']][['Condition']]
  variance <- ind_list[['Row']][['Variance']] |> as.numeric()
  sample.size <- ind_list[['Row']][['Sample.Size']] |> as.integer()
  working.params.df <- ind_list[['Means']]
  logic <- ind_list[['Logic']]
  
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
            paste0('beta_montecarlo/',condition,
                   '_var_',variance,
                   '_n_',sample.size,
                   '.csv'),
            row.names = F)
  
  prog <- prog + 1
  setTxtProgressBar(pb, prog)
}

# Pull out the select betas to make it easy to copy into Prism for plotting
samples_select <- c(3, 5, 30)

for (condition in conditions){
  files <- paste0('beta_montecarlo/', condition, 
                  '_var_0.1_n_',samples_select, '.csv')
  
  df <- lapply(samples_select, function(sample) {
    f <- paste0('beta_montecarlo/', condition, 
                '_var_0.1_n_',sample, '.csv') 
    newColName <- paste0('Beta_n_', sample)
    dat <- read.csv(f) |> 
      dplyr::select(Beta) |> 
      rename_with(~newColName)
  }) |> bind_cols()
  
  write.csv(df, paste0('sim_data/', condition, 
                       '_betas_', 'n_3_5_30.csv'),
            row.names = FALSE)
}
