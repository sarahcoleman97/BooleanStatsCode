setwd("~/Library/CloudStorage/Box-Box/BooleanStats/BooleanStatsCode/fig_s10")

library(lme4)
library(dplyr)
library(multcomp)
library(tidyr)
library(broom)

# Define column names
colNames <- c("off", "on_1", "on_2", "on_3")

# List all the data in figure 3
csv_files <- list.files(path = '../fig_3/sim_data/', pattern = '*_simulated_data') 

# Some stuff for dataframe reformatting
n <- 3
group_mapping <- c('-/-', '+/-', '-/+','+/+')
sample_idx <- rep(1:n, times = length(group_mapping))

# Analyze the data with different transformations
csv_results <- lapply(csv_files, function(file){
  
  # Get out group info
  Group <- strsplit(file,"_")[[1]][2]
  
  # Read dataframe, mutate Result column different ways
  df <- read.csv(paste0('../fig_3/sim_data/', file), row.names = 1) |> 
    mutate(Result_turnon = Results/min(Results),
           Result_onratio = Results/max(Results))
  
  # Write to csv in wide format for easy plotting
  df_for_plot <- df |> 
    mutate(Group = group_mapping[group],
           Replicate = sample_idx) |> 
    dplyr::select(Group, Replicate, Results, Result_turnon, Result_onratio) |> 
    pivot_wider(names_from = Group, 
                values_from = c(Results, Result_turnon, Result_onratio)) |> 
    dplyr::select(-Replicate)
  
  write.csv(df_for_plot,
            paste0('transformed_sim_data/group_',Group,'_includes_all_transforms.csv'),
            row.names = FALSE)
  
  # Unmodified Result column
  mod1 <- lmer(Results ~ on + (1 | group), data = df)
  Singular <- isSingular(mod1)
  
  if (Singular){ # refit without random effects
    mod1 <- lm(Results ~ on, data = df)
  }
  mod1_test <- glht(mod1, linfct=c('on==0')) |>
    tidy(conf.int = TRUE) 
  
  # Turn on Result column
  mod2 <- lmer(Result_turnon ~ on + (1 | group), data = df)
  Singular <- isSingular(mod2)
  
  if (Singular){ # refit without random effects
    mod2 <- lm(Result_turnon ~ on, data = df)
  }
  mod2_test <- glht(mod2, linfct=c('on==0')) |>
    tidy(conf.int = TRUE) 
  
  # Ratio Result column
  mod3 <- lmer(Result_onratio ~ on + (1 | group), data = df)
  Singular <- isSingular(mod3)
  
  if (Singular){ # refit without random effects
    mod3 <- lm(Result_onratio ~ on, data = df)
  }
  mod3_test <- glht(mod3, linfct=c('on==0')) |>
    tidy(conf.int = TRUE)
  
  df_results <- list(mod1_test, mod2_test, mod3_test) |> 
    bind_rows() |> 
    mutate(Transformation = c('None', 'Turn on', 'On ratio'),
           Group = Group)
}) |> bind_rows() |> 
  dplyr::select(Group, Transformation, estimate, conf.low, conf.high, adj.p.value)

write.csv(csv_results, 'figs10stats.csv', row.names = F)
