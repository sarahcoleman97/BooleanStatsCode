setwd("~/Library/CloudStorage/Box-Box/BooleanStats/BooleanStatsCode/fig_5")

library(lme4)
library(dplyr)
library(multcomp)
library(readxl)
library(tidyr)
library(broom)

### NEED TO UPDATE THE BELOW WITH CORRECT DATA ###

# Define column names
colNames <- c("off", "on_1", "on_2", "on_3")

# First do the CSV files
csv_files <- list.files(path = 'data/', pattern = '*.csv') 

csv_results <- lapply(csv_files, function(file){
  filePath <- paste0('data/', file)
  fileName <- strsplit(file, '*.csv')[[1]]
  df <- read.csv(filePath)  |>  
    rename_with(~colNames) |> 
    pivot_longer(cols = all_of(colNames), 
                              names_to = 'group', 
                              values_to = 'value') |>  
    mutate(expected = case_when(grepl('on', group) ~ 1,
                                grepl('off', group) ~ 0,
                                TRUE ~ NA),
           group = as.factor(group))
  lm <- lmer(value ~ expected + (1|group), data = df)
  Singular <- isSingular(lm)
  if (Singular){ # refit without random effects
    lm <- lm(Results ~ on, data = working.df)
  }
  working.test <- glht(lm,linfct=c('expected==0')) |>
    tidy(conf.int = TRUE) |> 
    mutate(Gate = fileName,
           isSingular = Singular) |> 
    select(Gate, isSingular, estimate, conf.low, conf.high, adj.p.value)
}) |> bind_rows()

# Now the XLSX files
# each is formatted differently, treat sep
xl_df <- read_xlsx('data/GMP5-022.xlsx', sheet = 'REU', skip = 3)
df_list <- list(xl_df[,3:6], xl_df[,8:11], xl_df[,13:16], 
             xl_df[,18:21], xl_df[,23:26], xl_df[,28:31])
gate_names <- c("Trad", "Med", "Med.Low", "Low", "LuxR", "Cons")

xlsx_results_1 <- lapply(1:length(df_list), function(i){
  df <- df_list[[i]]
  gate <- gate_names[i]
  working.df <- df |> 
    rename_with(~colNames) |> 
    pivot_longer(cols = all_of(colNames), 
                 names_to = 'group', 
                 values_to = 'value') |>  
    mutate(expected = case_when(grepl('on', group) ~ 1,
                                grepl('off', group) ~ 0,
                                TRUE ~ NA),
           group = as.factor(group))
  lm <- lmer(value ~ expected + (1|group), data = working.df)
  Singular <- isSingular(lm)
  if (Singular){ # refit without random effects
    lm <- lm(value ~ expected, data = working.df)
  }
  working.test <- glht(lm,linfct=c('expected==0')) |>
    tidy(conf.int = TRUE) |> 
    mutate(isSingular = Singular,
           Gate = gate) |> 
    select(Gate, isSingular, estimate, conf.low, conf.high, adj.p.value)
}) |> bind_rows()

# The second XLSX file
df <- read_xlsx('data/New.Lux.1000.xlsx', sheet = 'Sheet1')
gate <- 'New.Lux.1000'
working.df <- df |> 
  rename_with(~colNames) |> 
  pivot_longer(cols = all_of(colNames), 
               names_to = 'group', 
               values_to = 'value') |>  
  mutate(expected = case_when(grepl('on', group) ~ 1,
                              grepl('off', group) ~ 0,
                              TRUE ~ NA),
         group = as.factor(group))
lm <- lmer(value ~ expected + (1|group), data = working.df)
Singular <- isSingular(lm)
if (Singular){ # refit without random effects
  lm <- lm(value ~ expected, data = working.df)
}
xlsx_results_2 <- glht(lm,linfct=c('expected==0')) |>
  tidy(conf.int = TRUE) |> 
  mutate(isSingular = Singular,
         Gate = gate) |> 
  select(Gate, isSingular, estimate, conf.low, conf.high, adj.p.value)

# Combine results into one dataframe
all_results <- list(csv_results, xlsx_results_1, xlsx_results_2) |> 
  bind_rows()

# Write results to file
write.csv(all_results, 'fig5stats.csv')
