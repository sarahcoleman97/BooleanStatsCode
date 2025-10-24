setwd("~/Library/CloudStorage/Box-Box/BooleanStats/BooleanStatsCode/fig_5")

library(lme4)
library(dplyr)
library(multcomp)
library(readxl)
library(tidyr)
library(broom)

# Reading the .xslx files
files <- list.files(path = 'data/', pattern = '*.xlsx')
gate_names <- strsplit(files, '.xlsx')

df_list <- lapply(files, function(file){
  if (file == "AlexisnewGate.xlsx"){
    df <- read_xlsx(paste0('data/',file), skip = 1) # remove on/off indiciation in top column
  } else {
    df <- read_xlsx(paste0('data/',file))
  }
})

# order is aTc, OC6, IPTG
colNameToLogics <- list('-/-/-' = '000',
                        'aTc' = '100',
                        'OC6' = '010',
                        'aTc/OC6' = '110',
                        'aTc/IPTG' = '101',
                        'OC6/IPTG' = '011',
                        'IPTG' = '001',
                        'aTc/OC6/IPTG' = '111')

conditionLogics <- list('000' = 1, #-/-/-
                       '100' = 0, # aTc
                       '010' = 0, # OC6
                       '110' = 0, # aTc/OC6
                       '101' = 1, # aTc/IPTG
                       '011' = 1, # OC6/IPTG
                       '001' = 1, # IPTG
                       '111' = 1) # aTc/OCT/IPTG

xlsx_results <- lapply(1:length(df_list), function(i){
  df <- df_list[[i]]
  gate <- gate_names[[i]]
  working.df <- df |> 
    rename_with(~ unlist(colNameToLogics[.x])) |> 
    pivot_longer(cols = all_of(unname(unlist(colNameToLogics))), 
                 names_to = 'group', 
                 values_to = 'value') |>  
    mutate(expected = unlist(conditionLogics)[group],
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
    dplyr::select(Gate, isSingular, estimate, conf.low, conf.high, adj.p.value)
}) |> bind_rows()

# Write results to file
write.csv(xlsx_results, 'fig5stats.csv', row.names = FALSE)
