setwd("~/Library/CloudStorage/Box-Box/BooleanStats/BooleanStatsCode/fig_5")

library(readxl)
library(dplyr)
library(tidyr)
library(broom)
library(lme4)
library(multcomp)

# A note on normalization: some gates are REU (240-245) and the rest are not.
# I make a condition here to normalize them correctly

reu_gates <- c(240:245)

# A second note on interpretation: for all these gates, the question that is 
# being answered is whether the effect of on (beta) is statistically different 
# from zero. The CIs CANNOT be used to compare whether the effect of on (beta)
# is statistically different from another effect. 

#################################
### Reading three input gates ###
#################################

###
# PART 1: Gina's scrape
###

df <- read_excel("scraped_data/ScrapewithReplicates.xlsx", sheet = '3 input') 

# Assign sub-dfs (for lapply analysis)
df_list <- list(df[,1:6], 
                df[,8:13])

# Define column names
colNames <- c("gate", "in1", "in2", "in3", "Result", "on")

# Analysis loop
results_3in <- lapply(df_list, function(df_temp){

  df_temp <- df_temp |> 
    rename_with(~colNames)
  gate <- (df_temp |> pull(gate) |> na.omit())[1]
  df_temp <- df_temp |> 
    fill(in1, in2, in3, on, .direction = "down") |> 
    mutate(group = paste0(in1, in2, in3) |> as.factor())
  
  # Normalize by ON (only for NON REU gates)
  if (gate %in% reu_gates){
    df_norm <- df_temp # no normalization done
    isREU <- TRUE
  } else {
    avg_on <- df_temp |> filter(on == 1) |> pull(Result) |> mean(na.rm = TRUE)
    df_norm <- df_temp |> mutate(Result = Result/avg_on)
    isREU <- FALSE
  }
  
  # Fit model
  lm <- lmer(Result ~ on + (1 | group), data = df_norm)
  
  # Results
  Singular <- isSingular(lm)
  
  if (Singular){ # refit without random effects
    lm <- lm(Result ~ on, data = df_norm)
  }
  
  working.test <- glht(lm,linfct=c('on==0')) |>
    tidy(conf.int = TRUE) |> 
    mutate(Gate = gate,
           isSingular = Singular,
           isREU = isREU)  
}) |> bind_rows() |> 
  mutate(Type = '3 input',
         Gate = as.character(Gate)) |> 
  select(Gate, Type, isREU, isSingular, estimate, conf.low, conf.high, adj.p.value)

###
# PART 2: Voigt SI (necessary bc different data structure)
### 

voigt_path <- "scraped_data/data_s41589-024-01730-1.xlsx"
voigt_gates <- c('sc3', 'sc7', 'sc8', 'sc35', 'sc38', 'sc39', 'sc40')
voigt_df_list <- lapply(voigt_gates, function(voigt_gate){
  df <- read_excel(voigt_path, sheet = voigt_gate, skip = 1)
  df <- df |> 
    mutate(Gate = voigt_gate)
})
names(voigt_df_list) <- voigt_gates

# YFP/RFP gate logic
gate_logic <- list('sc3' = list(
                     'YFP' = list('-/-/-'=0,'+/-/-'=1,'-/+/-'=1,'-/-/+'=0,
                                '+/+/-'=1,'+/-/+'=0,'-/+/+'=0,'+/+/+'=0)),
                   'sc7' = list( # YFP/RFP
                     'YFP' = list('-/-/-'=0,'+/-/-'=1,'-/+/-'=0,'-/-/+'=0,
                                  '+/+/-'=1,'+/-/+'=0,'-/+/+'=0,'+/+/+'=0),
                     'RFP' = list('-/-/-'=0,'+/-/-'=0,'-/+/-'=0,'-/-/+'=0,
                                  '+/+/-'=1,'+/-/+'=0,'-/+/+'=0,'+/+/+'=1)), 
                   'sc8' = list(
                     'YFP' = list('-/-/-'=0,'+/-/-'=1,'-/+/-'=1,'-/-/+'=0,
                                '+/+/-'=1,'+/-/+'=0,'-/+/+'=0,'+/+/+'=0)),
                   'sc35' = list(
                     'YFP' = list('-/-/-'=0,'+/-/-'=0,'-/+/-'=1,'-/-/+'=0,
                                 '+/+/-'=0,'+/-/+'=0,'-/+/+'=0,'+/+/+'=0)),
                   'sc38' = list(
                     'YFP' = list('-/-/-'=0,'+/-/-'=1,'-/+/-'=1,'-/-/+'=0,
                                 '+/+/-'=1,'+/-/+'=0,'-/+/+'=0,'+/+/+'=0)),
                   'sc39' = list(
                     'YFP' = list('-/-/-'=0,'+/-/-'=0,'-/+/-'=1,'-/-/+'=0,
                                 '+/+/-'=0,'+/-/+'=1,'-/+/+'=1,'+/+/+'=1)),
                   'sc40' = list(
                     'YFP' = list('-/-/-'=0,'+/-/-'=0,'-/+/-'=1,'-/-/+'=0,
                                 '+/+/-'=0,'+/-/+'=1,'-/+/+'=1,'+/+/+'=1)))
###########
### YFP ###
###########

yfp_cols <- c(1,2) # Note these are the same in every sheet

results_3in_y <- lapply(voigt_df_list, function(voigt_df){
  
  # Get the gate
  gate <- voigt_df |> pull(Gate) |> unique()
  
  # Get only the YFP columns and format correctly
  df <- voigt_df[,c(1,2)]
  colnames(df) = c('group', 'Result')
  df <- df |> 
    fill(group, .direction = "down") 
  
  # Make the logic
  logic <- gate_logic[[gate]][['YFP']]
  logic_df <- data.frame(
    group = as.factor(names(logic)),
    on = unlist(logic),
    row.names = NULL
  ) 
  
  # Merge to get final df
  df_norm <- df %>% 
    right_join(logic_df, by = 'group', relationship = 'many-to-one') |> 
    mutate(group = group |> as.factor())
  
  # Fit model
  lm <- lmer(Result ~ on + (1 | group), data = df_norm)
  
  # Results
  Singular <- isSingular(lm)
  
  if (Singular){ # refit without random effects
    lm <- lm(Result ~ on, data = df_norm)
  }
  
  test <- glht(lm,linfct=c('on==0')) |>
    tidy(conf.int = TRUE) |> 
    mutate(Gate = gate,
           isSingular = Singular)  
}) |> bind_rows() |> 
  mutate(Type = '3 input Voigt YFP',
         isREU = TRUE,
         Gate = as.character(Gate)) |> 
  select(Gate, Type, isREU, isSingular, estimate, conf.low, conf.high, adj.p.value)

###########
### RFP ###
###########

gate <- 'sc7'
cols <- c(1,4)
df <- (voigt_df_list[[gate]])[,cols]
colnames(df) = c('group', 'Result')
df <- df |> 
  fill(group, .direction = "down") 

# Make the logic
logic <- gate_logic[[gate]][['RFP']]
logic_df <- data.frame(
  group = as.factor(names(logic)),
  on = unlist(logic),
  row.names = NULL
) 

# Merge to get final df
df_norm <- df %>% 
  right_join(logic_df, by = 'group', relationship = 'many-to-one') |> 
  mutate(group = group |> as.factor())

# Fit model
lm <- lmer(Result ~ on + (1 | group), data = df_norm)
glht <- lm |> glht(linfct=c('on==0'))
confint <- glht |> confint()

# Results
Singular <- isSingular(lm)

if (Singular){ # refit without random effects
  lm <- lm(Result ~ on, data = df_norm)
}

results_3in_r <- glht(lm,linfct=c('on==0')) |>
  tidy(conf.int = TRUE) |> 
  mutate(Gate = gate,
         isSingular = Singular) |> 
  mutate(Type = '3 input Voigt RFP',
         isREU = TRUE,
         Gate = as.character(Gate)) |> 
  select(Gate, Type, isREU, isSingular, estimate, conf.low, conf.high, adj.p.value)

################################
### Reading four input gates ###
################################

df <- read_excel("scraped_data/ScrapewithReplicates.xlsx", sheet = '4 input') 

# Assign sub-dfs (for lapply analysis)
df_list <- list(df[,1:7], df[,8:14], df[,15:21],df[,23:29], df[,31:37]) # df[,39:45] not used
results_4in <- lapply(df_list, function(df_temp){
  colnames(df_temp) <- c("gate", "in1", "in2", "in3", "in4", "Result", "on")
  gate <- (df_temp |> pull(gate) |> na.omit())[1]
  df_temp <- df_temp |> 
    fill(in1, in2, in3, in4, on, .direction = "down") |> 
    mutate(group = paste0(in1, in2, in3, in4) |> as.factor()) |> 
    drop_na(Result) # if there are any NA in result it should be dropped
  
  # Normalize by ON (only for NON REU gates)
  if (gate %in% reu_gates){
    df_norm <- df_temp # no normalization done
    isREU <- TRUE
  } else {
    avg_on <- df_temp |> filter(on == 1) |> pull(Result) |> mean(na.rm = TRUE)
    df_norm <- df_temp |> mutate(Result = Result/avg_on)
    isREU <- FALSE
  }
  
  # Fit model
  lm <- lmer(Result ~ on + (1 | group), data = df_norm)
  
  # Results
  Singular <- isSingular(lm)
  
  if (Singular){ # refit without random effects
    lm <- lm(Result ~ on, data = df_norm)
  }
  
  working.test <- glht(lm,linfct=c('on==0')) |>
    tidy(conf.int = TRUE) |> 
    mutate(Gate = gate,
           isSingular = Singular,
           isREU = isREU)  
}) |> bind_rows() |> 
  mutate(Type = '4 input',
         Gate = as.character(Gate)) |> 
  select(Gate, Type, isREU, isSingular, estimate, conf.low, conf.high, adj.p.value)

##############################
### Reading six input gate ###
##############################

df <- read_excel("scraped_data/ScrapewithReplicates.xlsx", sheet = '6 input') 
colnames(df) <- c("gate", "in1", "in2", "in3", "in4", "in5", "in6", "Result", "on")
df <- df |> 
  mutate(group = paste0(in1, in2, in3, in4, in5, in6) |> as.factor())
gate <- (df |> pull(gate) |> na.omit())[1]

# Normalize by ON (only for NON REU gates)
if (gate %in% reu_gates){
  df_norm <- df # no normalization done
  isREU <- TRUE
} else {
  avg_on <- df |> filter(on == 1) |> pull(Result) |> mean(na.rm = TRUE)
  df_norm <- df |> mutate(Result = Result/avg_on)
  isREU <- FALSE
}

# Fit model
lm <- lmer(Result ~ on + (1 | group), data = df_norm)

# Results
Singular <- isSingular(lm)

if (Singular){ # refit without random effects
  lm <- lm(Result ~ on, data = df_norm)
}

results_6in <- glht(lm,linfct=c('on==0')) |>
  tidy(conf.int = TRUE) |> 
  mutate(Gate = gate,
         Type = '6 input',
         isSingular = Singular,
         isREU = isREU,
         Gate = as.character(Gate)) |> 
  select(Gate, Type, isREU, isSingular, estimate, conf.low, conf.high, adj.p.value)

################################
### Reading adder 3 in 3 out ###
################################

df <- read_excel("scraped_data/ScrapewithReplicates.xlsx", sheet = 'Adder 3 in 3 out') 
colnames(df) <- c("gate", "in1", "in2", "in3", "Result1", "Result2", "Result3", "on1", "on2", "on3")
gate <- (df |> pull(gate) |> na.omit())[1]
df <- df |> 
  mutate(group = paste0(in1, in2, in3) |> as.factor())

# Normalize by ON (only for NON REU gates)
if (gate %in% reu_gates){
  df_norm <- df # no normalization done
  isREU <- TRUE
} else {
  on_1 <- df |> filter(on1 == 1) |> pull(Result1) |> mean(na.rm = TRUE)
  on_2 <- df |> filter(on2 == 1) |> pull(Result2) |> mean(na.rm = TRUE)
  on_3 <- df |> filter(on3 == 1) |> pull(Result3) |> mean(na.rm = TRUE)
  df_norm <- df |> 
    mutate(Result1 = Result1/on_1,
           Result2 = Result2/on_2,
           Result3 = Result3/on_3)
  isREU <- FALSE
}

# Fitting model 1
lm1 <- lmer(Result1 ~ on1 + (1 | group), data = df_norm)
Singular1 <- isSingular(lm1)

if (Singular1){ # refit without random effects
  lm1 <- lm(Result1 ~ on1, data = df_norm)
}
add1 <- glht(lm1,linfct=c('on1==0')) |>
  tidy(conf.int = TRUE) |> 
  mutate(Gate = gate,
         Type = 'Adder 3 in 3 out, output 1',
         isSingular = Singular1,
         isREU = isREU)

# Fitting model 2
lm2 <- lmer(Result2 ~ on2 + (1 | group), data = df_norm)
Singular2 <- isSingular(lm2)

if (Singular2){ # refit without random effects
  lm2 <- lm(Result2 ~ on2, data = df_norm)
}
add2 <- glht(lm2,linfct=c('on2==0')) |>
  tidy(conf.int = TRUE) |> 
  mutate(Gate = gate,
         Type = 'Adder 3 in 3 out, output 2',
         isSingular = Singular2,
         isREU = isREU)

# Fitting model 3
lm3 <- lmer(Result3 ~ on3 + (1 | group), data = df_norm)
Singular3 <- isSingular(lm3)

if (Singular3){ # refit without random effects
  lm3 <- lm(Result3 ~ on3, data = df_norm)
}
add3 <- glht(lm3,linfct=c('on3==0')) |>
  tidy(conf.int = TRUE) |> 
  mutate(Gate = gate,
         Type = 'Adder 3 in 3 out, output 3',
         isSingular = Singular3,
         isREU = isREU)

results_3in3out <- list(add1, add2, add3) |> 
  bind_rows() |> 
  mutate(Gate = as.character(Gate)) |> 
  select(Gate, Type, isREU, isSingular, estimate, conf.low, conf.high, adj.p.value)

##########################
### Reading half adder ###
##########################

df <- read_excel("scraped_data/ScrapewithReplicates.xlsx", sheet = 'Half adder') 

# Assign sub-dfs (for lapply analysis)
df_list <- list(df[,1:7], df[,9:15], df[,17:23])
results_half_add <- lapply(df_list, function(df_temp){
  colnames(df_temp) <- c("gate", "in1", "in2", "Result1", "Result2", "on1", "on2")
  gate <- (df_temp |> pull(gate) |> na.omit())[1]
  df_temp <- df_temp |> 
    fill(in1, in2, on1, on2, .direction = "down") |> 
    mutate(group = paste0(in1, in2) |> as.factor())
  
  # Normalize by ON (only for NON REU gates)
  if (gate %in% reu_gates){
    df_norm <- df_temp # no normalization done
    isREU <- TRUE
  } else {
    on_1 <- df_temp |> filter(on1 == 1) |> pull(Result1) |> mean(na.rm = TRUE)
    on_2 <- df_temp |> filter(on2 == 1) |> pull(Result2) |> mean(na.rm = TRUE)
    df_norm <- df_temp |> 
      mutate(Result1 = Result1/on_1,
             Result2 = Result2/on_2)
    isREU <- FALSE
  }
  
  # Fitting model 1
  lm1 <- lmer(Result1 ~ on1 + (1 | group), data = df_norm)
  Singular1 <- isSingular(lm1)
  
  if (Singular1){ # refit without random effects
    lm1 <- lm(Result1 ~ on1, data = df_norm)
  }
  half_add1 <- glht(lm1,linfct=c('on1==0')) |>
    tidy(conf.int = TRUE) |> 
    mutate(Gate = gate,
           Type = 'Half adder, output 1',
           isSingular = Singular1,
           isREU = isREU)
  
  # Fitting model 2
  lm2 <- lmer(Result2 ~ on2 + (1 | group), data = df_norm)
  Singular2 <- isSingular(lm2)
  
  if (Singular2){ # refit without random effects
    lm2 <- lm(Result2 ~ on2, data = df_norm)
  }
  half_add2 <- glht(lm2,linfct=c('on2==0')) |>
    tidy(conf.int = TRUE) |> 
    mutate(Gate = gate,
           Type = 'Half adder, output 2',
           isSingular = Singular2,
           isREU = isREU)
  
  result <- list(half_add1, half_add2) |> 
    bind_rows()
}) |> bind_rows() |> 
  mutate(Gate = as.character(Gate)) |> 
  select(Gate, Type, isREU, isSingular, estimate, conf.low, conf.high, adj.p.value)
  

#####################
### Reading adder ###
#####################

df <- read_excel("scraped_data/ScrapewithReplicates.xlsx", sheet = 'Adder') 

# Assign sub-dfs (for lapply analysis)
df_list <- list(df[,1:8], df[,10:17], df[,19:26])
results_add <- lapply(df_list, function(df_temp){
  colnames(df_temp) <- c("gate", "in1", "in2", "in3", "Result1", "Result2", "on1", "on2")
  gate <- (df_temp |> pull(gate) |> na.omit())[1]
  df_temp <- df_temp |> 
    fill(in1, in2, in3, on1, on2, .direction = "down") |> 
    mutate(group = paste0(in1, in2, in3) |> as.factor())
  
  # Normalize by ON (only for NON REU gates)
  if (gate %in% reu_gates){
    df_norm <- df_temp # no normalization done
    isREU <- TRUE
  } else {
    on_1 <- df_temp |> filter(on1 == 1) |> pull(Result1) |> mean(na.rm = TRUE)
    on_2 <- df_temp |> filter(on2 == 1) |> pull(Result2) |> mean(na.rm = TRUE)
    df_norm <- df_temp |> 
      mutate(Result1 = Result1/on_1,
             Result2 = Result2/on_2)
    isREU <- FALSE
  }
  
  # Fitting model 1
  lm1 <- lmer(Result1 ~ on1 + (1 | group), data = df_norm)
  Singular1 <- isSingular(lm1)
  
  if (Singular1){ # refit without random effects
    lm1 <- lm(Result1 ~ on1, data = df_norm)
  }
  full_add1 <- glht(lm1,linfct=c('on1==0')) |>
    tidy(conf.int = TRUE) |> 
    mutate(Gate = gate,
           Type = 'Full adder, output 1',
           isSingular = Singular1,
           isREU = isREU)
  
  # Fitting model 2
  lm2 <- lmer(Result2 ~ on2 + (1 | group), data = df_norm)
  Singular2 <- isSingular(lm2)
  
  if (Singular2){ # refit without random effects
    lm2 <- lm(Result2 ~ on2, data = df_norm)
  }
  full_add2 <- glht(lm2,linfct=c('on2==0')) |>
    tidy(conf.int = TRUE) |> 
    mutate(Gate = gate,
           Type = 'Full adder, output 2',
           isSingular = Singular2,
           isREU = isREU)
  
  result <- list(full_add1, full_add2) |> 
    bind_rows()
}) |> bind_rows() |> 
  mutate(Gate = as.character(Gate)) |> 
  select(Gate, Type, isREU, isSingular, estimate, conf.low, conf.high, adj.p.value)

# Writing a single csv file for the whole thing

out_file <- list(results_3in, results_3in_y, results_3in_r, results_3in3out,
                 results_4in, results_6in, results_add, results_half_add) |> 
  bind_rows() |> 
  arrange(Gate)

write.csv(out_file, 'fig5_scraped_stats.csv')
