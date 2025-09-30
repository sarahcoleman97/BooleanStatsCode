setwd("~/Library/CloudStorage/Box-Box/BooleanStats/BooleanStatsCode/fig_5")

library(readxl)
library(tidyr)
library(dplyr)
library(broom)
library(lme4)
library(emmeans)
library(multcomp)

reu_gates <- c(240:245)

################################
### Reading adder 3 in 3 out ###
################################

# One of the only ones in REUs (aka comparable)
# There is only one gate in this part of the analysis!

df <- read_excel("scraped_data/ScrapewithReplicates.xlsx", sheet = 'Adder 3 in 3 out') 
colnames(df) <- c("gate", "in1", "in2", "in3", 
                  "Result1", "Result2", "Result3", 
                  "on1", "on2", "on3")
gate <- (df |> pull(gate) |> na.omit())[1]
df <- df |> 
  mutate(group = paste0(in1, in2, in3) |> as.factor())

# No normalization done (because REU)
isREU <- TRUE

# Pivot the df
df_for_multi <- df |> 
  pivot_longer(cols = c(Result1, Result2, Result3), 
               names_prefix = 'Result',
               names_to= 'Output_Type',
               values_to = 'Result') |> 
  pivot_longer(cols = c(on1, on2, on3),
               names_prefix = 'on',
               names_to = 'On_Number',
               values_to = 'on') |> 
  filter(Output_Type == On_Number) |> 
  mutate(Output_Type = Output_Type |> as.factor(),
         group = paste0(group, Output_Type) |> as.factor())
  

# Fit model
lm1 <- lmer(Result ~ on * Output_Type + (1 | group), data = df_for_multi)

# Check for singular fit, refit if so
Singular <- isSingular(lm1)
if (Singular){
  lm1 <- lm(Result ~ on * Output_Type, data = df_for_multi)
}

# Getting specific pairwise comparison of interest (with int)
emm <- emmeans(lm1, ~ Output_Type * on)
contrast_individ_betas <- emm |> contrast(method = "revpairwise", 
                                      by = c("Output_Type"), 
                                      adjust = 'bonferroni') |> 
  tidy(conf.int = TRUE) |> 
  mutate(Gate = gate, 
         Type = 'Adder 3 in 3 out',
         contrast = paste0('on-off within Output_Type',Output_Type)) |> 
  rename(adj.p.value = p.value)

contrast_across_betas <- emm |> contrast(method = list("Output2(on-off)_minus_Output1(on-off)" = c(1, -1, 0, -1, 1, 0),
                                                      "Output3(on-off)_minus_Output1(on-off)" = c(1, 0, -1, -1, 0, 1),
                                                      "Output3(on-off)_minus_Output2(on-off)" = c(0, 1, -1, 0, -1, 1)),
                                        adjust = 'bonferroni') |> 
  tidy(conf.int = TRUE) |> 
  mutate(Gate = gate, 
         Type = 'Adder 3 in 3 out') 

adder3in3out_results <- list(contrast_individ_betas, contrast_across_betas) |> 
  bind_rows() |> 
  mutate(isSingular = Singular,
         isREU = TRUE,
         Gate = as.character(Gate)) |> 
  select(Gate, Type, contrast, isREU, isSingular, estimate, conf.low, conf.high, adj.p.value)

##########################
### Reading half adder ###
##########################

df <- read_excel("scraped_data/ScrapewithReplicates.xlsx", sheet = 'Half adder') 

# Assign sub-dfs (for lapply analysis)
df_list <- list(df[,1:7], df[,9:15], df[,17:23])
half_adder_results <- lapply(df_list, function(df_temp){
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
  
  # Pivot the df
  df_for_multi <- df_norm |> 
    pivot_longer(cols = c(Result1, Result2), 
                 names_prefix = 'Result',
                 names_to= 'Output_Type',
                 values_to = 'Result') |> 
    pivot_longer(cols = c(on1, on2),
                 names_prefix = 'on',
                 names_to = 'On_Number',
                 values_to = 'on') |> 
    filter(Output_Type == On_Number) |> 
    mutate(Output_Type = Output_Type |> as.factor(),
           group = paste0(group, Output_Type) |> as.factor())
  
  # Fit model
  lm1 <- lmer(Result ~ on * Output_Type + (1 | group), data = df_for_multi)
  
  # Check for singular fit, refit if so
  Singular <- isSingular(lm1)
  if (Singular){
    lm1 <- lm(Result ~ on * Output_Type, data = df_for_multi)
  }
  
  # Getting specific pairwise comparison of interest (with int)
  emm <- emmeans(lm1, ~ Output_Type * on)
  contrast_individ_betas <- emm |> contrast(method = "revpairwise", 
                                            by = c("Output_Type"), 
                                            adjust = 'bonferroni') |> 
    tidy(conf.int = TRUE) |> 
    mutate(Gate = gate, 
           Type = 'Half adder',
           contrast = paste0('on-off within Output_Type',Output_Type)) |> 
    rename(adj.p.value = p.value)
  
  contrast_across_betas <- emm |> contrast(method = list("Output2(on-off)_minus_Output1(on-off)" = c(1, -1, -1, 1)),
                                           adjust = 'bonferroni') |> 
    tidy(conf.int = TRUE) |> 
    mutate(Gate = gate, 
           Type = 'Half adder') |> 
    rename(adj.p.value = p.value)
  
  halfadder <- list(contrast_individ_betas, contrast_across_betas) |> 
    bind_rows() |> 
    mutate(isSingular = Singular,
           isREU = isREU,
           Gate = as.character(Gate)) |> 
    select(Gate, Type, contrast, isREU, isSingular, estimate, conf.low, conf.high, adj.p.value)
}) |> bind_rows()


#####################
### Reading adder ###
#####################

df <- read_excel("scraped_data/ScrapewithReplicates.xlsx", sheet = 'Adder') 

# Assign sub-dfs (for lapply analysis)
df_list <- list(df[,1:8], df[,10:17], df[,19:26])
full_adder_results <- lapply(df_list, function(df_temp){
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
  
  # Pivot the df
  df_for_multi <- df_norm |> 
    pivot_longer(cols = c(Result1, Result2), 
                 names_prefix = 'Result',
                 names_to= 'Output_Type',
                 values_to = 'Result') |> 
    pivot_longer(cols = c(on1, on2),
                 names_prefix = 'on',
                 names_to = 'On_Number',
                 values_to = 'on') |> 
    filter(Output_Type == On_Number) |> 
    mutate(Output_Type = Output_Type |> as.factor(),
           group = paste0(group, Output_Type) |> as.factor())
  
  # Fit model
  lm1 <- lmer(Result ~ on * Output_Type + (1 | group), data = df_for_multi)
  
  # Check for singular fit, refit if so
  Singular <- isSingular(lm1)
  if (Singular){
    lm1 <- lm(Result ~ on * Output_Type, data = df_for_multi)
  }
  
  # Getting specific pairwise comparison of interest (with int)
  emm <- emmeans(lm1, ~ Output_Type * on)
  contrast_individ_betas <- emm |> contrast(method = "revpairwise", 
                                            by = c("Output_Type"), 
                                            adjust = 'bonferroni') |> 
    tidy(conf.int = TRUE) |> 
    mutate(Gate = gate, 
           Type = 'Full adder',
           contrast = paste0('on-off within Output_Type',Output_Type)) |> 
    rename(adj.p.value = p.value)
  
  contrast_across_betas <- emm |> contrast(method = list("Output2(on-off)_minus_Output1(on-off)" = c(1, -1, -1, 1)),
                                           adjust = 'bonferroni') |> 
    tidy(conf.int = TRUE) |> 
    mutate(Gate = gate, 
           Type = 'Full adder') |> 
    rename(adj.p.value = p.value)
  
  fulladder <- list(contrast_individ_betas, contrast_across_betas) |> 
    bind_rows() |> 
    mutate(isSingular = Singular,
           isREU = isREU,
           Gate = as.character(Gate)) |> 
    select(Gate, Type, contrast, isREU, isSingular, estimate, conf.low, conf.high, adj.p.value)
}) |> bind_rows()

#########################
### Reading Voigt sc7 ###
#########################

# Define logic
logic_yfp <- list('-/-/-'=0,'+/-/-'=1,'-/+/-'=0,'-/-/+'=0,
                  '+/+/-'=1,'+/-/+'=0,'-/+/+'=0,'+/+/+'=0) 
logic_rfp <- list('-/-/-'=0,'+/-/-'=0,'-/+/-'=0,'-/-/+'=0,
                  '+/+/-'=1,'+/-/+'=0,'-/+/+'=0,'+/+/+'=1)

# This is REU
isREU <- TRUE

# Make dataframes
logic_df_yfp <- data.frame(
  group = as.factor(names(logic_yfp)),
  on_YFP = unlist(logic_yfp),
  row.names = NULL
) 
logic_df_rfp <- data.frame(
  group = as.factor(names(logic_rfp)),
  on_RFP = unlist(logic_rfp),
  row.names = NULL
) 

# Read dataframe
voigt_path <- "scraped_data/data_s41589-024-01730-1.xlsx"
gate <- 'sc7'
df <- (read_excel(voigt_path, sheet = gate, skip = 1))[,c(1, 2, 4)] 
colnames(df) <- c('group', 'YFP', 'RFP')
df <- df |> 
  fill(group, .direction = "down")

df <- df |> 
  right_join(logic_df_yfp, by = 'group', relationship = 'many-to-one') |> 
  right_join(logic_df_rfp, by = 'group', relationship = 'many-to-one')

df_for_multi <- df |> # no normalization
  pivot_longer(cols = c(YFP, RFP), 
               names_to= 'Output_Type',
               values_to = 'Result') |> 
  pivot_longer(cols = c(on_YFP, on_RFP),
               names_prefix = 'on_',
               names_to = 'On_Number',
               values_to = 'on') |> 
  filter(Output_Type == On_Number) |> 
  mutate(Output_Type = Output_Type |> as.factor(),
         group = paste0(group, Output_Type) |> as.factor())

# Fit model
lm1 <- lmer(Result ~ on * Output_Type + (1 | group), data = df_for_multi)

# Check for singular fit, refit if so
Singular <- isSingular(lm1)
if (Singular){
  lm1 <- lm(Result ~ on * Output_Type, data = df_for_multi)
}

# Getting specific pairwise comparison of interest (with int)
emm <- emmeans(lm1, ~ Output_Type * on)
contrast_individ_betas <- emm |> contrast(method = "revpairwise", 
                                          by = c("Output_Type"), 
                                          adjust = 'bonferroni') |> 
  tidy(conf.int = TRUE) |> 
  mutate(Gate = gate, 
         Type = 'Full adder',
         contrast = paste0('on-off within Output_Type',Output_Type)) |> 
  rename(adj.p.value = p.value)

contrast_across_betas <- emm |> contrast(method = list("OutputYFP(on-off)_minus_OutputRFP(on-off)" = c(1, -1, -1, 1)),
                                         adjust = 'bonferroni') |> 
  tidy(conf.int = TRUE) |> 
  mutate(Gate = gate, 
         Type = 'Full adder') |> 
  rename(adj.p.value = p.value)

voigt <- list(contrast_individ_betas, contrast_across_betas) |> 
  bind_rows() |> 
  mutate(isSingular = Singular,
         isREU = isREU,
         Gate = as.character(Gate)) |> 
  select(Gate, Type, contrast, isREU, isSingular, estimate, conf.low, conf.high, adj.p.value)

# Writing and saving everything
out_data <- list(adder3in3out_results, half_adder_results, 
                 full_adder_results, voigt) |> bind_rows() |> 
  arrange(Gate)

write.csv(out_data, 'fig5_scraped_multiinput_stats.csv', row.names = FALSE)
