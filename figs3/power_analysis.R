setwd("~/Library/CloudStorage/Box-Box/BooleanStats/BooleanStatsCode/figs3")

source("simulate_data.R") # allows to pass a variance input, fraction of the mean response.

set.seed(2024) 

library('lme4')
library('effects')
library('broom')
library('multcomp')

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
new.cond <- read.csv("simulation_conditions.csv",row.names=1,header=TRUE)
conditions <- new.cond |> row.names()
conditions_logics <- c(conditions[1:8], "OR") # specific to this case
samples <- 2:30
n.mc <- 10 # number of monte carlo simulations
variance.values <- 1:10 * 0.1
or.logic <- c(0, 1, 1, 1)

conditions_idx <- 1:length(conditions)
variance_idx <- 1:length(variance.values)
samples_idx <- 1:length(samples)
mc_idx <- 1:n.mc

result_top <- parallel::mclapply(conditions_idx, function(i){
  condition <- conditions[i]
  logic <- logic_list[[conditions_logics[[i]]]]
  
  result_n1 <- parallel::mclapply(variance_idx, function(m){
    variance <- variance.values[m]
    beta.store <- matrix(nrow=length(samples),ncol=n.mc) # init
    p.value <- matrix(nrow=length(samples),ncol=n.mc) # init
    singular.store <- matrix(nrow=length(samples),ncol=n.mc) # init

    result_n2 <- parallel::mclapply(samples_idx, function(j){ 
      n <- samples[j]
      
      result_n3 <- parallel::mclapply(mc_idx, function(k){ 
        working.df <- simulate_data(new.cond[i,], 
                                 n, 
                                 logic, 
                                 paste0('sim_data/',
                                        condition),
                                 variance = variance,
                                 write.file = FALSE)
        working.lmer <- lmer(Results ~ on + (1|group), data = working.df)
        Singular <- isSingular(working.lmer)
    
        if (Singular){ # refit without random effects
          working.lmer <- lm(Results ~ on, data = working.df)
          beta.estimate <- coef(working.lmer)[['on']]
        } else {
          beta.estimate <- coef(working.lmer)[['group']][['on']][1]
        }
        
        answer <- summary(glht(working.lmer,linfct=c('on==0')))
        beta.store[j,k] <- beta.estimate
        p.value[j,k] <- answer$test$pvalues[1]
        singular.store[j,k] <- Singular
      })
    })
    
    # Write outputs
    beta.store <- data.frame(beta.store, row.names = samples)
    colnames(beta.store) <- mc_idx
    beta.store <- t(beta.store)
    write.csv(beta.store,paste0('beta_montecarlo/',condition,
                                '_betas_var_',variance,'.csv'))
      
    p.value <- data.frame(p.value, row.names = samples)
    colnames(p.value) <- mc_idx
    p.value <- t(p.value)
    write.csv(p.value,paste0('beta_montecarlo/',condition,
                             '_pvals_var_',variance,'.csv'))
    
    singularities <- data.frame(singular.store, row.names = samples)
    colnames(singularities) <- mc_idx
    singularities <- t(singularities)
    write.csv(singularities,paste0('beta_montecarlo/',condition,
                                   '_singularities_var_',variance,'.csv'))
  })
})
