simulate_data <- function(working.params.df, sample.size, logic, fileNamePrefix, write.file = TRUE) {
  # This is for figure 3
  library(lme4)
  
  working.state.names <- row.names(working.params.df)

  for (i in 1:length(working.state.names)){
    for (j in 1:length(sample.size)){
      # get out loop stuff
      state <- working.state.names[i]
      n <- sample.size[j]
      
      # simulate data
      var <- 0.1 * (working.params.df |> unlist())
      Results <- lapply(names(var), function(name){
        sim <- rnorm(n, working.params.df[[name]], var[[name]])
      }) |> unlist()

      # grouping factors
      group <- rep(c(1:4), each = n)
      on <- c(rep(0, times = (sum(or.logic == 0) * 3)), 
              rep(1, times = (sum(or.logic == 1) * 3)))
      
      # make df
      working.df <- data.frame(Results, on, group)
      
      if (write.file) {
        write.csv(working.df,paste0(fileNamePrefix,'_sample_size_',n,'_simulated_data.csv'))
      }
    }
  }
  
  return(working.df)
}
