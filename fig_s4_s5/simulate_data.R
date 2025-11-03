simulate_data <- function(working.params.df, sample.size, logic, fileNamePrefix, variance = 0.1, write.file = TRUE) {
  # this is for figure s3
  library(lme4)
  
  working.state.names <- row.names(working.params.df)
  
  for (i in 1:length(working.state.names)){
    for (j in 1:length(sample.size)){
      # get out loop stuff
      state <- working.state.names[i]
      n <- sample.size[j]
      
      # simulate data
      var <- variance * (working.params.df |> unlist())
      Results <- lapply(names(var), function(name){
        sim <- rnorm(n, working.params.df[[name]], var[[name]])
      }) |> unlist()
      
      # grouping factors
      group <- rep(c(1:4), each = n) |> as.factor()
      on <- rep(logic, each = n)
      
      # make df
      working.df <- data.frame(Results, on, group)
      
      if (write.file) {
        write.csv(working.df,paste0(fileNamePrefix,'_sample_size_',n,'_simulated_data.csv'))
      }
    }
  }
  
  return(working.df)
}
