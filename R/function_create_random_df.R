# function_create_random_df.R
#Generate df of random values for model predictions
create_random_df <- function(n_rows, n_cols, dist_func = rnorm, ...) {
  # Create a list of random columns
  random_columns <- lapply(1:n_cols, function(x) dist_func(n_rows, ...))
  # Convert the list to a data frame and set column names
  random_df <- as.data.frame(random_columns)
  colnames(random_df) <- paste0("V", 1:n_cols)
  return(random_df)
}
