#function_count_formula_vars.R
# Function to count the number of variables in a model formula
extract_formula_vars <- function(formula) {
  # Extract term labels
  term_labels <- attr(terms(formula), "term.labels")
  # Get the number of variables
  num_vars <- length(term_labels)
  # Return as a list
  return(list(variable_names = term_labels, number_of_variables = num_vars))
}
