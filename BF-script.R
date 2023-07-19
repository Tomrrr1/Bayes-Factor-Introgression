#!/usr/bin/env Rscript

# Load necessary packages

library(stats) 


# Parse the command line arguments 

args <- commandArgs(trailingOnly = TRUE)

if (any(commandArgs(trailingOnly = TRUE) %in% "--help")) {

    cat("
\033[32mThis script calculates Bayes Factors using either the BF_Gamma or BF_Beta function.
  
Usage:
 	Rscript BF-script.R [function] [alpha] [beta] [epsilon] [mcmc_file_path] [column_indices]
  
Options:
        	function        Either BF_Gamma or BF_Beta.
        	alpha           Numeric value for the alpha parameter.
        	beta            Numeric value for the beta parameter.
        	epsilon         Numeric value for the epsilon parameter.
        	file            Path to the MCMC sample file.
        	column_indices  Column indices, separated by spaces. Ranges can be specified using ':'.

	If the BF_Gamma option is selected then 'alpha' is the shape and 'beta' is the rate of the prior distribution.
	If the BF_Beta option is selected then both 'alpha' and 'beta' are shape parameters of the prior distribution.
  
Example:
  	Rscript BF-script.R BF_Gamma 2 10 0.1 /path/to/mcmc_file.txt 10 11:13\033[39m
    ")
    
    quit(save="no", status=0)
}


function_name <- args[1]
alpha <- as.numeric(args[2])
beta <- as.numeric(args[3])
epsilon <- as.numeric(args[4])
file <- args[5]



# Handle column indices

column_args <- args[6:length(args)]

column_index <- integer(0)  # Empty integer vector to hold indices

for (arg in column_args) {

  if (grepl(":", arg)) {

    range <- as.integer(strsplit(arg, ":")[[1]])

    column_index <- c(column_index, seq(range[1], range[2]))

  } else {

    column_index <- c(column_index, as.integer(arg))

  }

}


# Read in the MCMC sample file

mcmc_file <- read.table(file, sep="\t", header=TRUE)


# Define a helper function that computes the Bayes Factor for a single column index

compute_result_gamma <- function(i, mcmc_file, epsilon, alpha, beta) {

  mig_rate = mcmc_file[, i]
  
  prior_prop = pgamma(epsilon, shape=alpha, rate=beta) # proportion of prior < epsilon (gamma prior)
  
  posterior_prop = mean(mig_rate < epsilon) # proportion of posterior mcmc samples < epsilon
  
  result = prior_prop / posterior_prop
  
  return(result)
}



BF_Gamma <- function(alpha, beta, epsilon, mcmc_file, column_index) {
  
  # Check for valid alpha, beta, epsilon
  if(alpha <= 0 || beta <= 0 || epsilon <= 0) {
    
    stop("The parameters 'alpha', 'beta' and 'epsilon' must be greater than zero.")
    
  }

  # Check if any indices are out-of-bounds
  if (!all(column_index %in% 1:ncol(mcmc_file))) {
    
    stop("One or more column indices are out of bounds.")
    
  }
  
  # Compute result for each column index and store in a data frame
  results_df = data.frame(Bayes_Factor = sapply(column_index, compute_result_gamma, mcmc_file, epsilon, alpha, beta))
  
  # Generate the row names
  row_names <- paste0(colnames(mcmc_file[column_index]), " (", as.character(epsilon), ")")
  
  # Change the ".." to "->". This is specific to the BPP convention of naming columns
  row_names <- sapply(row_names, function(row) gsub("\\.\\.", " -> ", row))
  
  # Assign the row names to the results_df
  rownames(results_df) <- row_names


  # Define output file path
  output_file_path <- paste0("BF_Gamma_results_", Sys.Date(), "_", ".txt")

  # Write the results to a text file
  write.table(results_df, file = output_file_path, sep = "\t", quote = FALSE, row.names = TRUE)
  
  return(results_df)
  
}


compute_result_beta <- function(i, mcmc_file, epsilon, alpha, beta) {
  
  mig_rate = mcmc_file[, i]
  
  prior_prop = pbeta(epsilon, shape1 = alpha, shape2 = beta) # proportion of prior < epsilon
  
  posterior_prop = mean(mig_rate < epsilon) # proportion of posterior < epsilon
  
  result = prior_prop / posterior_prop
  
  return(result)
}

BF_Beta <- function(alpha, beta, epsilon, mcmc_file, column_index) {
  
  # Check for valid alpha, beta, epsilon
  if(alpha <= 0 || beta <= 0 || epsilon <= 0) {
    
    stop("The parameters 'alpha', 'beta' and 'epsilon' must be non-zero.")
    
  }

  # Check if any indices are out-of-bounds
  if (!all(column_index %in% 1:ncol(mcmc_file))) {
    
    stop("One or more column indices are out of bounds.")
    
  }
  
  
  # Compute result for each column index and store in a data frame
  results_df = data.frame(Bayes_Factor = sapply(column_index, compute_result_beta, mcmc_file, epsilon, alpha, beta))
  
  # Generate the row names
  row_names <- paste0(colnames(mcmc_file[column_index]), " (", as.character(epsilon), ")")
  
  # Change the ".." to "->". This is specific to the BPP convention of naming columns
  row_names <- sapply(row_names, function(row) gsub("\\.\\.", " -> ", row))
  
  # Assign the row names to the results_df
  rownames(results_df) <- row_names
  
  # Define output file path
  output_file_path <- paste0("BF_Beta_results_", Sys.Date(), "_", ".txt")

  # Write the results to a text file
  write.table(results_df, file = output_file_path, sep = "\t", quote = FALSE, row.names = TRUE)
  
  return(results_df)
  
}



# Call the appropriate function based on the function_name argument

if (function_name == "BF_Gamma") {

  cat("\n")
  print(BF_Gamma(alpha, beta, epsilon, mcmc_file, column_index))
  cat("\n")

} else if (function_name == "BF_Beta") {

  cat("\n")
  print(BF_Beta(alpha, beta, epsilon, mcmc_file, column_index))
  cat("\n")

} else {
 
  cat("\n Unknown function:", function_name)
  cat("\n")

}


