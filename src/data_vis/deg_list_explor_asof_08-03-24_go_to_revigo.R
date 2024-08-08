# Define the directory paths
dir_path <- "./reports/local_only/deg_result_clustering/go_term_v3"
output_dir_path <- "./reports/local_only/deg_result_clustering/revigo_output"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir_path)) {
  dir.create(output_dir_path)
}

# Get a list of all .txt files in the directory
file_list <- list.files(path = dir_path, pattern = "\\.txt$", full.names = TRUE)

# Function to read the file into a dataframe
read_file <- function(file_path) {
  read.delim(file_path, stringsAsFactors = FALSE)
}

# Function to ensure 'Term' column is of type character
convert_term_to_character <- function(df) {
  df$Term <- as.character(df$Term)
  df
}

# Function to modify the 'Term' column to keep only the strings in front of '~'
modify_term_column <- function(df) {
  df$Term <- sapply(df$Term, function(x) sub("~.*", "", x))
  df
}

# Function to sort the dataframe by the 'FDR' column
sort_by_fdr <- function(df) {
  df[order(df$FDR), ]
}

# Function to select only the 'Term' and 'FDR' columns
select_columns <- function(df) {
  df[, c("Term", "FDR")]
}

# Function to save the resulting dataframe as a new file
save_dataframe <- function(df, original_file_path, output_dir_path) {
  new_file_name <- paste0("revigo_", basename(original_file_path))
  new_file_path <- file.path(output_dir_path, new_file_name)
  write.table(df, file = new_file_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  new_file_name
}

# Main function to process each file
process_file <- function(file_path) {
  df <- read_file(file_path)
  df <- convert_term_to_character(df)
  df <- modify_term_column(df)
  df <- sort_by_fdr(df)
  df <- select_columns(df)
  save_dataframe(df, file_path, output_dir_path)
}

# Apply the function to all files and store the result in a list of new file names
new_file_list <- sapply(file_list, process_file)

# Print the names of the new files
print(new_file_list)
