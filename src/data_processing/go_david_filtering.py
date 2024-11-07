import pandas as pd
import os
import sys
from glob import glob

# Check if the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python3 process_files.py <input_dir> <output_dir>")
    sys.exit(1)

# Get input and output directories from command-line arguments
input_dir = sys.argv[1]
output_dir = sys.argv[2]
os.makedirs(output_dir, exist_ok=True)

# Get a list of all tab-delimited files in the input directory
file_paths = glob(os.path.join(input_dir, '*.txt'))

# Process each file
for file_path in file_paths:
    # Load the file
    df = pd.read_csv(file_path, delimiter='\t')
    
    # 1) Split the 'Term' column at '~' to create 'GO_ID' and 'Description' columns
    df[['GO_ID', 'Description']] = df['Term'].str.split('~', expand=True)
    
    # 2) Filter rows where FDR < 0.05
    filtered_df = df[df['FDR'] < 0.05]
    
    # 3) Sort by 'Fold Enrichment' in descending order
    filtered_df = filtered_df.sort_values(by='Fold Enrichment', ascending=False)
    
    # Define the output file name based on the original file name
    file_name = os.path.basename(file_path).replace('.txt', '_filtered.txt')
    output_file = os.path.join(output_dir, file_name)
    
    # 4) Save the resulting table to the output directory
    filtered_df.to_csv(output_file, sep='\t', index=False)
    print(f"Processed file saved to: {output_file}")
