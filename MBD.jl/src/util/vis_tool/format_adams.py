import argparse
import json
import pandas as pd
import os

# Determine the directory where the script is located
script_dir = os.path.dirname(__file__)

# Set up argument parser
parser = argparse.ArgumentParser(description='Process text files into a CSV based on JSON mapping.')
parser.add_argument('--column_mapping_file', type=str, default=os.path.join(script_dir, 'column_mapping.json'), help='Path to the column mapping JSON file.')
parser.add_argument('--file_list_file', type=str, default=os.path.join(script_dir, 'file_list.json'), help='Path to the JSON file with the list of text files.')
parser.add_argument('--output_file', type=str, default='adams_output.csv', help='Path to the JSON file with the list of text files.')

# Parse arguments
args = parser.parse_args()

# Load column mapping
with open(os.path.join(script_dir, args.column_mapping_file), 'r') as file:
    column_mapping = json.load(file)

# Load list of text files
with open(os.path.join(script_dir, args.file_list_file), 'r') as file:
    file_list = json.load(file)['files']

# Initialize an empty list to store DataFrames
dfs = []

# Process each file
for file_name in file_list:
    file_path = os.path.join(os.path.dirname(args.file_list_file), file_name)
    # Read the file, skipping initial rows that don't contain data
    data = pd.read_csv(file_path, sep="\t", skiprows=[0,1,2,3,4])

    # Keep only columns that are in the column mapping
    data = data[[col for col in data.columns if col in column_mapping]]

    # Rename columns based on mapping
    data.rename(columns=column_mapping, inplace=True)

    # Add the DataFrame to the list
    dfs.append(data)

# Concatenate all DataFrames horizontally
combined_data = pd.concat(dfs, axis=1)

# Remove duplicate columns after concatenation
combined_data = combined_data.loc[:,~combined_data.columns.duplicated()]

# Write the combined data to a CSV file
output_file = os.path.join(script_dir, args.output_file)
combined_data.to_csv(output_file, index=False)
