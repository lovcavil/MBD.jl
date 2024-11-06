import os
import pandas as pd

def csv_folder_to_excel_with_styles(folder_path, output_excel):
    combined_df = pd.DataFrame()  # Initialize an empty DataFrame to hold all CSV data
    
    # Loop through all files in the folder
    for filename in os.listdir(folder_path):
        # Only process CSV files
        if filename.endswith(".csv"):
            # Create the full file path
            file_path = os.path.join(folder_path, filename)
            
            # Read the CSV file
            df = pd.read_csv(file_path)
            
            # Add a column to indicate the source CSV file for easier tracking in the combined sheet
            df['Source'] = filename
            
            # Append the data to the combined DataFrame
            combined_df = pd.concat([combined_df, df], ignore_index=True)
    
    # Write the combined DataFrame to a single sheet in the Excel file
    with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
        # Write the combined data into one sheet
        combined_df.to_excel(writer, sheet_name='Combined', index=False)
        
        # Get the XlsxWriter workbook and worksheet objects to apply styles
        workbook  = writer.book
        worksheet = writer.sheets['Combined']
        
        # Apply autofilter to all columns based on the first row
        worksheet.autofilter(0, 0, combined_df.shape[0], combined_df.shape[1] - 1)
        
        # Conditional formatting to apply color bars to time columns
        # Assuming 'Avg Time (s)' and 'All Times' contain time-related data
        if 'Avg Time (s)' in combined_df.columns:
            avg_time_col_idx = combined_df.columns.get_loc('Avg Time (s)')
            worksheet.conditional_format(1, avg_time_col_idx, combined_df.shape[0], avg_time_col_idx, 
                                         {'type': 'data_bar', 'bar_color': '#FFA07A'})  # Light salmon color

        if 'All Times' in combined_df.columns:
            all_times_col_idx = combined_df.columns.get_loc('All Times')
            worksheet.conditional_format(1, all_times_col_idx, combined_df.shape[0], all_times_col_idx, 
                                         {'type': 'data_bar', 'bar_color': '#87CEFA'})  # Light sky blue color
    
    print(f"All CSV files from {folder_path} have been combined into {output_excel} in a single sheet with filters and styles applied.")

# Example usage
folder_path = "D:\OneDrive\Articles/10.Working\[D21][20211009]ContactMechanics\MBD.jl\MBD.jl\src\dejl option comp3/2024-10-27-17-37-11"  # replace with the actual folder path
output_excel = "D:\OneDrive\Articles/10.Working\[D21][20211009]ContactMechanics\MBD.jl\MBD.jl\src\dejl option comp3/2024-10-27-17-37-11/combined_output_with_styles.xlsx"  # replace with desired Excel output filename

csv_folder_to_excel_with_styles(folder_path, output_excel)
