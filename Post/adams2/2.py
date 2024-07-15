import pandas as pd
import re

# Function to read the structured file and create a DataFrame
def read_structured_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Extract CHANNELNAME
    channel_line = next(line for line in lines if line.startswith('CHANNELNAME'))
    channel_names = re.findall(r"'(.*?)'", channel_line)

    # Extract the data lines (after the END marker)
    data_lines = lines[lines.index('END\n') + 1:]

    # Remove anything after '&' in the data lines
    data_lines = [line.split('&')[0].strip() for line in data_lines]

    # Create a DataFrame from the data lines
    data = [list(map(float, re.split(r'\s+', line))) for line in data_lines]
    df = pd.DataFrame(data, columns=channel_names)

    return df

# Read the file and create the DataFrame
file_path = r'V:\tmp\A009.BYD_CONTACTMODEL3\RUN\203.asc'  # Replace with your file path
df = read_structured_file(file_path)
print( df)
