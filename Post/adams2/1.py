import pandas as pd
import re

# Define the file content as a string for demonstration purposes
file_content = """
BEGIN
COLUMNWIDTH = [  15,  15]
COLUMNOFFSET = [   0,   2]
CHANNELNAME = ['CPU FUNCTION_MEA_CPU','HSIZE FUNCTION_MEA_HSIZE']
LENGTH = [      1528,      1528]
UNIT = ['NO UNITS','NO UNITS']
MINIMUM = [-3.5625000e+000,-1.0000000e-003]
MAXIMUM = [+3.5625000e+000,+1.0000000e-003]
START = [+0.0000000e+000,+0.0000000e+000]
DELTA = [+1.0000000e-003,+1.0000000e-003]
# descriptive data
END
+3.1217558e-002  +0.0000000e+000  
+4.6771951e-002  +3.3332314e-004  
+4.6771951e-002  +7.3332927e-004  
+4.6771951e-002  +5.3331093e-004  
+4.6771951e-002  +8.4999390e-004  
+4.6771951e-002  +2.9998168e-004  
"""

# Parse the file content
lines = file_content.strip().split('\n')

# Extract CHANNELNAME
channel_line = next(line for line in lines if line.startswith('CHANNELNAME'))
channel_names = re.findall(r"'(.*?)'", channel_line)

# Extract the data lines (after the END marker)
data_lines = lines[lines.index('END') + 1:]

# Remove anything after '&' in the data lines
data_lines = [line.split('&')[0].strip() for line in data_lines]

# Create a DataFrame from the data lines
data = [list(map(float, re.split(r'\s+', line))) for line in data_lines]
df = pd.DataFrame(data, columns=channel_names)


print(df)
