import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.signal import savgol_filter

def load_test(file):
    # Define the column names (equivalent to VariableNames in MATLAB)
    column_names = [
        "StrainB1_1", "StrainB2_1", "StrainB2_2", "StrainB2_3", "StrainB2_4", "StrainB3_1", "StrainB3_2", 
        "StrainB3_3", "StrainB3_4", "accA1_X", "accA1_Y", "accA1_Z", "accA2_X", "accA2_Y", "accA2_Z", 
        "accA3_X", "accA3_Y", "accA3_Z", "accA4_X", "accA4_Y", "accA4_Z", "accA5_X", "accA5_Y", "accA5_Z", 
        "accA6_X", "accA6_Y", "accA6_Z", "accA7_X", "accA7_Y", "accA7_Z", "canpress"
    ]
    # Specify the column types (equivalent to VariableTypes in MATLAB)
    column_types = {name: 'float64' for name in column_names}
    
    # Read the file, starting from line 27 and using space and comma as delimiters
    df = pd.read_csv(file, delim_whitespace=True, skiprows=26, names=column_names, dtype=column_types, engine='python')
    
    return df

def post_test(fn):
    test = load_test(fn)
    delta = 1.9531250e-3
    t = np.arange(delta, delta * (len(test) + 1), delta)
    
    # Add time vector `t` to the DataFrame
    test['t'] = t
    
    return test

# Initialize variables

label_data = {
    (24, 1): 0.86, (24, 2): 0.80, (24, 3): 0.82, (24, 4): 1.00, (24, 5): 1.04, (24, 6): 1.20, (24, 7): 1.21,
    (25, 1): 1.21, (25, 2): 1.17, (25, 3): 1.06, (25, 4): 0.94, (25, 5): 0.89, (25, 6): 0.68, (25, 7): 0.95, (25, 8): 0.73,
    (26, 1): 0.80, (26, 2): 0.87, (26, 3): 0.87, (26, 4): 0.92, (26, 5): 0.99, (26, 6): 1.17, (26, 7): 0.80,
    (27, 1): 0.69, (27, 2): 0.74, (27, 3): 0.92, (27, 4): 1.05, (27, 5): 1.00,
    (28, 1): 1.16, (28, 2): 0.82, (28, 3): 0.90, (28, 4): 1.04,
    (29, 1): 0.56, (29, 2): 0.67, (29, 3): 0.80, (29, 4): 0.78, (29, 5): 0.77, (29, 6): 0.95, (29, 7): 0.95, (29, 8): 1.03,
    (30, 1): 0.95, (30, 2): 1.14, (30, 3): 0.73, (30, 4): 0.92, (30, 5): 0.97,
    #(31, 1): 1.22, (31, 2): 1.27, (31, 3): 1.32, (31, 4): 1.35, (31, 5): 1.38, (31, 6): 1.54,
    #(32, 1): 1.73, (32, 2): 1.66, (32, 3): 1.77,  (32, 5): 1.68,#(32, 4): 1.87,
}

# Calculate min and max values for normalization
min_val = min(label_data.values())
max_val = max(label_data.values())

# Normalize the label data to the range [0, 1]
normalized_label_data = {key: (val - min_val) / (max_val - min_val) for key, val in label_data.items()}

dforcesme = {}
dforcesme0 = []
count = 1

# Create the figure and set up the subplot
plt.figure(facecolor='white')
ax = plt.subplot(1, 1, 1)
ax.grid(True)

# Use the normalized data directly with the colormap
cmap = cm.get_cmap('plasma')
cmap = cm.get_cmap('autumn_r')
# Iterate over the ranges
for i1 in range(24, 31):
    for i2 in range(1, 9):
        if i1 == 24 and i2 == 6:
            continue
        if i1 == 31 and i2 == 1:
            continue
        if i1 == 32 and i2 == 4:
            continue
        b = f"D:/D11.Proj/A009.BYD_R2/20220710/DATA_PROC/CUT/CLOSE_NORMAL/MR_door_{i1}{i2}.asc"
        sub_dir_name = f"{i1}_{i2}"

        if (i1, i2) in normalized_label_data:
            norm_val = normalized_label_data[(i1, i2)]
            color = cmap(norm_val)  # Map the normalized value to a color

            if os.path.isfile(b):
                # Load the data
                test = post_test(b)

                x = test['t'].to_numpy()
                y = test['canpress'].to_numpy()

                aaa = y > 10
                V1 = np.argmax(aaa) if np.any(aaa) else None
                M = np.max(y)
                V2 = np.argmax(y)
                V3 = len(y) - np.argmax(aaa[::-1]) - 1 if np.any(aaa) else None

                if V1 is not None and V3 is not None:
                    dforcesme[(i1, i2)] = [x[V1], x[V2], x[V3], M]
                    dforcesme0.append([x[V1], x[V2], x[V3], M])
                    count += 1

                    # Apply smoothing to y
                    # Ensure window_length is less than or equal to the length of y
                    window_length = min(17, len(y) if len(y) % 2 != 0 else len(y) - 1)
                    if window_length < 5:
                        window_length = 5  # Minimum window length
                    elif window_length % 2 == 0:
                        window_length += 1  # Make it odd

                    y_smoothed = savgol_filter(y, window_length=window_length, polyorder=3)

                    # Plotting with color according to the normalized value
                    st = 1
                    en = min(V3 + 40, len(x))  # Adjust end index to avoid out of bounds
                    ax.plot(x[st:en], y_smoothed[st:en], color=color)

# Add a colorbar to the figure using the original normalization
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=min_val, vmax=max_val))
sm.set_array([])  # Only needed for ScalarMappable
plt.colorbar(sm, ax=ax, label='Velocity')
ax.set_xlabel('t (s)')
ax.set_ylabel('Driving Forces (N)')
ax.set_title('Smoothed Pressure Over Time with Velocity Color Mapping')
plt.show()

# Convert dforcesme0 to a numpy array for further processing if needed
dforcesme0 = np.array(dforcesme0)
