import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
def test(solver):
    # Load the true data
    df_true = pd.read_csv('./SavedResult/solvercompare/VCABM_high_run1.csv')
    # Prepare the interpolators
    interpolators = {}
    min_time, max_time = df_true['Time'].min(), df_true['Time'].max()
    for column in df_true.columns[1:]:
        interpolators[column] = interp1d(df_true['Time'], df_true[column], kind='linear', bounds_error=False, fill_value="extrapolate")

    # Load the run data
    #solver = "DP5"
    df_run = pd.read_csv(f'./SavedResult/solvercompare/{solver}.csv')
    times = df_run['Time']

    # Evaluate the interpolators at the valid solver times
    valid_times = times[(times >= min_time) & (times <= max_time)]
    true_ys = np.array([interpolators[col](valid_times) for col in interpolators]).T
    run_ys = df_run.loc[df_run['Time'].isin(valid_times), df_true.columns[1:]].values
    rmse = np.sqrt(np.mean((run_ys - true_ys)**2))
    print(f"RMSE: {rmse}")
    # Make sure errors array starts with an initial value (e.g., 0) to match times[0]
    # Calculate the differences for each time point and each variable
    differences = run_ys - true_ys
    
    # Define sets of colors and line styles
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # 7 basic color abbreviations
    colors = [
    (16/256, 70/256, 128/256),   
    (49/256, 124/256, 183/256),    
    (109/256, 173/256, 209/256),     
    (182/256, 215/256, 232/256),       
    #(233/256, 241/256, 244/256),       
    #(251/256, 227/256, 213/256),       
    (246/256, 178/256, 147/256),
    (220/256, 109/256, 87/256),
    (183/256, 34/256, 48/256),
    #(109/256, 48/256, 31/256)
]
    line_styles = ['-', '--', '-.', ':', (0, (1, 10)), (0, (5, 10)), (0, (3, 5, 1, 5)), 
                (0, (5, 1)), (0, (3, 10, 1, 10)), (0, (3, 1, 1, 1))]  # 10 different line styles

    # Create a list of combined color and line style for up to 70 combinations
    from itertools import product
    style_combinations = list(product(colors, line_styles))

    # Plotting configuration
    plt.figure(figsize=(8, 3))

    # Example data generation

    # Plot each variable's differences using a unique style combination
    for i, column in enumerate(df_true.columns[1:70]):
        plt.plot(valid_times, differences[:, i], label=f'Var {i+1}', color=style_combinations[i][0], linestyle=style_combinations[i][1])

    plt.title(f'Time vs Differences for 70 Variables RMSE: {rmse}')
    plt.xlabel('Time')
    plt.ylabel('Difference')
    #plt.legend()
    plt.grid(True, which="both", ls="-")
    plt.savefig(f'./SavedResult\solvercompare\{solver}_tc{len(times)}_rsme{rmse}.png', dpi=500)
    
    # Plotting configuration
    plt.figure(figsize=(8, 3))

    # Example data generation

    # Plot each variable's differences using a unique style combination
    for i, column in enumerate(df_true.columns[1:70]):
        plt.plot(valid_times, differences[:, i], label=f'Var {i+1}', color=style_combinations[i][0], linestyle=style_combinations[i][1])

    plt.title(f'Time vs Differences for 70 Variables (Log Scale) RMSE: {rmse}')
    plt.xlabel('Time')
    plt.ylabel('Difference')
    plt.yscale('log')
    #plt.legend()
    plt.grid(True, which="both", ls="-")
    plt.savefig(f'./SavedResult\solvercompare\{solver}_tc{len(times)}_rsme{rmse}_log.png', dpi=500)
    plt.show()
strs=['VCAB3_run1','VCAB4_run1','VCAB5_run1',
    'VCABM3_run1','VCABM4_run1','VCABM5_run1','VCABM_run1','VCABM_high12_run1']
for str in strs:
    test(str)