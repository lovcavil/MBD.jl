import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
def test(path,solver,color):
    df_run=pd.read_csv(f'{path}{solver}.csv')
    times=df_run['Time']
    
    # Plotting configuration
    # Example data generation
    time_steps = times.diff()
    # Plot each variable's differences using a unique style combination
    plt.plot(times[1:], time_steps[1:], color=color)

def run(path):
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
    (109/256, 48/256, 31/256),
    (109/256,109/256,109/256)
]
    
    #DP5=Dict(:alg => OrdinaryDiffEq.DP5())
    #Tsit5=Dict(:alg => OrdinaryDiffEq.Tsit5())
    #RKO65 =Dict(:alg => OrdinaryDiffEq.RKO65()) #fix
    #TanYam7 =Dict(:alg => OrdinaryDiffEq.TanYam7())
    #DP8=Dict(:alg => OrdinaryDiffEq.DP8())
    #TsitPap8=Dict(:alg => OrdinaryDiffEq.TsitPap8())
    
    
    str=["contact_off_Euler_2024-06-06_13-17-51_run1",
         "contact_off_Euler_2024-06-06_14-53-25_run1"]
    plt.figure(figsize=(8, 6))
    i=1
    for s in str:
        test(path,s,colors[i])
        i=i+1
        
    plt.title('Time vs Differences for 70 Variables (Log Scale)')
    plt.xlabel('Time')
    plt.ylabel('Difference')
    plt.legend(str)
    plt.grid(True, which="both", ls="-")    
    plt.savefig(f'{path}timestep_diff_contact_off_Euler.png', dpi=500)

    plt.show()

path='./csv/solvercompare/'
run(path)