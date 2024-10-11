import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
def test(solver,color):
    df_run=pd.read_csv(f'./SavedResult/solvercompare/{solver}.csv')
    times=df_run['Time']
    
    # Plotting configuration
    # Example data generation
    time_steps = times.diff()
    # Plot each variable's differences using a unique style combination
    plt.plot(times[1:], time_steps[1:], color=color)

def run():
    colors = [
    (16/256, 70/256, 128/256),   
    (49/256, 124/256, 183/256),    
    (109/256, 173/256, 209/256),     
    (182/256, 215/256, 232/256),       
    (233/256, 241/256, 244/256),       
    #(251/256, 227/256, 213/256),       
    (246/256, 178/256, 147/256),
    (220/256, 109/256, 87/256),
    (183/256, 34/256, 48/256),
    (109/256, 48/256, 31/256),
    (109/256,109/256,109/256)
]
    str=["Midpoint_run1","Heun_run1","Ralston_run1","RK4_run1",
         "BS3_run1","OwrenZen3_run1","OwrenZen4_run1","OwrenZen5_run1",'VCABM_high12_run1']
    #str=['Midpoint_run1','VCABM_high12_run1']
    plt.figure(figsize=(8, 6))
    i=1
    for s in str:
        test(s,colors[i])
        i=i+1
        
    plt.title('Time vs Differences for 70 Variables (Log Scale)')
    plt.xlabel('Time')
    plt.ylabel('Difference')
    plt.legend(str)
    plt.grid(True, which="both", ls="-")    
    plt.savefig(f'./SavedResult\solvercompare/timestep_diff_ERK2.png', dpi=500)

    plt.show()    
run()