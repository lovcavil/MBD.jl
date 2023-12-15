import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
df = pd.read_csv("1slider/slider.tab", sep='\t')

# Plotting
def d(label,pname):
    plt.figure(figsize=(16,10))
    plt.plot(df['.MODEL_1.Last_Run.PART_sl_XFORM.TIME'],df[label], label=label)
    plt.savefig(f"{pname}/{label}.png")  # Saves the plot as a PNG file
    plt.close()  # Close the plot to free up memory
    
    

pname="1slider"
for label in df.columns:
    d(label,pname)


