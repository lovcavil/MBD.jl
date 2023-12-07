import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
# dfad = pd.read_csv("1slider/cr.tab", sep='\t')
# dfad = pd.read_csv("1slider/slider.tab", sep='\t')
# dfad = pd.read_csv("1slider/link.tab", sep='\t')
# Read the CSV file
dfjl = pd.read_csv("1slider/data.csv")

# Plotting
def d(labeljl,tad,yad,pname,i0):
    plt.figure(figsize=(16,10))
    plt.plot(dfjl['t'],dfjl[labeljl], label=labeljl)
    plt.plot(tad,yad+i0)
    plt.savefig(f"{pname}/{labeljl}.png")  # Saves the plot as a PNG file
    plt.close()  # Close the plot to free up memory
    
pname="1slider"

dfad = pd.read_csv("1slider/slider.tab", sep='\t')
t,y=dfad['.MODEL_1.Last_Run.PART_4_XFORM.TIME'],dfad['.MODEL_1.Last_Run.PART_4_XFORM.X']
d('x3',t,y,pname,0.2) 


# dfad = pd.read_csv("1slider/cr.tab", sep='\t')
# t,y=dfad['.MODEL_1.Last_Run.cr_XFORM.TIME'],dfad['.MODEL_1.Last_Run.cr_XFORM.Y']
# d('z1',t,y,pname) 


dfad = pd.read_csv("1slider/link.tab", sep='\t')
t,y=dfad['.MODEL_1.Last_Run.PART_5_XFORM.TIME'],dfad['.MODEL_1.Last_Run.PART_5_XFORM.Z']
d('z2',t,y,pname,0.1) 
    