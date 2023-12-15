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
    plt.savefig(f"{pname}/{labeljl}-303-ad-jl.png")  # Saves the plot as a PNG file
    plt.close()  # Close the plot to free up memory
    
pname="1slider"

# dfad = pd.read_csv("1slider/slider.tab", sep='\t')
# t,y=dfad['.MODEL_1.Last_Run.PART_sl_XFORM.TIME'],dfad['.MODEL_1.Last_Run.PART_sl_XFORM.X']
# d('x3',t,y,pname,0.2) 

# dfad = pd.read_csv("1slider/slider.tab", sep='\t')
# t,y=dfad['.MODEL_1.Last_Run.PART_sl_XFORM.TIME'],dfad['.MODEL_1.Last_Run.PART_sl_XFORM.VX']
# d('xd3',t,y,pname,0.0) 

# dfad = pd.read_csv("1slider/cr.tab", sep='\t')
# t,y=dfad['.MODEL_1.Last_Run.cr_XFORM.TIME'],dfad['.MODEL_1.Last_Run.cr_XFORM.Y']
# d('z1',t,y,pname) 


# dfad = pd.read_csv("1slider/link.tab", sep='\t')
# t,y=dfad['.MODEL_1.Last_Run.PART_lk_XFORM.TIME'],dfad['.MODEL_1.Last_Run.PART_lk_XFORM.VZ']
# d('z2',t,y,pname,0.1) 
    
    
    

dfmt = pd.read_csv("1slider/t.csv")
dfmx = pd.read_csv("1slider/x3.csv")
# Plotting

def draw_m_res(label,pname):
    plt.figure(figsize=(16,10))
    plt.plot(dfmt['t'],dfmx[label], label=label)
    plt.savefig(f"{pname}/{label}m.png")  # Saves the plot as a PNG file
    plt.close()  # Close the plot to free up memory


for label in dfmx.columns:
    draw_m_res(label,"1slider")    
pass

pname="1slider"


#_______________________________________________________________________
fig=plt.figure(figsize=(16,10))
dfad = pd.read_csv("1slider/slider.tab", sep='\t')
t,y=dfad['.MODEL_1.Last_Run.PART_sl_XFORM.TIME'],dfad['.MODEL_1.Last_Run.PART_sl_XFORM.X']
plt.plot(t,y+0.2, label='x3-adams',marker = "*")
plt.plot(dfjl['t'],dfjl['x3'], label='x3-jl',marker = "o")
plt.plot(dfmt['t'],dfmx['x3'], label='x3-mat',marker = "x")
plt.legend()
plt.savefig(f"{pname}/compare-x3.png")  # Saves the plot as a PNG file
plt.close()    

fig=plt.figure(figsize=(16,10))
dfad = pd.read_csv("1slider/slider.tab", sep='\t')
t,y=dfad['.MODEL_1.Last_Run.PART_sl_XFORM.TIME'],dfad['.MODEL_1.Last_Run.PART_sl_XFORM.VX']
plt.plot(t,y, label='x3d-adams',marker = "*")
plt.plot(dfjl['t'],dfjl['xd3'], label='x3d-jl',marker = "o")
dfmx = pd.read_csv("1slider/x3d.csv")
plt.plot(dfmt['t'],dfmx['x3d'], label='x3d-mat',marker = "x")
plt.legend()
plt.savefig(f"{pname}/compare-xd3.png")  # Saves the plot as a PNG file
plt.close()    