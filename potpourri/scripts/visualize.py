import pandas as pd
from matplotlib import pyplot as plt
import numpy as np



xls = pd.ExcelFile("/home/bengisu/potpourri/potpourri/results/results_multiperiod_2.xlsx")
df = pd.read_excel(xls, 'summary')

time_periods = df["Time Periods"].to_numpy()

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(df["Solar Generation(MW)"].to_numpy(), label="Solar Generation(MW)",linewidth=2)

ax.plot(df["Conventional generation (MW)"].to_numpy(), label="Conventional generation (MW)",linewidth=2)

ax.plot(df["Demand (MW)"].to_numpy(), label="Demand (MW)",linewidth=2)

# Set x-axis tick labels
ax.set_xticks(range(len(time_periods)))
ax.set_xticklabels(time_periods, rotation=45)

ax.set_xlim(0)
ax.set_ylim(0)

ax.grid(True)

plt.legend()
ax.set_ylabel('[MW]')
plt.show()
'''
'''
#xls = pd.ExcelFile("/home/bengisu/potpourri/potpourri/results/results_multiperiod.xlsx")
#df = pd.read_excel(xls, 'battery')
#Namen = df['name'].unique().tolist()
linestyle = ["dashdot", "dotted", "solid"]
i = 0
fig, ax = plt.subplots(figsize=(10, 6))
paths = ["/home/bengisu/potpourri/potpourri/results/results_multiperiod_exact.xlsx",
         "/home/bengisu/potpourri/potpourri/results/results_multiperiod_simplified.xlsx",
         "/home/bengisu/potpourri/potpourri/results/results_multiperiod_relaxed.xlsx"]
         #"/home/bengisu/potpourri/potpourri/results/results_multiperiod_extended.xlsx"]
for path in paths:
    xls = pd.ExcelFile(path)
    df = pd.read_excel(xls, 'battery')
    name = "LV1.101_Load_2_bat"
    model = path.rsplit("/",1)[1]
    line = linestyle[i]
    i = i+1
    ax.plot(df[df['name']==name]['pChar(MW)'].to_numpy(), linewidth=2, label="pChar : " + model, linestyle = line, color = "green")
    ax.plot(df[df['name']==name]['pDis(MW)'].to_numpy()*(-1), linewidth=2, label="pDis: " + model, linestyle = line, color = "red")

    

# Show the plot
#plt.tight_layout()
plt.legend()
ax.set_ylabel('[MW]')
plt.show()