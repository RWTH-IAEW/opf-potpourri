import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

xls = pd.ExcelFile("/home/bengisu/potpourri/potpourri/results/results_multiperiod_exact.xlsx")
df = pd.read_excel(xls, 'battery')
Namen = df['name'].unique().tolist()
Namen = ["LV1.101_Load_2_bat"]
fig, ax = plt.subplots(figsize=(10, 6))
for name in Namen:
    ax.plot(df[df['name']==name]['pChar(MW)'].to_numpy(), linewidth=2, label="pChar_exact : " + name, linestyle = "dashdot", color = "green")
    ax.plot(df[df['name']==name]['pDis(MW)'].to_numpy(), linewidth=2, label="pDis_exact: " + name, linestyle = "dashdot", color = "red")
    

xls = pd.ExcelFile("/home/bengisu/potpourri/potpourri/results/results_multiperiod_simplified.xlsx")
df = pd.read_excel(xls, 'battery')
#Namen = df['name'].unique().tolist()
#Namen = ["LV1.101_Load_1_bat"]

for name in Namen:
    ax.plot(df[df['name']==name]['pChar(MW)'].to_numpy(), linewidth=6, label="pChar_simplified : " + name, linestyle = "dotted", color = "green", alpha=0.4)
    ax.plot(df[df['name']==name]['pDis(MW)'].to_numpy(), linewidth=6, label="pDis_simplified: " + name, linestyle = "dotted", color = "red", alpha=0.4)
    
xls = pd.ExcelFile("/home/bengisu/potpourri/potpourri/results/results_multiperiod_extended.xlsx")
df = pd.read_excel(xls, 'battery')
#Namen = df['name'].unique().tolist()
#Namen = ["LV1.101_Load_1_bat"]

for name in Namen:
    ax.plot(df[df['name']==name]['pChar(MW)'].to_numpy(), linewidth=2, label="pChar_extended : " + name, linestyle = "dashed", color = "green")
    ax.plot(df[df['name']==name]['pDis(MW)'].to_numpy(), linewidth=2, label="pDis_extended: " + name, linestyle = "dashed", color = "red")


xls = pd.ExcelFile("/home/bengisu/potpourri/potpourri/results/results_multiperiod_relaxed.xlsx")
df = pd.read_excel(xls, 'battery')
#Namen = df['name'].unique().tolist()
#Namen = ["LV1.101_Load_1_bat"]

for name in Namen:
    ax.plot(df[df['name']==name]['pChar(MW)'].to_numpy(), linewidth=2, label="pChar_relaxed : " + name, linestyle = "dashed", color = "green")
    ax.plot(df[df['name']==name]['pDis(MW)'].to_numpy(), linewidth=2, label="pDis_relaxed: " + name, linestyle = "dashed", color = "red")
    



# Show the plot
#plt.tight_layout()
plt.legend()
ax.set_ylabel('[MW]')
plt.show()