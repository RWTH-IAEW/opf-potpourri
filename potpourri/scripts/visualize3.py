
import pandas as pd
from matplotlib import pyplot as plt

# Create a single figure and axes object
fig, ax = plt.subplots(figsize=(10, 6))

# Define the Excel files and sheets
files = [
    ("/home/bengisu/potpourri/potpourri/results/results_multiperiod_exact.xlsx", "battery"),
    ("/home/bengisu/potpourri/potpourri/results/results_multiperiod_simplified.xlsx", "battery"),
    ("/home/bengisu/potpourri/potpourri/results/results_multiperiod_extended.xlsx", "battery")
]

# Iterate over each file and plot the data
for file_idx, (file, sheet) in enumerate(files):
    xls = pd.ExcelFile(file)
    df = pd.read_excel(xls, sheet)
    Namen = df['name'].unique().tolist()

    for name_idx, name in enumerate(Namen):
        pChar_label = f"pChar{file_idx+1}_{name_idx+1}_{sheet}"
        pDis_label = f"pDis{file_idx+1}_{name_idx+1}_{sheet}"
        
        ax.plot(df[df['name']==name]['pChar(MW)'].to_numpy(), linewidth=2, label=pChar_label, linestyle="dashdot")
        ax.plot(df[df['name']==name]['pDis(MW)'].to_numpy(), linewidth=2, label=pDis_label, linestyle="dotted")

# Set labels and legend
ax.set_ylabel('[MW]')
plt.legend()

# Show the plot
plt.show()



