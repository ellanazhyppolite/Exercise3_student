
#%%


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os
import glob
import re
import math
import subprocess

# Parameters
repertoire = '/Users/ella/Documents/EPFL/BA4/PHYSNUM/physnum_ex3/probleme'
executable = '/engine'
input_filename = 'configuration.in.example' # Strictly no longer needed, but we keep it for now to avoid having to change the code in engine.cpp


input_parameters = {
    'numBodies': 2,
    'timeScheme': 0, # 0 pour dt fixe, 1 pour schema adaptatif
    'sampling': 1,
    'tEnd': 10,
    'dt': 0.01,
    'tolerance': 1e-6,
    'G': 1.0,

    #atmosphere
    'useAtmosphere' : 0, # 0 si pas d'athmosphere, 1 si atm
    'roh0' : 0,
    'atmosphereScale': 1.0,
    'dragArea': 0.01,       #je crois que c'est l'aire de l'ombre de dragBody
    'dragCoefficient': 1.0,
    'dragBody': 1,  #indice du corps qui rentre dans l'atm
    'dragCenterBody': 0, #indice du corps qui a une atm      !! indices commencent a 0

    # corps 1
    'm1': 1.0,
    'r1': 0.1,
    'x1': 0.0,
    'y1': 0.0,
    'vx1': 0.0,
    'vy1': 0.0,

    # corps 2
    'm2': 1.0,
    'r2': 0.1,
    'x2': 1.0,
    'y2': 0.0,
    'vx2': 0.0,
    'vy2': 1.0

    #on peut ajouter d'autre corps de la même manière !!!!il faut qu'il y en ai autant que numbodies!!!!
}

# -------------------------------------------------

# Updated from last time, the code below can now be used to scan any parameter, just make sure to update the paramstr and the variable_array accordingly

paramstr = 'dt' # The parameter to scan, must be one of the keys in input_parameters


variable_array = [200,100,50,20,10,1]  # Example values for the parameter scan
# vu que je commprends pas ce qui se passe, mettre la valeur de tEnd ici aussi please

input_parameters['tEnd'] = 172800 #2 jours

outstr = f"syst_gravitationnel_atm_{input_parameters['useAtmosphere']:.2g}_nBodies_{input_parameters['numBodies']:.2g}_q32a"   # CHANGER LA FIN SELON QUESTION

# -------------------------------------------------
# Create output directory (2 significant digits)
# -------------------------------------------------
#outdir = f"Scan_{paramstr}_{outstr}"
outdir = outstr
os.makedirs(outdir, exist_ok=True)
print("Saving results in:", outdir)


for i in range(len(variable_array)):

    # Copy parameters and overwrite scanned one
    params = input_parameters.copy()
    params[paramstr] = variable_array[i]

    output_file = f"{outstr}_{paramstr}_{variable_array[i]}.txt"
    output_path = os.path.join(outdir, output_file)

    # Build parameter string
    param_string = " ".join(f"{k}={v:.15g}" for k, v in params.items())

    cmd = (
        f"{repertoire}{executable} {input_filename} "
        f"{param_string} output={output_path}"
    )

    print(cmd)
    subprocess.run(cmd, shell=True)
    print("Done.")


output_file = f"{outstr}.txt"

# ============================================================
# USER SETTINGS
# ============================================================

folder = outdir  


## a quoi ca sert
plot_layout = {
    "theta_time": True,
    "phase_space": False,
    "energy": True,
    "real_space": False,
    "power": True,
    "energy_balance": True
}

# ============================================================
# Output folder
# ============================================================

fig_dir = os.path.join(folder, "figures")
os.makedirs(fig_dir, exist_ok=True)

# ============================================================
# Scan files
# ============================================================


files = sorted(glob.glob(os.path.join(folder, "*.txt")))

datasets = []
param_values = []
param_name = None

for f in files:

    name = os.path.basename(f)      # remove path
    name = name[:-4]                # remove ".txt"

    parts = name.split("_")

    param_name = parts[-2]          # scanned parameter
    value = float(parts[-1])        # parameter value

    data = np.loadtxt(f)

    datasets.append(data)
    param_values.append(value)

print(f"Found {len(datasets)} datasets.")

# Sort datasets
order = np.argsort(param_values)
param_values = np.array(param_values)[order]
datasets = [datasets[i] for i in order]


# ============================================================
# Helper: colored line
# ============================================================

def colored_line(x, y, t, ax, vmin=None, vmax=None):

    points = np.array([x, y]).T.reshape(-1,1,2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, cmap='viridis')
    lc.set_array(t)

    if vmin is not None and vmax is not None:
        lc.set_clim(vmin, vmax)

    lc.set_linewidth(2)

    ax.add_collection(lc)
    ax.autoscale()

    return lc

# ============================================================
# Shared time range for color maps
# ============================================================

tmin = min(data[:,0].min() for data in datasets)
tmax = max(data[:,0].max() for data in datasets)

# ============================================================
# Axis layout helper
# ============================================================

def get_axes(plot_key, title):

    if plot_layout[plot_key]:

        fig, ax = plt.subplots()
        axes = [ax]*len(datasets)

    else:

        n = len(datasets)
        ncols = min(3, n)
        nrows = math.ceil(n/3)

        fig, axarr = plt.subplots(nrows, ncols,
                                  figsize=(5*ncols,4*nrows))

        axes = np.array(axarr).reshape(-1)

        for j in range(n, len(axes)):
            fig.delaxes(axes[j])

        axes = axes[:n]

    fig.suptitle(title)

    return fig, axes


cmap = plt.get_cmap("tab10")


# ============================================================
# Plot 1 : Position en fonction du temps
# ============================================================



# %%
