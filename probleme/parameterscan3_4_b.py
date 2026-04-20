
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
repertoire = '/Users/ella/Documents/EPFL/BA4/PHYSNUM/Exercise3_student/probleme'
#repertoire = '/Users/eldidi/Desktop/physnum/physnum_ex3/physnum_ex3/probleme'
executable = '/engine'
input_filename = 'configuration.in.example' # Strictly no longer needed, but we keep it for now to avoid having to change the code in engine.cpp


input_parameters = {
    'numBodies': 2,
    'timeScheme': 0, # 0 pour dt fixe, 1 pour schema adaptatif
    'sampling': 1, # mettre 0 pour ecrire tout les pas
    #'tEnd': 172800,  # 2 jours
    'tEnd' : 2.333e+6 * 10, # 27 jours.
    'dt': 1,
    'tolerance': 1e-6,
    'G': 6.6740e-11,

    #atmosphere
    'useAtmosphere' : 0, # 0 si pas d'athmosphere, 1 si atm
    'roh0' : 1.2,
    'atmosphereScale': 7238.2,
    'dragArea': 19.7923475393,       #je crois que c'est l'aire de l'ombre de dragBody
    'dragCoefficient': 0.3,
    'dragBody': 1,  #indice du corps qui rentre dans l'atm
    'dragCenterBody': 0, #indice du corps qui a une atm      !! indices commencent a 0

    # Terre
    'm1': 5.972e24,
    'r1': 6378.1e3, # en mètres
    'x1': -4.67e6,
    'y1': 0.0, # mètres
    'vx1': 0.0,
    'vy1': -12.4,

    #Lune
    'm2': 7.348e22,
    'r2': 1737.4e3,
    'x2': 3.8e8,
    'y2': 0.0,  # mètres
    'vx2' : 0.0,
    'vy2': 1010,

    # Artémis
    #'m2': 8500,
    #'r2': 2.26, # rayon de la sonde
    #'x2': 314159e3, # en mètres
    #'y2': 0, # mètres
    #'vx2': -1178.0, # v_circulaire = sqrt(GM/r) ≈ 1018 m/s
    #'vy2': 226.0, #en m/s
    #on peut ajouter d'autre corps de la même manière !!!!il faut qu'il y en ai autant que numbodies!!!!
}


# -------------------------------------------------

paramstr = 'timeScheme' # The parameter to scan, must be one of the keys in input_parameters

variable_array = [1]
print(variable_array)

outstr = f"syst_gravitationnel_dt{input_parameters['timeScheme']:.2g}_atm{input_parameters['useAtmosphere']:.2g}_nBodies{input_parameters['numBodies']:.2g}_q34a"   # CHANGER LA FIN SELON QUESTION

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

    lc = LineCollection(segments, cmap='plasma')
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
# Plot 1 : position terre artemis
# ============================================================


#%%

n = len(datasets)
ncols = min(3, n)
nrows = math.ceil(n / ncols)

fig, axarr = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows))
axes = np.array(axarr).reshape(-1)

for j in range(n, len(axes)):
    fig.delaxes(axes[j])


for i, (data, pval) in enumerate(zip(datasets, param_values)):
    ax = axes[i]

    t   = data[:, 0]
    x_T = data[:, 1]
    y_T = data[:, 2]
    x_L = data[:, 5]
    y_L = data[:, 6]
    Em = data[:, 9]
    p = data[:, 10]
    d_TL = data[:, 11]

    ax.plot(x_T, y_T, color='steelblue', lw=1.5, label='Terre')
    ax.plot(x_L, y_L, color='tomato',    lw=1.5, label='Lune')
    ax.plot(x_T[0], y_T[0], 'bo', ms=5)
    ax.plot(x_L[0], y_L[0], 'ro', ms=5)

    ax.ticklabel_format(style='sci', scilimits=(0,0), axis='both')
    ax.set_title(f"{param_name} = {pval:.3g}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect('equal')
    ax.legend(fontsize=7)

    print(t[-1])

fig.suptitle("Position Terre & Lune", fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, "position_terre_lune_3_4_b.png"), dpi=300)


plt.show()

plt.figure()
plt.plot(t, Em)
plt.xlabel('Temps (s)')
plt.ylabel('Energie mécanique')
plt.grid()

plt.figure()
plt.plot(t, d_TL)
plt.xlabel('Temps (s)')
plt.ylabel('Distance Terre-Lune (m)')
plt.grid()

plt.figure()
plt.plot(t, p)
plt.xlabel('Temps (s)')
plt.ylabel('Norme de la qtt de mouv')
plt.grid()

plt.show()

print("Variation relative distance :", (np.max(d_TL) - np.min(d_TL)) / np.mean(d_TL))
    # %%
