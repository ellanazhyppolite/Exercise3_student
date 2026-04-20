
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
    'numBodies': 3,
    'timeScheme': 1, # 0 pour dt fixe, 1 pour schema adaptatif
    'sampling': 0, ##mettre 0 pour ecrire tout les pas
    'tEnd': 172800*5,  # 2 jours
    #'tEnd' : 2.333e+6 * 10, # 27 jours.
    'dt': 1,
    'tolerance': 1e-6,
    'G': 6.6740e-11,
    'theta_v' : 0.0,
    'q3_5' : 1.0,       ## mettre 1 si on veut faire le swipe sur theta

    #atmosphere
    'useAtmosphere' : 0, ## 0 si pas d'athmosphere, 1 si atm
    'roh0' : 1.2,
    'atmosphereScale': 7238.2,
    'dragArea': 19.7923475393,       #je crois que c'est l'aire de l'ombre de dragBody
    'dragCoefficient': 0.3,
    'dragBody': 1,  #indice du corps qui rentre dans l'atm
    'dragCenterBody': 0, #indice du corps qui a une atm      !! indices commencent a 0

    # Terre
    'm1': 5.972e24,
    'r1': 6378.1e3, # en mètres
    'x1': -4676244.537,
    'y1': 0.0, # mètres
    'vx1': 0.0,
    'vy1': -12.45670891,

    # Artémis
    'm2': 8500,
    'r2': 2.26, # rayon de la sonde
    'x2': 314159e3 - 4676244.537, # en mètres
    'y2': 0, # mètres
    'vx2': 1.2e3,           ## metter ici norme vitesse artemis
    'vy2': 0.0, #en m/s
    
    #Lune
    'm3': 7.348e22,
    'r3': 1737.4e3,
    'x3': 380071755.5,
    'y3': 0.0,  # mètres
    'vx3' : 0.0,
    'vy3': 1012.095037,
    
}


# -------------------------------------------------

paramstr = 'theta_v' # The parameter to scan, must be one of the keys in input_parameters

variable_array = np.linspace(1.6, 2.96, 100)
print(variable_array)

outstr = f"syst_gravitationnel_theta_q35a"   # CHANGER LA FIN SELON QUESTION

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

    param_name = r"$\theta$"          # scanned parameter
    value = float(parts[-1])        # parameter value

    data = np.loadtxt(f)

    datasets.append(data)
    param_values.append(value)

print(f"Found {len(datasets)} datasets.")

# Sort datasets
order = np.argsort(param_values)
param_values = np.array(param_values)[order]
datasets = [datasets[i] for i in order]





from scipy.interpolate import UnivariateSpline
from scipy.optimize import brentq
import numpy as np

def find_extremum(t, y, i_min, window=5, kind='min'):
    """
    Find a continuous extremum near discrete index i_min.
  
    Parameters
    ----------
    t, y     : full time and value arrays from RK4
    i_min    : index of the discrete extremum (from argmin/argmax)
    window   : number of points on each side to fit the spline on
    kind     : 'min' or 'max'
    """
    # Extract local window — enough points for a good spline, not too many
    lo = max(0, i_min - window)
    hi = min(len(t) - 1, i_min + window)
    t_loc = t[lo:hi+1]
    y_loc = y[lo:hi+1]

# Ensure strictly increasing t
    t_loc, unique_idx = np.unique(t_loc, return_index=True)
    y_loc = y_loc[unique_idx]

# Adjust spline degree safely
    k = min(5, len(t_loc) - 1)

    spl = UnivariateSpline(t_loc, y_loc, k=k, s=0)
    dspl = spl.derivative()

    # Find root of derivative in the interval
    # brentq requires a sign change — check it exists
    da = dspl(t_loc[0])
    db = dspl(t_loc[-1])
    if da * db > 0:
        # No sign change found — fall back to denser window or return discrete min
        #raise ValueError(f"No sign change in derivative over window [{t_loc[0]}, {t_loc[-1]}]")
        t_extremum = 0.0
        y_extremum = 0.0
    else :
        t_extremum = brentq(dspl, t_loc[0], t_loc[-1])
        y_extremum = spl(t_extremum)
    return float(t_extremum), float(y_extremum)


# ============================================================
# Plot hmin for each theta
# ============================================================


hmins = []

n = len(datasets)
ncols = min(3, n)
nrows = math.ceil(n / ncols)

fig, axarr = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows))
axes = np.array(axarr).reshape(-1)

for j in range(n, len(axes)):
    fig.delaxes(axes[j])


for i, (data, pval) in enumerate(zip(datasets, param_values)):
    ax = axes[i]

    print(i, type(data), np.shape(data))


    if( i >= 4):
        t   = data[:-3, 0]
        x_T = data[:-3, 1]
        y_T = data[:-3, 2]
        x_A = data[:-3, 5]
        y_A = data[:-3, 6]
        x_L = data[:-3, 9]
        y_L = data[:-3, 10]
    else:
        t   = data[:, 0]
        x_T = data[:, 1]
        y_T = data[:, 2]
        x_A = data[:, 5]
        y_A = data[:, 6]
        x_L = data[:, 9]
        y_L = data[:, 10]

    r = np.sqrt((x_T - x_A)**2 + (y_T - y_A)**2)
    _, h_min = find_extremum(t, r, np.argmin(r))
    hmins.append(h_min)

    ax.plot(x_T, y_T, color='steelblue', lw=1.5, label='Terre')
    ax.plot(x_L, y_L, color='tomato',    lw=1.5, label='Lune')
    ax.plot(x_A, y_A, color='k',    lw=1.5, label='Artemis')
    ax.plot(x_T[0], y_T[0], 'bo', ms=5)
    ax.plot(x_L[0], y_L[0], 'ro', ms=5)
    ax.plot(x_A[0], y_A[0], 'ko', ms=5)

    ax.ticklabel_format(style='sci', scilimits=(0,0), axis='both')
    ax.set_title(f"{param_name} = {pval:.3g}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend(fontsize=7)

    R_earth = 6378e3
    if h_min < R_earth:
        print("Collision!")

    print(t[-1])

fig.suptitle("Position Terre & Lune & Artemis", fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, "position_terre_lune_3_4_b.png"), dpi=300)
plt.show()
plt.hlines(10000, param_values[0],param_values[-1], colors='k', linestyles='--')
plt.plot(variable_array, hmins, '.b')
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "position_terre_lune_3_4_b.png"), dpi=300)


plt.show()

    #%%
