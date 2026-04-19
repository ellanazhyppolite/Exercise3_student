#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os
import glob
import re
import math
import subprocess
import matplotlib.colors as colors



repertoire = '/Users/eldidi/Desktop/Exercise3_student/probleme'
executable = '/engine'
input_filename = 'configuration.in.example'

input_parameters = {
    'numBodies': 2,
    'timeScheme': 1, # 1 pour adaptatif
    'sampling': 0, 
    'tEnd': 172800,  
    'dt': 0.5,
    'tolerance': 1e-6,
    'G': 6.6740e-11,

    # Atmosphere
    'useAtmosphere' : 1, 
    'roh0' : 1.2,
    'atmosphereScale': 7238.2,
    'dragArea': 19.7923475393,       
    'dragCoefficient': 0.3,
    'dragBody': 1,  
    'dragCenterBody': 0, 

    # Terre
    'm1': 5.972e24,
    'r1': 6378.1e3, 
    'x1': 0.0, 'y1': 0, 'vx1': 0.0, 'vy1': 0.0,

    # Artémis 
    'm2': 8500,
    'r2': 2.26, 
    'x2': 314159e3, 
    'y2': 0,  
    'vx2': -1178.501,  
    #'vy2': 226.13115,
    'vy2': 235.0,  
}


paramstr = 'tolerance' 
variable_array = [1, 0.1, 0.01, 0.001, 1e-4, 1e-5]

outstr = f"syst_gravitationnel_q34a"
outdir = outstr
os.makedirs(outdir, exist_ok=True)



for val in variable_array:
    params = input_parameters.copy()
    params[paramstr] = val
    output_path = os.path.join(outdir, f"run_{paramstr}_{val}.txt")
    param_string = " ".join(f"{k}={v:.15g}" for k, v in params.items())
    cmd = f"{repertoire}{executable} {input_filename} {param_string} output={output_path}"
    print(f"Running: {paramstr}={val}")
    subprocess.run(cmd, shell=True)



files = sorted(glob.glob(os.path.join(outdir, "*.txt")))
datasets = []
param_values = []

for f in files:
    val = float(re.findall(r"[-+]?\d*\.\d+|\d+", f.split('_')[-1])[0])
    datasets.append(np.loadtxt(f))
    param_values.append(val)


order = np.argsort(param_values)
param_values = np.array(param_values)[order]
datasets = [datasets[i] for i in order]


a_max_list = []
p_max_list = []
h_min_list = []

print("\n--- ANALYSE DES MAXIMA (Run le plus précis) ---")

for i, data in enumerate(datasets):
    t = data[:, 0]

    x, y = data[:, 5], data[:, 6]
    vx, vy = data[:, 7], data[:, 8]
    

    dist = np.sqrt(x**2 + y**2)
    h = dist - input_parameters['r1']
    h_min_list.append(np.min(h))
    

    dvx = np.diff(vx) / np.diff(t)
    dvy = np.diff(vy) / np.diff(t)
    accel = np.sqrt(dvx**2 + dvy**2)
    a_max_list.append(np.max(accel))
    

    v_mag = np.sqrt(vx**2 + vy**2)
    rho = input_parameters['roh0'] * np.exp(-h / input_parameters['atmosphereScale'])
    f_drag = 0.5 * rho * v_mag**2 * input_parameters['dragArea'] * input_parameters['dragCoefficient']
    power = f_drag * v_mag
    p_max_list.append(np.max(power))

print(f"Accélération max : {a_max_list[0]/9.81:.2f} g")
print(f"Puissance traînée max : {p_max_list[0]/1e6:.2f} MW")
print(f"Altitude min : {h_min_list[0]:.2f} m")

# ============================================================
# plot 1: trajectoires 
# ============================================================

fig_all, ax_all = plt.subplots(figsize=(12, 8))
fig_all.subplots_adjust(left=0.15, right=0.85) 

norm = colors.LogNorm(vmin=min(param_values), vmax=max(param_values))
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)

for i, (data, pval) in enumerate(zip(datasets, param_values)):
    ax_all.plot(data[:, 5], data[:, 6], color=sm.to_rgba(pval), lw=1, alpha=0.7)

terre = plt.Circle((0, 0), input_parameters['r1'], color='blue', alpha=0.2, label='Terre')
ax_all.add_artist(terre)
ax_all.set_aspect('equal')
ax_all.set_xlabel("x (m)")
ax_all.set_ylabel("y (m)")
ax_all.set_title(f"Trajectoires vs {paramstr}")
ax_all.ticklabel_format(style='sci', scilimits=(0,0), axis='both')
cbar = fig_all.colorbar(sm, ax=ax_all, fraction=0.046, pad=0.04)
cbar.set_label(f'Valeur de {paramstr}')
plt.show()

# ============================================================
# plot 2 : conv
# ============================================================

# Calcul de l'erreur relative par rapport au run le plus précis (index 0)
ref_a = a_max_list[0]
ref_p = p_max_list[0]

err_a = [abs(a - ref_a)/ref_a for a in a_max_list[1:]]
err_p = [abs(p - ref_p)/ref_p for p in p_max_list[1:]]

plt.figure(figsize=(10, 5))
plt.loglog(param_values[1:], err_a, 'bo-', label='Erreur Accélération Max')
plt.loglog(param_values[1:], err_p, 'ro-', label='Erreur Puissance Max')
plt.xlabel("Tolérance (epsilon)")
plt.ylabel("Erreur Relative")
plt.title("Étude de convergence (Pas Adaptatif)")
plt.grid(True, which="both", ls="--")
plt.legend()
plt.show()