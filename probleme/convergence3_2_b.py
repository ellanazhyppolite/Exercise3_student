
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
#repertoire = '/Users/ella/Documents/EPFL/BA4/PHYSNUM/physnum_ex3/probleme'
repertoire = '/Users/eldidi/Desktop/Exercise3_student/probleme'
executable = '/engine'
input_filename = 'configuration.in.example' # Strictly no longer needed, but we keep it for now to avoid having to change the code in engine.cpp


input_parameters = {
    'numBodies': 2,
    'timeScheme': 0, # 0 pour dt fixe, 1 pour schema adaptatif
    'sampling': 100, # mettre 0 pour ecrire tout les pas
    'tEnd': 172800,  # 2 jours


   

    'dt': 0.5,
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
    'x1': 0.0,
    'y1': 0,  # mètres
    'vx1': 0.0,
    'vy1': 0.0,

    # Artémis 
    'm2': 8500,
    'r2': 2.26, # rayon de la sonde
    'x2': 314159e3, # en mètres
    'y2': 0,  # mètres
    'vx2': -1178.0,  # v_circulaire = sqrt(GM/r) ≈ 1018 m/s
    'vy2': 226.0,  #en m/s


    #on peut ajouter d'autre corps de la même manière !!!!il faut qu'il y en ai autant que numbodies!!!!
}


# -------------------------------------------------



# Updated from last time, the code below can now be used to scan any parameter, just make sure to update the paramstr and the variable_array accordingly


paramstr = 'dt'
variable_array = 2**np.linspace(-4,0,20)  # CHANGER LA VALEUR SELON QUESTION

# vu que je commprends pas ce qui se passe, mettre la valeur de tEnd ici aussi please

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
'''


print(f"tEnd = {input_parameters['tEnd']:.2e}")
print(f"y2   = {input_parameters['y2']:.2e}")
print(f"vx2  = {input_parameters['vx2']:.2e}")
#T_orbital = 2 * np.pi * input_parameters['y2'] / input_parameters['vx2']

#print(f"T orbital = {T_orbital:.2e} s")

#print(f"tEnd / T  = {input_parameters['tEnd'] / T_orbital:.2f} orbits")

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
    x_A = data[:, 5]
    y_A = data[:, 6]
   

   

    



    ax.plot(x_T, y_T, color='steelblue', lw=1.5, label='Terre')
    ax.plot(x_A, y_A, color='tomato', lw=1.5, label='Artemis')
    ax.plot(x_T[0], y_T[0], 'bo', ms=5)
    ax.plot(x_A[0], y_A[0], 'ro', ms=5)

    ax.ticklabel_format(style='sci', scilimits=(0,0), axis='both')
    ax.set_title(f"{param_name} = {pval:.3g}")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.legend(fontsize=7)

    

fig.suptitle("Position Terre & Artemis", fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, "position_terre_artemis_q34a.png"), dpi=300)

print(np.shape(x_A))
print(np.shape(y_A))

plt.show()
'''
# %%
#==============================================================
# Plot 2 : minimum h et maximum v 
#==============================================================
h_min = []
v_max = []



for data in datasets:
    t = data[:, 0]
    x_A = data[:, 5]
    y_A = data[:, 6]
    vx_A = data[:, 7]
    vy_A = data[:, 8]


    

    R_T = input_parameters['r1']

    d = np.sqrt(x_A**2 + y_A**2)
    
    h =  d - R_T
    v = np.sqrt(vx_A**2 + vy_A**2)

    # --- Interpolation pour h_min ---
    i = np.argmin(h)
    if 0 < i < len(h) - 1:
        y1, y2, y3 = h[i-1], h[i], h[i+1]
      
        h_min_val = y2 - (y3 - y1)**2 / (8 * (y3 - 2*y2 + y1))
    else:
        h_min_val = h[i]

    # --- Interpolation pour v_max ---
    j = np.argmax(v)
    if 0 < j < len(v) - 1:
        v1, v2, v3 = v[j-1], v[j], v[j+1]
        # Formule identique pour le maximum
        v_max_val = v2 - (v3 - v1)**2 / (8 * (v3 - 2*v2 + v1))
    else:
        v_max_val = v[j]

    


    h_min.append(h_min_val)
    v_max.append(v_max_val)
 
    '''
    h_min.append(h.min())  # apres faire avec interpolation sur les 3 points minimaux plutot  
    v_max.append(v.max())
    '''


print("h_min:", h_min)
print("v_max:", v_max)

fig, axs = plt.subplots(2, 1)

axs[0].plot(param_values, h_min, marker='o', color='blue')
axs[0].set_title("h_min")

axs[1].plot(param_values, v_max, marker='s', color='red')
axs[1].set_title("v_max")




h_ref = h_min[0]
v_ref = v_max[0]


dts_plot = param_values[1:]
err_h = [abs(h - h_ref) / h_ref for h in h_min[1:]]
err_v = [abs(v - v_ref) / v_ref for v in v_max[1:]]


plt.figure(figsize=(10, 5))



h_ref = h_min[0]
v_ref = v_max[0]


dts_plot = param_values[1:]
err_h = [abs(h - h_ref) / h_ref for h in h_min[1:]]
err_v = [abs(v - v_ref) / v_ref for v in v_max[1:]]




plt.figure(figsize=(12, 6))


plt.subplot(1, 2, 1)


log_dt = np.log(dts_plot)
log_err_h = np.log(err_h)
p_h, b_h = np.polyfit(log_dt, log_err_h, 1)

plt.loglog(dts_plot, err_h, 'bo', label='Données $h_{min}$')


plt.loglog(dts_plot, np.exp(b_h) * dts_plot**p_h, 'b-', 
           label=f'Fit linéaire (pente = {p_h:.2f})')


pente_4_h = [err_h[0] * (d / dts_plot[0])**4 for d in dts_plot]
plt.loglog(dts_plot, pente_4_h, 'k--', alpha=0.5, label='Théorie (pente 4)')

plt.xlabel('dt [s]')
plt.ylabel('Erreur relative')
plt.title(f'Convergence $h_{min}$ (Ordre exp: {p_h:.2f})')
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.legend()


plt.subplot(1, 2, 2)


log_err_v = np.log(err_v)
p_v, b_v = np.polyfit(log_dt, log_err_v, 1)


plt.loglog(dts_plot, err_v, 'ro', label='Données $v_{max}$')


plt.loglog(dts_plot, np.exp(b_v) * dts_plot**p_v, 'r-', 
           label=f'Fit linéaire (pente = {p_v:.2f})')


pente_4_v = [err_v[0] * (d / dts_plot[0])**4 for d in dts_plot]
plt.loglog(dts_plot, pente_4_v, 'k--', alpha=0.5, label='Théorie (pente 4)')

plt.xlabel('dt [s]')
plt.ylabel('Erreur relative')
plt.title(f'Convergence $v_{max}$ (Ordre exp: {p_v:.2f})')
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.legend()

plt.tight_layout()
plt.show()