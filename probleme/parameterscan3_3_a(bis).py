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
from scipy.interpolate import UnivariateSpline
from scipy.optimize import brentq



repertoire = '/Users/eldidi/Desktop/Exercise3_student/probleme'
executable = '/engine'
input_filename = 'configuration.in.example'

input_parameters = {
    'numBodies': 2,
    'timeScheme': 1, # 1 pour adaptatif
    'sampling': 0, 
    'tEnd': 172800,  
    'dt': 0.1,

    'tolerance': 1e-6,
    'G': 6.6740e-11,

    # Atmosphere
    'useAtmosphere' : 1, 
    'rho0' : 1.2,
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
    'vx2': -1178.511308,  
    'vy2': 226.0776347,
   

    
}


paramstr = 'tolerance' 
#variable_array = [1, 0.1, 0.01, 0.001, 1e-4, 1e-5]

variable_array = [1e-5]


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

    # Fit spline (k=4 gives a differentiable cubic-like derivative, k=5 is smoother)
    # s=0 forces interpolation through all points (appropriate for RK4 output)
    spl = UnivariateSpline(t_loc, y_loc, k=5, s=0)
    dspl = spl.derivative()

    # Find root of derivative in the interval
    # brentq requires a sign change — check it exists
    da = dspl(t_loc[0])
    db = dspl(t_loc[-1])
    if da * db > 0:
        # No sign change found — fall back to denser window or return discrete min
        raise ValueError(f"No sign change  derivative over window [{t_loc[0]}, {t_loc[-1]}]")

    t_extremum = brentq(dspl, t_loc[0], t_loc[-1])
    y_extremum = spl(t_extremum)
    return float(t_extremum), float(y_extremum)







data = datasets[0]
t = data[:, 0]

x, y = data[:, 7], data[:, 8]
vx, vy = data[:, 9], data[:, 10]
a_x, a_y = data[:, 11], data[:, 12]
    
a = np.sqrt(a_x**2 +a_y**2)/9.81 # pour avoir un multiplwes de g

h = np.sqrt(x**2 + y**2) - input_parameters['r1']


v = np.sqrt(vx**2 + vy**2)
rho_p = input_parameters['rho0'] * np.exp(-h / input_parameters['atmosphereScale'])
f_drag_p = 0.5 * rho_p * v**2 * input_parameters['dragArea'] * input_parameters['dragCoefficient']
power_p = (f_drag_p * v) / 1e6  # Conversion en MW


t_max_a, a_max = find_extremum(t,a,np.argmax(a),window = 5,kind='max')

t_max_p, p_max = find_extremum(t,power_p,np.argmax(power_p),window = 5,kind='max')

print(t_max_a,a_max)
print(t_max_p,p_max)



fig, (ax1, ax2,ax3) = plt.subplots(3, 1, figsize=(10, 10), sharex=True)


# --- Graphique 1 : Accélération ---
ax1.plot(t, a, color='tab:red', label='Profil Accélération')

ax1.scatter(t_max_a, a_max, color='black', zorder=5, label=f'Max: {a_max:.2f} g')
ax1.set_ylabel("Accélération (g)")
ax1.set_title("Analyse de la rentrée atmosphérique")
ax1.grid(True, alpha=0.3)
ax1.legend()

# --- Graphique 2 : Puissance ---
ax2.plot(t, power_p, color='tab:orange', label='Profil Puissance')

ax2.scatter(t_max_p, p_max, color='black', zorder=5, label=f'Max: {p_max:.2f} MW')
ax2.set_ylabel("Puissance dissipée (MW)")
ax2.set_xlabel("Temps (s)")
ax2.grid(True, alpha=0.3)
ax2.legend()

# Zoom sur la phase active de rentrée
ax2.set_xlim(t_max_a - 3000, t_max_a + 3000) 


ax3.plot(x,y)
ax3.grid(True,alpha=0.3)

plt.tight_layout()
plt.show()
    
    # %%
