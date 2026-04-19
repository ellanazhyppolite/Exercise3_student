#%%
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import subprocess

# ============================================================
# CONFIGURATION ET SIMULATION
# ============================================================

repertoire = '/Users/eldidi/Desktop/Exercise3_student/probleme'
executable = '/engine'
input_filename = 'configuration.in.example'

outdir = "adaptive_scan_q34a"
os.makedirs(outdir, exist_ok=True)

input_parameters = {
    'numBodies': 2,
    'timeScheme': 1,   
    'sampling': 1,      # Mis à 1 pour une précision maximale sur h_min
    'tEnd': 172800,
    'dt': 0.5,
    'tolerance': 1e-6,
    'G': 6.6740e-11,
    'm1': 5.972e24,
    'r1': 6378.1e3,
    'x2': 314159e3,
    'y2': 0,
    'vx2': -1178.501,
    'vy2': 226.13115,
  
}
# Valeurs d'epsilon à tester
eps_values = np.array([1.0, 0.1, 0.01, 0.001, 1e-4, 1e-5])

for eps in eps_values:
    params = input_parameters.copy()
    params['tolerance'] = eps
    output_file = os.path.join(outdir, f"run_eps_{eps}.txt")
    
    param_string = " ".join(f"{k}={v:.15g}" for k, v in params.items())
    cmd = f"{repertoire}{executable} {input_filename} {param_string} output={output_file}"
    
    print(f"Running simulation for eps = {eps}...")
    subprocess.run(cmd, shell=True)

# ============================================================
# CHARGEMENT ET TRI (CRUCIAL)
# ============================================================

files = glob.glob(os.path.join(outdir, "run_eps_*.txt"))
datasets = []
eps_list = []

for f in files:
    eps = float(os.path.basename(f).replace("run_eps_", "").replace(".txt", ""))
    eps_list.append(eps)
    datasets.append(np.loadtxt(f))

# Conversion en arrays
eps_list = np.array(eps_list)

# Tri : eps_list[0] sera le plus petit (ex: 1e-5), eps_list[-1] le plus grand (1.0)
order = np.argsort(eps_list)
eps_list = eps_list[order]
datasets = [datasets[i] for i in order]

print("Epsilons triés (du plus précis au moins précis) :", eps_list)

# ============================================================
# ANALYSE PHYSIQUE (h_min, v_max, N_steps)
# ============================================================

R_T = input_parameters['r1']
h_min_list = []
v_max_list = []
N_steps = []

for data in datasets:
    t = data[:,0]
    x_T, y_T = data[:,1], data[:,2]
    x_A, y_A = data[:,5], data[:,6]
    vx_A, vy_A = data[:,7], data[:,8]

    dist = np.sqrt((x_A-x_T)**2 + (y_A-y_T)**2)
    h = dist - R_T
    v = np.sqrt(vx_A**2 + vy_A**2)

    # Interpolation quadratique pour h_min
    idx_h = np.argmin(h)
    if 0 < idx_h < len(h) - 1:
        y1, y2, y3 = h[idx_h-1], h[idx_h], h[idx_h+1]
        h_min_list.append(y2 - (y3 - y1)**2 / (8*(y3 - 2*y2 + y1)))
    else:
        h_min_list.append(h[idx_h])

    # Interpolation quadratique pour v_max
    idx_v = np.argmax(v)
    if 0 < idx_v < len(v) - 1:
        v1, v2, v3 = v[idx_v-1], v[idx_v], v[idx_v+1]
        v_max_list.append(v2 - (v3 - v1)**2 / (8*(v3 - 2*v2 + v1)))
    else:
        v_max_list.append(v[idx_v])
        
    N_steps.append(len(t))

h_min = np.array(h_min_list)
v_max = np.array(v_max_list)

# ============================================================
# CALCUL DES ERREURS (Référence = run le plus précis)
# ============================================================

# La référence est le premier élément car on a trié par eps croissant
h_ref = h_min[0]
v_ref = v_max[0]

# Erreur relative par rapport au run le plus précis
err_h = np.abs(h_min[1:] - h_ref) / np.abs(h_ref)
err_v = np.abs(v_max[1:] - v_ref) / np.abs(v_ref)
eps_plot = eps_list[1:]

# ============================================================
# GRAPHIQUES
# ============================================================

# 1. Étude de Convergence
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.loglog(eps_plot, err_h, 'bo-', lw=1.5, label="Erreur sur h_min")
plt.xlabel(r"Tolérance $\epsilon$")
plt.ylabel("Erreur relative")
plt.title("Convergence de l'altitude minimale")
plt.grid(True, which="both", ls="--", alpha=0.5)

plt.subplot(1, 2, 2)
plt.loglog(eps_plot, err_v, 'ro-', lw=1.5, label="Erreur sur v_max")
plt.xlabel(r"Tolérance $\epsilon$")
plt.ylabel("Erreur relative")
plt.title("Convergence de la vitesse maximale")
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout()
plt.show()

# 2. Corrélation Delta t vs Distance (Axe Double)
# On utilise le run le plus précis [0]
data_best = datasets[0]
t_b = data_best[:, 0]
dist_b = np.sqrt((data_best[:,5]-data_best[:,1])**2 + (data_best[:,6]-data_best[:,2])**2)
dt_b = np.diff(t_b)

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.set_xlabel('Temps [s]')
ax1.set_ylabel('Pas de temps Δt [s]', color='tab:blue')
ax1.plot(t_b[1:], dt_b, color='tab:blue', label='Pas de temps Δt')
ax1.set_yscale('log')
ax1.tick_params(axis='y', labelcolor='tab:blue')

ax2 = ax1.twinx()
ax2.set_ylabel('Distance Terre-Artemis [m]', color='tab:red')
ax2.plot(t_b, dist_b, color='tab:red', ls='--', alpha=0.7, label='Distance')
ax2.tick_params(axis='y', labelcolor='tab:red')

plt.title("Lien entre proximité orbitale et adaptation du pas de temps")
plt.show()

# ============================================================
# COMPARAISON ET GAIN
# ============================================================

N_fixed = int(input_parameters['tEnd'] / input_parameters['dt'])

print("\n" + "="*50)
print(f"{'Epsilon':>10} | {'Nb Pas':>10} | {'Gain vs Fixe':>12} | {'h_min [km]':>12}")
print("-" * 50)
for i in range(len(eps_list)):
    gain = N_fixed / N_steps[i]
    print(f"{eps_list[i]:10.1e} | {N_steps[i]:10d} | {gain:11.1f}x | {h_min[i]/1e3:12.2f}")
print("="*50)