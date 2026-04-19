
#%%
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import subprocess

# ============================================================
# CONFIG
# ============================================================

repertoire = '/Users/eldidi/Desktop/Exercise3_student/probleme'
executable = '/engine'
input_filename = 'configuration.in.example'

outdir = "adaptive_scan_q34a"
os.makedirs(outdir, exist_ok=True)

input_parameters = {
    
    'numBodies': 2,
    'timeScheme': 1,   
    'sampling': 10,
    'tEnd': 172800,

    'dt': 0.5,
    'tolerance': 1e-6,

    'm1': 5.972e24,
    'r1': 6378.1e3,

    'x2': 314159e3,
    'y2': 0,
    'vx2': -1178.501,
    'vy2': 226.13115,
}

# ============================================================
# EPSILON STUDY
# ============================================================

eps_values = np.array([1, 0.1, 0.01, 0.001, 1e-4])

for eps in eps_values:

    params = input_parameters.copy()
    params['tolerance'] = eps

    output_file = os.path.join(outdir, f"run_eps_{eps}.txt")

    param_string = " ".join(f"{k}={v:.15g}" for k, v in params.items())

    cmd = f"{repertoire}{executable} {input_filename} {param_string} output={output_file}"

    subprocess.run(cmd, shell=True)

# ============================================================
# LOAD DATA
# ============================================================

files = sorted(glob.glob(os.path.join(outdir, "*.txt")))

datasets = []
eps_list = []

for f in files:

    data = np.loadtxt(f)

    name = os.path.basename(f)
    eps = float(name.split("_")[-1].replace(".txt", ""))

    datasets.append(data)
    eps_list.append(eps)

order = np.argsort(eps_list)
eps_list = np.array(eps_list)[order]
datasets = [datasets[i] for i in order]

# ============================================================
# ANALYSIS h_min / v_max
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

    d = np.sqrt(x_A**2 + y_A**2)
    h = d - R_T
    v = np.sqrt(vx_A**2 + vy_A**2)

    # h_min
    i = np.argmin(h)
    if 0 < i < len(h) - 1:
        y1, y2, y3 = h[i-1], h[i], h[i+1]
        h_min = y2 - (y3 - y1)**2 / (8*(y3 - 2*y2 + y1))
    else:
        h_min = h[i]

    # v_max
    j = np.argmax(v)
    if 0 < j < len(v) - 1:
        v1, v2, v3 = v[j-1], v[j], v[j+1]
        v_max = v2 - (v3 - v1)**2 / (8*(v3 - 2*v2 + v1))
    else:
        v_max = v[j]

    h_min_list.append(h_min)
    v_max_list.append(v_max)
    N_steps.append(len(t))

h_min = np.array(h_min_list)
v_max = np.array(v_max_list)
eps_list = np.array(eps_list)

# ============================================================
# ERRORS
# ============================================================

h_ref = np.min(h_min)
v_ref = np.min(v_max)

err_h = np.abs(h_min - h_ref) / h_ref
err_v = np.abs(v_max - v_ref) / v_ref

# ============================================================
# plot 1: conv pas adaptatif
# ============================================================

plt.figure(figsize=(12,5))

plt.subplot(1,2,1)
plt.loglog(eps_list, err_h, 'bo-')
plt.xlabel("epsilon")
plt.ylabel("erreur relative")
plt.title("Convergence h_min")
plt.grid(True, which="both")

plt.subplot(1,2,2)
plt.loglog(eps_list, err_v, 'ro-')
plt.xlabel("ε")
plt.ylabel("erreur relative")
plt.title("Convergence v_max")
plt.grid(True, which="both")

plt.tight_layout()
plt.show()

# ============================================================
# plot 2: coût
# ============================================================

plt.figure()
plt.loglog(eps_list, N_steps, 'ko-')
plt.xlabel("ε")
plt.ylabel("nombre de pas")
plt.title("Coût computationnel adaptatif")
plt.grid(True, which="both")
plt.show()

# ============================================================
# plot 3: évolution de dt au cours du temps
# ============================================================

data = datasets[len(datasets)//2]

t = data[:,0]
dt_inst = np.diff(t)

plt.figure()
plt.plot(t[1:], dt_inst)
plt.xlabel("t")
plt.ylabel("Δt")
plt.title("Évolution du pas de temps adaptatif")
plt.grid()
plt.show()

# ============================================================
# plot 4: distance Terre-Artemis
# ============================================================

x_T, y_T = data[:,1], data[:,2]
x_A, y_A = data[:,5], data[:,6]

dist = np.sqrt((x_A-x_T)**2 + (y_A-y_T)**2)

plt.figure()
plt.plot(t, dist)
plt.xlabel("t")
plt.ylabel("distance Terre-Artemis")
plt.title("Distance orbitale")
plt.grid()
plt.show()

# ============================================================
# fixed vs adaptive 
# ============================================================

dt_fixed = input_parameters['dt']
tEnd = input_parameters['tEnd']

N_fixed = int(tEnd / dt_fixed)
N_adapt = np.mean(N_steps)

print("\n===== COMPARISON =====")
print("Fixed     :", N_fixed)
print("Adaptive  :", N_adapt)

