#%%
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess

# --- 1. CONFIGURATION DES CHEMINS ---
repertoire = '/Users/eldidi/Desktop/Exercise3_student/probleme'
executable = '/engine'
input_filename = 'configuration.in.example'
outdir = os.path.expanduser("~/Desktop/lagrange_D_results")
os.makedirs(outdir, exist_ok=True)
output_path = os.path.join(outdir, "results_L1_instability.txt")

# --- 2. PARAMÈTRES PHYSIQUES ET CALCUL DES POINTS ---
G, m1, m2 = 6.67430e-11, 5.972e24, 7.342e22
x1_val, x2_val = -4668434.6, 379731560.0
d_TL = abs(x2_val - x1_val)
Omega = np.sqrt(G * (m1 + m2) / d_TL**3)
mu = m2 / (m1 + m2)

# Calcul des positions des points colinéaires (Approximations de Hill)
r_L1 = d_TL * (1 - (mu/3)**(1/3))
x_L1 = x1_val + r_L1

r_L2 = d_TL * (1 + (mu/3)**(1/3))
x_L2 = x1_val + r_L2

r_L3 = d_TL * (1 + 5*mu/12)
x_L3 = x1_val - r_L3

print(f"Points de Lagrange trouvés (10^6 m) :")
print(f"L1: {x_L1/1e6:.2f} | L2: {x_L2/1e6:.2f} | L3: {x_L3/1e6:.2f}")

# --- 3. SIMULATION DE L'INSTABILITÉ (SUR L1) ---
# On place Artemis à L1 avec une perturbation de 500km pour voir la chute
perturb = 500e3 
x3_init = x_L1 + perturb
y3_init = perturb

params = {
    'numBodies': 3, 'timeScheme': 1, 'sampling': 50,
    'tEnd': 2592000, # 30 jours suffisent pour voir l'instabilité
    'dt': 0.1, 'tolerance': 1e-7, 'G': G,
    'm1': m1, 'x1': x1_val, 'y1': 0.0, 'vx1': 0.0, 'vy1': Omega * x1_val,
    'm2': m2, 'x2': x2_val, 'y2': 0.0, 'vx2': 0.0, 'vy2': Omega * x2_val,
    'm3': 1000.0, 'x3': x3_init, 'y3': y3_init,
    'vx3': -Omega * y3_init, 'vy3': Omega * x3_init,
}

param_str = " ".join(f"{k}={v:.15g}" for k, v in params.items())
cmd = f"{repertoire}{executable} {input_filename} {param_str} output={output_path}"

print("\nExécution du moteur C++ pour tester L1...")
subprocess.run(cmd, shell=True, capture_output=True, text=True)

# --- 4. CHARGEMENT ET TRANSFORMATION ---
data = np.genfromtxt(output_path, invalid_raise=False)
data = data[~np.isnan(data).any(axis=1)]

t = data[:, 0]
x3, y3 = data[:, 9], data[:, 10]
x1, y1 = data[:, 1], data[:, 2]

# Transformation vers le référentiel tournant R'
x3_rot = x3 * np.cos(Omega * t) + y3 * np.sin(Omega * t)
y3_rot = -x3 * np.sin(Omega * t) + y3 * np.cos(Omega * t)

# --- 5. GÉNÉRATION DES 3 GRAPHES SÉPARÉS ---

# Graphe 1 : Vue Inertielle
plt.figure(figsize=(7, 7))
plt.plot(data[:, 1]/1e6, data[:, 2]/1e6, label='Terre')
plt.plot(data[:, 5]/1e6, data[:, 6]/1e6, label='Lune')
plt.plot(x3/1e6, y3/1e6, label='Artemis (L1 dévié)')


plt.axis('equal'); plt.legend(); plt.grid(True)

# Graphe 2 : Démonstration de l'instabilité (Référentiel Tournant)
plt.figure(figsize=(8, 6))
plt.plot(x3_rot/1e3, y3_rot/1e3, color='red', label='Chute d\'Artemis')
plt.scatter(x_L1/1e3, 0, color='black', marker='x', s=100, label='Point L1 stable théorique')

plt.xlabel("x' (km)"); plt.ylabel("y' (km)")
plt.axis('equal'); plt.legend(); plt.grid(True)

# Graphe 3 : Éloignement par rapport à la Terre
dist_T = np.sqrt((x3-x1)**2 + (y3-y1)**2)
plt.figure(figsize=(10, 4))
plt.plot(t/86400, dist_T/1e3, color='purple')
plt.axhline(y=r_L1/1e3, color='black', linestyle='--', label='Distance L1 théorique')

plt.xlabel("Temps (jours)"); plt.ylabel("Distance (km)")
plt.legend(); plt.grid(True)

plt.show()
# %%
