#%%

import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess

# =============================================================================
# 1. CONFIGURATION DES CHEMINS ET PARAMÈTRES
# =============================================================================
# Change ces chemins selon ton installation
repertoire = '/Users/eldidi/Desktop/Exercise3_student/probleme'
executable = '/engine'
input_filename = 'configuration.in.example'
outdir = os.path.expanduser("~/Desktop/lagrange_simulation_results")
os.makedirs(outdir, exist_ok=True)
output_path = os.path.join(outdir, "results.txt")

# Paramètres physiques
G  = 6.67430e-11
m1 = 5.972e24    # Terre
m2 = 7.342e22    # Lune
x1_val = -4668434.6 
x2_val = 379731560.0

# Calcul dynamique
d_TL = abs(x2_val - x1_val)
Omega = np.sqrt(G * (m1 + m2) / d_TL**3)

# Position L4 théorique
x3_theo = (d_TL / 2) * (m1 - m2) / (m1 + m2)
y3_theo = (np.sqrt(3) / 2) * d_TL

# --- QUESTION (c) : On place Artemis à 1000km de L4 ---
perturbation = 1000e3 # 1000 km
x3_init = x3_theo + perturbation
y3_init = y3_theo + perturbation

# =============================================================================
# 2. PRÉPARATION ET EXÉCUTION DU MOTEUR
# =============================================================================
params = {
    'numBodies': 3, 'timeScheme': 1, 'sampling': 100,
    'tEnd': 31557600, 'dt': 0.1, 'tolerance': 1e-7, 'G': G,
    'm1': m1, 'x1': x1_val, 'y1': 0.0, 'vx1': 0.0, 'vy1': Omega * x1_val,
    'm2': m2, 'x2': x2_val, 'y2': 0.0, 'vx2': 0.0, 'vy2': Omega * x2_val,
    'm3': 1000.0, 'x3': x3_init, 'y3': y3_init,
    'vx3': -Omega * y3_init, 'vy3': Omega * x3_init,
}

param_str = " ".join(f"{k}={v:.15g}" for k, v in params.items())
cmd = f"{repertoire}{executable} {input_filename} {param_str} output={output_path}"

print("Simulation en cours...")
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

# =============================================================================
# 3. CHARGEMENT ET TRAITEMENT DES DONNÉES
# =============================================================================
if result.returncode == 0:
    # Chargement robuste (ignore les lignes corrompues ou NaNs)
    data = np.genfromtxt(output_path, invalid_raise=False)
    data = data[~np.isnan(data).any(axis=1)]

    t = data[:, 0]
    t_jours = t / 86400
    # Positions Inertielles (x, y) pour les 3 corps
    pos = {
        'Terre':  (data[:, 1], data[:, 2]),
        'Lune':   (data[:, 5], data[:, 6]),
        'Artemis':(data[:, 9], data[:, 10])
    }

    # Calcul des distances
    d_AT = np.sqrt((pos['Artemis'][0] - pos['Terre'][0])**2 + (pos['Artemis'][1] - pos['Terre'][1])**2)
    d_AL = np.sqrt((pos['Artemis'][0] - pos['Lune'][0])**2 + (pos['Artemis'][1] - pos['Lune'][1])**2)
    ecart_relatif = np.abs(d_AT - d_AL) / d_TL * 100

    # Transformation vers le Référentiel Tournant R' (Question c)
    def to_rot(x, y, t):
        xr = x * np.cos(Omega * t) + y * np.sin(Omega * t)
        yr = -x * np.sin(Omega * t) + y * np.cos(Omega * t)
        return xr, yr

    x3_r, y3_r = to_rot(pos['Artemis'][0], pos['Artemis'][1], t)

    # =========================================================================
    # 4. GÉNÉRATION DES GRAPHES (SÉPARÉS)
    # =========================================================================
    
    # --- GRAPHE 1 : Trajectoires Inertielles ---
    plt.figure(figsize=(7, 7))
    plt.plot(pos['Terre'][0]/1e6, pos['Terre'][1]/1e6, label='Terre')
    plt.plot(pos['Lune'][0]/1e6, pos['Lune'][1]/1e6, label='Lune')
    plt.plot(pos['Artemis'][0]/1e6, pos['Artemis'][1]/1e6, label='Artemis')
    plt.title("Trajectoires dans le Référentiel Inertiel")
    plt.xlabel("x (10^6 m)"); plt.ylabel("y (10^6 m)")
    plt.axis('equal'); plt.legend(); plt.grid(True)

    # --- GRAPHE 2 : Référentiel Tournant R' (La Stabilité) ---
    plt.figure(figsize=(7, 7))
    plt.plot(x3_r/1e3, y3_r/1e3, color='green', label='Artemis (L4 + 1000km)')
    plt.scatter(x3_theo/1e3, y3_theo/1e3, color='red', marker='x', label='Point L4 fixe')
    plt.title("Mouvement dans le Référentiel Tournant Terre-Lune")
    plt.xlabel("x' (km)"); plt.ylabel("y' (km)")
    plt.axis('equal'); plt.legend(); plt.grid(True)

    # --- GRAPHE 3 : Écart Relatif ---
    plt.figure(figsize=(10, 4))
    plt.plot(t_jours, ecart_relatif, color='orange')
    plt.title(f"Écart relatif de distance (Libration) - Moyenne: {np.mean(ecart_relatif):.4f}%")
    plt.xlabel("Temps (jours)"); plt.ylabel("Écart (%)")
    plt.grid(True)

    plt.show()
    print("Graphiques générés avec succès.")
else:
    print("Erreur moteur:", result.stderr)