
#%%

import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess

# --- CONFIGURATION DES CHEMINS ---
repertoire = '/Users/eldidi/Desktop/Exercise3_student/probleme'
executable = '/engine'
input_filename = 'configuration.in.example'
outdir = 'lagrange_L4_simulation'
os.makedirs(outdir, exist_ok=True)
output_path = os.path.join(outdir, "lagrange_L4_results.txt")

# --- 1. PARAMÈTRES DE BASE ---
G  = 6.67430e-11
m1 = 5.972e24    
m2 = 7.342e22    
x1_val = -4668434.6 
x2_val = 379731560.0

# --- 2. CALCUL DYNAMIQUE (Auto-cohérence) ---
d_TL = abs(x2_val - x1_val)
Omega = np.sqrt(G * (m1 + m2) / d_TL**3)

x3_theo = (d_TL / 2) * (m1 - m2) / (m1 + m2)
y3_theo = (np.sqrt(3) / 2) * d_TL

# Correctif pour éviter la division par zéro dans le moteur C++
x3_init = x3_theo + 10.0
y3_init = y3_theo + 10.0

# --- 3. DICTIONNAIRE POUR LE MOTEUR C++ ---
input_parameters = {
    'numBodies': 3,
    'timeScheme': 1,
    'sampling': 100,      
    'tEnd': 31557600,     
    'dt': 0.1,            
    'tolerance': 1e-7,    
    'G': G,
    'useAtmosphere': 0, 

    'm1': m1, 'r1': 6371000, 'x1': x1_val, 'y1': 0.0, 
    'vx1': 0.0, 'vy1': Omega * x1_val,

    'm2': m2, 'r2': 1737000, 'x2': x2_val, 'y2': 0.0,  
    'vx2': 0.0, 'vy2': Omega * x2_val,

    'm3': 1000.0, 'r3': 1.0, 'x3': x3_init, 'y3': y3_init,  
    'vx3': -Omega * y3_init, 'vy3':  Omega * x3_init,  
}

# --- 4. EXÉCUTION ---
param_string = " ".join(f"{k}={v:.15g}" for k, v in input_parameters.items())
cmd = f"{repertoire}{executable} {input_filename} {param_string} output={output_path}"

print(f"Lancement de la simulation...")
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

# --- 5. ANALYSE ET GÉNÉRATION DE 3 GRAPHES SÉPARÉS ---
if result.returncode == 0:
    try:
        data = np.genfromtxt(output_path, invalid_raise=False)
        data = data[~np.isnan(data).any(axis=1)] 
        
        if len(data) < 2:
            raise ValueError("Données insuffisantes.")

        t = data[:, 0]
        t_jours = t / 86400
        x1_p, y1_p = data[:, 1], data[:, 2]
        x2_p, y2_p = data[:, 5], data[:, 6]
        x3_p, y3_p = data[:, 9], data[:, 10]
        
        d_AT = np.sqrt((x3_p - x1_p)**2 + (y3_p - y1_p)**2)
        d_AL = np.sqrt((x3_p - x2_p)**2 + (y3_p - y2_p)**2)
        ecart_relatif = np.abs(d_AT - d_AL) / d_TL * 100

        # --- GRAPHE 1 : Trajectoires ---
        plt.figure(figsize=(8, 8))
        plt.plot(x1_p/1e6, y1_p/1e6, label='Terre', color='royalblue')
        plt.plot(x2_p/1e6, y2_p/1e6, label='Lune', color='gray')
        plt.plot(x3_p/1e6, y3_p/1e6, label='Artemis (L4)', color='limegreen', linewidth=2)

        plt.xlabel("x ($10^6$ m)")
        plt.ylabel("y ($10^6$ m)")
        plt.legend()
        plt.axis('equal')
        plt.grid(alpha=0.3)

        # --- GRAPHE 2 : Évolution des Distances ---
        plt.figure(figsize=(10, 5))
        plt.plot(t_jours, d_AT/1e3, label='Distance Terre-Artemis', color='royalblue')
        plt.plot(t_jours, d_AL/1e3, label='Distance Lune-Artemis', color='crimson', linestyle='--')

        plt.xlabel("Temps (jours)")
        plt.ylabel("Distance (km)")
        plt.legend()
        plt.grid(alpha=0.3)

        # --- GRAPHE 3 : Écart Relatif ---
        plt.figure(figsize=(10, 5))
        plt.plot(t_jours, ecart_relatif, color='darkorange', linewidth=1.5)

        plt.xlabel("Temps (jours)")
        plt.ylabel("Écart relatif (%)")
        plt.grid(True, alpha=0.3)
        
        plt.show()

        print(f"Simulation réussie sur {t_jours[-1]:.1f} jours.")

    except Exception as e:
        print(f"Erreur d'analyse : {e}")
else:
    print(f"Erreur moteur C++ : {result.stderr}")