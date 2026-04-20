#%%
import numpy as np
import matplotlib.pyplot as plt


G = 6.674e-11
mT = 5.972e24      
mL = 7.3477e22    
d = 384748e3       
M = mT + mL


x1 = -(mL / M) * d
x2 = (mT / M) * d
omega = np.sqrt(G * M / d**3)



limit = 1.5 * d
x = np.linspace(-limit, limit, 400)
y = np.linspace(-limit, limit, 400)
X, Y = np.meshgrid(x, y)

#calcul du potntiel effectif = Potentiel gravitationnel + potentiel centrifuge

r1 = np.sqrt((X - x1)**2 + Y**2)
r2 = np.sqrt((X - x2)**2 + Y**2)



eps = 1e-7# Pour éviter division par 0 aux centres jsp si ca change qlqc
V_grav = -G * mT / (r1 + eps) - G * mL / (r2 + eps)
V_cent = -0.5 * (omega**2) * (X**2 + Y**2)
V_total = V_grav + V_cent


plt.figure(figsize=(10, 8))


v_min, v_max = -1.8e6, -1.2e6 
levels = np.linspace(v_min, v_max, 40)

contour = plt.contour(X, Y, V_total, levels=levels, cmap='viridis')
plt.colorbar(contour, label="Potentiel effectif [J/kg]")


plt.plot(x1, 0, 'bo', markersize=8, label='Terre')
plt.plot(x2, 0, 'ro', markersize=4, label='Lune')


plt.plot([d/2 - (mL/M)*d], [np.sqrt(3)/2 * d], 'gx', label='L4/L5 (approx)')
plt.plot([d/2 - (mL/M)*d], [-np.sqrt(3)/2 * d], 'gx')


plt.xlabel("x' [m]")
plt.ylabel("y' [m]")
plt.legend()
plt.gca().set_aspect('equal')
plt.grid(alpha=0.3)


plt.savefig('potentiel_lagrange.pdf')
plt.show()