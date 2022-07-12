import pandas as pd
import matplotlib.pyplot as plt

def cargaDatos(filename):
    # Cargo los datos
    datos = pd.read_csv(filename, sep = '\s+', header = None)
    datos.columns = ['nPol', 'Ftot']
    return datos

# Cargo los datos
ftot = cargaDatos('F_tot.dat')

# Grafico
fig, ax = plt.subplots(tight_layout = True)
ax.scatter(ftot['nPol'], ftot['Ftot'])
# Settings del grafico
ax.set_xlabel('nPol', fontsize = 16)
ax.set_ylabel('$F_{TOT}/N$ ($k_{B}T$)', fontsize = 16)

plt.show()

print('El minimo en Ftot', min(ftot['Ftot']))
