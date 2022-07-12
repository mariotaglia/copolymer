import sys
import pandas as pd
import matplotlib.pyplot as plt

def cargaDatos(filename):
    # Cargo los datos
    datos = pd.read_csv(filename, sep = '\s+', header = None)
    datos.columns = ['iR', 'iZ', 'density']
    return datos

def formatThreeDigit(number):
    return '0'*(3 - len(number)) + number

# Creo la figura
fig, ax = plt.subplots(tight_layout = True)
colors = plt.cm.get_cmap('gnuplot2')
nColors = 4*len(sys.argv[1:]) + 1
colorCounter = 1

for i, iteration in enumerate(sys.argv[1:]):
   # Cargo los datos
   iterCode = formatThreeDigit(iteration)
   densityPol1_1 = cargaDatos('densitypolymer001.001.001.' + iterCode + '.dat')
   densityPol1_2 = cargaDatos('densitypolymer001.002.001.' + iterCode + '.dat')
   densityPol2_2 = cargaDatos('densitypolymer002.002.001.' + iterCode + '.dat')
   densitySolv = cargaDatos('densitysolvent.001.' + iterCode + '.dat')
   # Grafico
   ax.plot(densityPol1_1['iR'], densityPol1_1['density'], label = 'Pol 1, type 1, Iteration ' + iteration, c = colors(colorCounter/nColors))
   ax.plot(densityPol1_2['iR'], densityPol1_2['density'], label = 'Pol 1, type 2, Iteration ' + iteration, c = colors((colorCounter + 1)/nColors))
   ax.plot(densityPol2_2['iR'], densityPol2_2['density'], label = 'Pol 2, type 2, Iteration ' + iteration, c = colors((colorCounter + 2)/nColors))
   ax.plot(densitySolv['iR'], densitySolv['density'], label = 'Solvent, Iteration ' + iteration, c = colors((colorCounter + 3)/nColors))
   colorCounter += 3

#Settings del grafico
ax.set_title('density pol', fontsize = 16)
ax.set_xlabel('$x$ ($nm$)', fontsize = 16)
ax.set_ylabel('Density', fontsize = 16)
ax.set_xlim([0.1, 6.0])
ax.legend(fontsize = 8)

plt.show()
