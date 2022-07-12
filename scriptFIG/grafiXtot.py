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
nColors = 2*len(sys.argv[1:]) + 1
colorCounter = 1

for i, iteration in enumerate(sys.argv[1:]):
   # Cargo los datos
   iterCode = formatThreeDigit(iteration)
   xTotPol1 = cargaDatos('xdensitytota.001.001.' + iterCode + '.dat')
   xTotPol2 = cargaDatos('xdensitytota.002.001.' + iterCode + '.dat')
   # Grafico
   ax.plot(xTotPol1['iR'], xTotPol1['density'], label = 'Pol 1, Iteration ' + iteration, c = colors(colorCounter/nColors))
   ax.plot(xTotPol2['iR'], xTotPol2['density'], label = 'Pol 2, Iteration ' + iteration, c = colors((colorCounter + 1)/nColors))
   colorCounter += 1

#Settings del grafico
ax.set_title('x total', fontsize = 16)
ax.set_xlabel('$x$ ($nm$)', fontsize = 16)
ax.set_ylabel('Density polymer', fontsize = 16)
ax.set_xlim([0.1, 6.0])
ax.legend(fontsize = 8)

plt.show()
