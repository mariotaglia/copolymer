import sys
import pandas as pd
import matplotlib.pyplot as plt

def cargaDatos(filename):
    # Cargo los datos
    datos = pd.read_csv(filename, sep = '\s+', header = None)
    datos.columns = ['iR', 'iZ', 'fraction']
    return datos

def formatThreeDigit(number):
    return '0'*(3 - len(number)) + number

# Creo la figura
fig, ax = plt.subplots(tight_layout = True)
colors = plt.cm.get_cmap('gnuplot2')
nColors = 4*len(sys.argv[2:]) + 1
colorCounter = 1

# Elige COP/MOL
if sys.argv[1] == '-cop':
    fCdata, fNCdata, fASdata, fIONdata = 'fAmin001.001.', 'fCopANC001.001.', 'fASmol001.001.', 'fCopAion001.001.'
    labelC, labelNC, labelAS, labelION = '$f_{COP}^{C}$, Iteration ', '$f_{COP}^{NC}$, Iteration ', '$f_{AS, COP}^{mol}$, Iteration ', '$f_{COP}^{ion}$, Iteration '
elif sys.argv[1] == '-mol':
    fCdata, fNCdata, fASdata, fIONdata = 'fBHplus001.001.', 'fMolNC001.001.', 'fAScopA001.001.', 'fMolion001.001.'
    labelC, labelNC, labelAS, labelION = '$f_{MOL}^{C}$, Iteration ', '$f_{MOL}^{NC}$, Iteration ', '$f_{AS, MOL}^{cop}$, Iteration ', '$f_{MOL}^{ion}$, Iteration '

for i, iteration in enumerate(sys.argv[2:]):
   # Cargo los datos
   iterCode = formatThreeDigit(iteration)
   if sys.argv[1] == '-cop':
      fC = cargaDatos(fCdata + iterCode + '.dat')
   elif sys.argv[1] == '-mol': 
      fC = cargaDatos(fCdata + iterCode + '.dat')
   fNC = cargaDatos(fNCdata + iterCode + '.dat')
   fAS = cargaDatos(fASdata + iterCode + '.dat')
   fION = cargaDatos(fIONdata + iterCode + '.dat')
   # Grafico
   ax.plot(fC['iR'], fC['fraction'], label = labelC + iteration, c = colors(colorCounter/nColors))
   ax.plot(fNC['iR'], fNC['fraction'], label = labelNC + iteration, c = colors((colorCounter + 1)/nColors))
   ax.plot(fAS['iR'], fAS['fraction'], label = labelAS + iteration, c = colors((colorCounter + 2)/nColors))
   ax.plot(fION['iR'], fION['fraction'], label = labelION + iteration, c = colors((colorCounter + 3)/nColors))
   colorCounter += 4

#Settings del grafico
ax.set_xlabel('$x$ ($nm$)', fontsize = 16)
ax.set_ylabel('Fraction', fontsize = 16)
ax.set_xlim([0.1, 5.0])
ax.legend(fontsize = 12)

plt.show()
