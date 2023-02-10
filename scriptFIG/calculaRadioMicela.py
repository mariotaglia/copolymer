import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def cargaDatos(filename):
    # Cargo los datos
    datos = pd.read_csv(filename, sep = '\s+', header = None)
    datos.columns = ['iR', 'iZ', 'density']
    return datos

def formatThreeDigit(number):
    return '0'*(3 - len(number)) + number

def calculaRadioMicela(r, rho):
    rhoLiq = np.max(rho)
    rhoDeriv = np.gradient(rho, r)
    integral = np.trapz(r**3*rhoDeriv, r)
    radioMic = np.power((-1.0/rhoLiq)*integral, 1/3)
    return radioMic

# Creo la figura
fig, ax = plt.subplots(figsize = (9, 6), tight_layout = True)
colors = plt.cm.get_cmap('gnuplot2')
nColors = 3*len(sys.argv[1:]) + 1
colorCounter = 1

for i, iteration in enumerate(sys.argv[1:]):
   # Cargo los datos
   iterCode = formatThreeDigit(iteration)
   densityPol1_1 = cargaDatos('densitypolymer001.001.001.' + iterCode + '.dat')
   densityPol1_2 = cargaDatos('densitypolymer001.002.001.' + iterCode + '.dat')
   densityPol2_2 = cargaDatos('densitypolymer002.002.001.' + iterCode + '.dat')
   # Calculo del radio de la micela
   sumaDensidades = densityPol1_1['density'] + densityPol1_2['density'] + densityPol2_2['density']
   sumaDensidadesCargados = densityPol1_2['density'] + densityPol2_2['density']
   indexMaxNeutros = densityPol1_1['density'].argmax()
   radiosNeutros = densityPol1_1['iR'].iloc[indexMaxNeutros:]
   sumaDensidadesNeutros = densityPol1_1['density'].iloc[indexMaxNeutros:]
   radioMic = calculaRadioMicela(densityPol2_2['iR'], sumaDensidades)
   radioMicCargados = calculaRadioMicela(densityPol2_2['iR'], sumaDensidadesCargados)
   radioMicNeutros = calculaRadioMicela(radiosNeutros, sumaDensidadesNeutros)
   print('El radio de la micela de la iteracion', iteration, 'es', radioMic)
   print('El radio de la micela (CARGADOS) de la iteracion', iteration, 'es', radioMicCargados)
   print('El radio de la micela (NEUTROS) de la iteracion', iteration, 'es', radioMicNeutros)
   # Grafico
   ax.axvline(x = radioMic, c= 'k', lw = 3, ls = 'dashed')
   ax.plot(densityPol2_2['iR'], sumaDensidades, label = r'$\rho(r)$, Iteration ' + iteration, c = 'r', lw = 4)
   ax.plot(densityPol1_1['iR'], densityPol1_1['density'], label = 'Cop (Neutral), Iteration ' + iteration, c = colors(colorCounter/nColors), lw = 4)
   ax.plot(densityPol1_2['iR'], densityPol1_2['density'], label = 'Cop (Acid), Iteration ' + iteration, c = colors((colorCounter + 1)/nColors), lw = 4)
   ax.plot(densityPol2_2['iR'], densityPol2_2['density'], label = 'Mol (Basic), Iteration ' + iteration, c = colors((colorCounter + 2)/nColors), lw = 4)
   colorCounter += 2
   ax.axvline(x = radioMic, c= 'k', lw = 3, ls = 'dashed')
   ax.axvline(x = radioMicCargados, c= 'r', lw = 3, ls = 'dashed')
   ax.axvline(x = radioMicNeutros, c= 'b', lw = 3, ls = 'dashed')

#Settings del grafico
ax.set_xlabel('$r$ ($nm$)', fontsize = 30)
ax.set_ylabel('Density polymer', fontsize = 30)
ax.set_xlim([0.1, 140*0.2])
#ax.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0))
#ax.yaxis.get_offset_text().set_fontsize(20)
ax.tick_params(axis = 'both', labelsize = 20)
ax.legend(fontsize = 6)

plt.show()
