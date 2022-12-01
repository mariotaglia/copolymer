import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def cargaDatos(filename):
    # Cargo los datos
    datos = pd.read_csv(filename, sep = '\s+', header = None)
    datos.columns = ['iR', 'iZ', 'density']
    datos['r2*density'] = datos['density']*datos['iR']**2#*0.2**2
    return datos

def formatThreeDigit(number):
    return '0'*(3 - len(number)) + number

def calculoR(density):
    der = []
    r = []
    for i in range(0,139):
       bb = (density[i+1]-density[i])/0.2
       der.append(bb)
       bb = 0.2*(i+0.5)
       r.append(bb)
    inte = []
    r2 = []
    sumin = 0 
    for i in range(0,138):
       sumin =sumin +1/(min(density)-max(density))*(r[i+1]**3*der[i+1]+r[i]**3*der[i])/2*0.2
       inte.append(sumin)
       bb = 0.2*(i+0.5)
       r2.append(bb)
    return inte # r,der ,r2,inte
def buscoR(density):
    int1 = 0
    int2 = 0
    for i in range(0,139):
        int1 = int1 + (density[i+1]*((i+1)*0.2)**3+density[i]*((i)*0.2)**3)/2*0.2
        int2 = int2 + (density[i+1]*((i+1)*0.2)**2+density[i]*((i)*0.2)**2)/2*0.2
    return 4.0/3.0*int1/int2

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

   densityTOT  = densityPol1_1['density'] + densityPol1_2['density'] + densityPol2_2['density']
   ax.plot(densityPol2_2['iR'], densityTOT, label = 'Total, Iteration ' + iteration, c = colors((colorCounter )/nColors), lw = 4)

   inte = calculoR(densityTOT)
   #ax.plot(r, der, label = 'derivada, Iteration ' + iteration, c = colors((colorCounter + 1)/nColors), lw = 4)
   #ax.plot(r2, inte, label = 'Integral, Iteration ' + iteration, c = colors((colorCounter + 2)/nColors), lw = 4)
   ax.plot(densityPol2_2['iR'],densityPol1_1['density'], label = 'neutral, Iteration ' + iteration, c = colors((colorCounter + 1)/nColors), lw = 4)
   #plt.axvline(x = inte[-1]**(1/3), color = colors((colorCounter + 2)/nColors), linestyle='dashed',label = 'radio iso')
   rmario =buscoR(densityPol1_1['density'])
   plt.axvline(x = rmario, color = 'b', linestyle='dashed',label = 'radio')
   densityCORE=  densityPol1_2['density'] + densityPol2_2['density']
   inte_core = calculoR(densityCORE)
   rmario_core =buscoR(densityCORE)

   plt.axvline(x = inte_core[-1]**(1/3), color = 'b', linestyle='dashed',label = 'radio core')
   plt.axvline(x = rmario_core, color = 'r', linestyle='dashed',label = 'radio core mario')
   colorCounter += 4

#Settings del grafico
#ax.set_title('density pol', fontsize = 30)
ax.set_xlabel('$r$ ($nm$)', fontsize = 30)
ax.set_ylabel('Density polymer', fontsize = 30)
#ax.axvline(x = 5, c= 'k', lw = 3, ls = 'dashed')
ax.set_xlim([0.1, 140*0.2])
#ax.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0))
#ax.yaxis.get_offset_text().set_fontsize(20)
ax.tick_params(axis = 'both', labelsize = 20)
ax.legend(fontsize = 16)

plt.show()

print('el radio iso es ', inte[-1]**(1/3), 'nm')
print('el radio mario es ', rmario, 'nm')

print('el radio de core de  mario es ', rmario_core, 'nm')
print('el radio de la parte neutra es ', rmario-rmario_core, 'nm')


