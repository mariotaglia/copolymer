import sys
import pandas as pd
import numpy as np

def cargaDatos(filename):
    # Cargo los datos
    datos = pd.read_csv(filename, sep = '\s+', header = None)
    datos.columns = ['iR', 'iZ', 'density']
    datos['r2*density'] = datos['density']*datos['iR']**2#*0.2**2
    return datos

def formatThreeDigit(number):
    return '0'*(3 - len(number)) + number

for i, iteration in enumerate(sys.argv[1:]):
   # Cargo los datos
   iterCode = formatThreeDigit(iteration)
   xdensityPol1 = cargaDatos('xdensitytota.001.001.' + iterCode + '.dat')
   xdensityPol2 = cargaDatos('xdensitytota.002.001.' + iterCode + '.dat')
   # Encuentro el criterio de integracion
   corte = 0.001
   maxDensity1 = xdensityPol1['density'].max()*corte
   maxDensity2 = xdensityPol2['density'].max()*corte
   limite1 = xdensityPol1[xdensityPol1['density'] > maxDensity1]['iR'].max()
   limite2 = xdensityPol2[xdensityPol2['density'] > maxDensity2]['iR'].max()
   limite = np.max([limite1, limite2])
   print('Limite Pol 1:', limite1)
   print('Limite Pol 2:', limite2)
   # Integro numericamente hasta o desde el limite hallado
   integral1_bajo = 4*np.pi*np.trapz(xdensityPol1[xdensityPol1['iR'] <= limite]['r2*density'], xdensityPol1[xdensityPol1['iR'] <= limite]['iR'])
   integral1_alto = 4*np.pi*np.trapz(xdensityPol1[xdensityPol1['iR'] > limite]['r2*density'], xdensityPol1[xdensityPol1['iR'] > limite]['iR'])
   integral2_bajo = 4*np.pi*np.trapz(xdensityPol2[xdensityPol2['iR'] <= limite]['r2*density'], xdensityPol2[xdensityPol2['iR'] <= limite]['iR'])
   integral2_alto = 4*np.pi*np.trapz(xdensityPol2[xdensityPol2['iR'] > limite]['r2*density'], xdensityPol2[xdensityPol2['iR'] > limite]['iR'])

print('Integral hasta', limite, ': Cop:', integral1_bajo, ', Mol:', integral2_bajo)
print('Integral desde', limite, ': Cop:', integral1_alto, ', Mol:', integral2_alto)
print('Integral total:', integral1_bajo + integral1_alto + integral2_bajo + integral2_alto)

