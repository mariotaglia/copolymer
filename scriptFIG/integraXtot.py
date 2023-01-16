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
   # Calculo de la derivada
   dy, dx = np.gradient(xdensityPol1['density'] + xdensityPol2['density']), np.gradient(xdensityPol1['iR'])
   derivative = dy/dx
   sizeDeriv = derivative.size
   criterio = 1e-5
   for i in range(sizeDeriv):
       if np.abs(derivative[sizeDeriv - i - 16]) > criterio:
           argLimite = sizeDeriv - i - 16
           break
   limite = xdensityPol1['iR'].iloc[argLimite]
   # Integro numericamente hasta o desde el limite hallado
   integral1_bajo = 4*np.pi*np.trapz(xdensityPol1[xdensityPol1['iR'] <= limite]['r2*density'], xdensityPol1[xdensityPol1['iR'] <= limite]['iR'])
   integral1_alto = 4*np.pi*np.trapz(xdensityPol1[xdensityPol1['iR'] >= limite]['r2*density'], xdensityPol1[xdensityPol1['iR'] >= limite]['iR'])
   integral2_bajo = 4*np.pi*np.trapz(xdensityPol2[xdensityPol2['iR'] <= limite]['r2*density'], xdensityPol2[xdensityPol2['iR'] <= limite]['iR'])
   integral2_alto = 4*np.pi*np.trapz(xdensityPol2[xdensityPol2['iR'] >= limite]['r2*density'], xdensityPol2[xdensityPol2['iR'] >= limite]['iR'])

print('Integral hasta', limite, ': Cop:', integral1_bajo, ', Mol:', integral2_bajo)
print('Integral desde', limite, ': Cop:', integral1_alto, ', Mol:', integral2_alto)
print('Integral total:', integral1_bajo + integral1_alto + integral2_bajo + integral2_alto)

