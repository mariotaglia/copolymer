import pandas as pd
import os
import matplotlib.pyplot as plt

def extract_min_values(base_folder_path, subfolders):
   min_values = []
   for subfolder in subfolders:
       folder_path = os.path.join(base_folder_path, subfolder)
       file_path = os.path.join(folder_path, 'F_tot.dat')
       
       if os.path.exists(file_path):
           df = pd.read_csv(file_path, sep=r'\s+', header=None, names=['N', 'F'])
           min_value = df['F'].min()
           min_index = df['F'].idxmin()
           min_N = df.loc[min_index, 'N']
           min_values.append(min_value)
   return min_values

base_folder_path_curvatura0 = os.path.expanduser('~/projects/C16K3-COONH2/C16K3_newmodel_opcion1/curvatura0/')
base_folder_path_curvatura1 = os.path.expanduser('~/projects/C16K3-COONH2/C16K3_newmodel_opcion1/curvatura1/')
base_folder_path_curvatura2 = os.path.expanduser('~/projects/C16K3-COONH2/C16K3_newmodel_opcion1/curvatura2/')

subfolders = ['ph6', 'ph6.5', 'ph7', 'ph7.5', 'ph8', 'ph8.5', 'ph9', 'ph9.5', 'ph10', 'ph10.5', 'ph11', 'ph11.5', 'ph12']

min_values_curvatura0 = extract_min_values(base_folder_path_curvatura0, subfolders)
min_values_curvatura1 = extract_min_values(base_folder_path_curvatura1, subfolders)
min_values_curvatura2 = extract_min_values(base_folder_path_curvatura2, subfolders)

pH = [6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12]

plt.scatter(pH, min_values_curvatura0, color='blue', label='Lamelas')
plt.scatter(pH, min_values_curvatura1, color='red', label='Fibras')
plt.scatter(pH, min_values_curvatura2, color='green', label='Micelas')

plt.xlabel('pH')
plt.ylabel('F')
plt.legend()

plt.savefig('valores_minimos.png')
