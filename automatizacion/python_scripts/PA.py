import pandas as pd
import argparse
import os 
import numpy as np
import shutil
import re


#agregar BB correspondiente a cada aa
#npoorsv vpol lseg etc van aa sin repetir y BB sin repetir 

#ejecutar como python3 definitions_script.py [código de PA]
parser = argparse.ArgumentParser(description="Crear DEFINITIONS según código de péptidoanfifilo")
parser.add_argument("codigo", type=str, help="Código alfanumérico péptidoanfifilo")
args = parser.parse_args()
PA_codigo = args.codigo


aa = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "Z", "O", "B", "J"]
beads = [2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2]
aa_acidos = ["D", "C", "E", "Y"]
aa_basicos = ["R", "H", "K", "O", "B", "J"]
pka_lista = [3.9, 8.37, 4.07, 10.5]
pkb_lista  = [1.52, 7.96, 3.46, 3.5, 3.73, 4.57]
long_valor = sum([beads[aa.index(s)] for s in PA_codigo])
vpol = [0.048, 0.157, 0.076, 0.070, 0.069, 0.094, 0.103, None, 0.112, 0.126, 0.126, 0.128, 0.122, 0.149, 0.085, 0.047, 0.075, 0.186, 0.153, 0.099, 0.113, 0.086, 0.059, 0.032]
aa_beads = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'Ñ', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W'] ###modificar
nombres = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 'His', 'Ile','Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val', 'HC', "Orn", "Dab", "Dap"] ###faltaria BB que no puedo ponerlo

pepanf = [] 
pka = []
pkb = []

for i in PA_codigo:
    if i in aa:
       pepanf.append(i)

for i in pepanf:
    if i in aa_acidos:
        valor_pka = pka_lista[aa_acidos.index(i)]
        if valor_pka not in pka:
            pka.append(valor_pka)
    if i in aa_basicos:
        valor_pkb = pkb_lista[aa_basicos.index(i)]
        if valor_pkb not in pkb:
            pkb.append(valor_pkb)

#nbranches
nbranches = []
for i, j in enumerate(pepanf, start=1):
    if j in ("Z", "G"):
        pass
    else:
        nbranches.append(f"{i} 1")

#matrices diagonales
aa_sZ = []

for i in pepanf:
    if i == "Z":
        pass
    else:
        aa_sZ.append(i)

dimf = np.tril(np.ones((len(aa_sZ),len(aa_sZ)), dtype=int)) * 6
lseg = np.tril(np.ones((len(aa_sZ),len(aa_sZ)), dtype=int)) * 0.47
Utg = np.zeros((len(aa_sZ),2), dtype=int) 

#vpol
vpol_valores = []
for i, j in enumerate(pepanf, start=1):
    if j == "G":
        pass
    elif vpol[aa.index(j)] not in vpol_valores:
        vpol_valores.append(vpol[aa.index(j)])

#creo una carpeta con el nomnbre del codigo 
dir_general = f"{PA_codigo}"
os.makedirs(dir_general, exist_ok=True)

#carboxi y amino terminal
terminal = re.findall(r'\d+', PA_codigo)
pkbterminal = None
pkaterminal = None
for i in terminal:
    if i == "00":
        pass
    elif i == "01":
        pkbterminal = 4.5
    elif i == "10":
        pkaterminal = 4.5
    elif i == "11":
        pkbterminal = 4.5
        pkaterminal = 4.5   

"""
#structure
matriz_structure = np.zeros((long_valor,3))
for i in range(matriz_structure.shape[0]):
    for j in pepanf:
"""


#epsilon
epsilon_nombres = []
epsilon_beads = []
for i in pepanf:
    valor_nombre = nombres[aa.index(i)]
    if valor_nombre not in epsilon_nombres:
        epsilon_nombres.append(valor_nombre)
    valor_beads = aa_beads[aa.index(i)]
    if valor_beads not in epsilon_beads:
        epsilon_beads.append(valor_beads)

epsilon_texto = []
for k in range(len(epsilon_nombres)):
  for i, j in zip(epsilon_nombres, epsilon_beads):
    epsilon_texto.append(f"{k+2}:   {i}   {j}")



#creo las carpetas curvatura
for i in range(0,3):
    dir_curvatura = f"curvatura{i}"
    ruta_curvatura = os.path.join(dir_general, dir_curvatura)
    os.makedirs(ruta_curvatura, exist_ok=True)
    if i == 0:
        npol = f"0.01 0.01 2.5 0.01"
    if i == 1:
        npol = f"1 1 20 0.1"
    if i == 2:
        npol = f"1 1 50 0.5"

    definitions = f"""
Ncomp 1

npolratio 1.

curvature {i}
dimensions 140 1 1 1
cuantas
10000
long
{long_valor}

layersize 0.2 1

Npoorsv {len()}
vpol {len(nbranches)}
0.113
{"\n".join(str(j) for j in vpol_valores)}

rsalt 0.3 0.3

dimf {len(aa_sZ)}
{"\n".join(" ".join(str(j) for j in fila if j != 0) for fila in dimf)}

Nacids 1
{"\n".join(str(j) for j in pka)}
{pkaterminal if pkaterminal is not None else ''}
Nbasics 1
{"\n".join(str(j) for j in pkb)}
{pkbterminal if pkbterminal is not None else ''}

PBCflag 1
infile 0
flagkai 0

npol {npol}

Xulimit 5
lseg {len(aa_sZ)}
{"\n".join(" ".join(str(j) for j in fila if j != 0) for fila in lseg)}

lsegkai {len(aa_sZ)}
{"\n".join(" ".join(str(j) for j in fila if j != 0) for fila in lseg)}

Utg {len(aa_sZ)}
{"\n".join(" ".join(f"{j:.2f}" for j in fila) for fila in Utg)}

csalt 0.1
pHbulk 7
dielP 3.0

nbranches
{len(nbranches)}
{"\n".join(str(j) for j in nbranches)}

saveflag 0
"""

    ruta_definitions = os.path.join(ruta_curvatura, "DEFINITIONS.txt")
    with open(ruta_definitions, 'w', encoding='utf-8') as f:
        f.write(definitions)

    structure = f"""

"""

    ruta_structure = os.path.join(ruta_curvatura, "structure.001.in")
    with open(ruta_structure, 'w', encoding='utf-8') as f:
        f.write(structure)


    epsilon = f"""


0:00    W   P41:    BB  Nda{" ".join(epsilon_texto)}
"""

    ruta_epsilon = os.path.join(ruta_curvatura, "epsilon.in")
    with open(ruta_epsilon, 'w', encoding='utf-8') as f:
        f.write(epsilon)

