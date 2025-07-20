import pandas as pd
import argparse
import os
import numpy as np
import shutil
import re


#ejecutar como python3 PA.py [código de PA]
parser = argparse.ArgumentParser(description="Crear DEFINITIONS según código de péptidoanfifilo")
parser.add_argument("codigo", type=str, help="Código alfanumérico péptidoanfifilo")
args = parser.parse_args()
PA_codigo = args.codigo


aa = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "Z"]
beads = [2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1]
aa_acidos = ["D", "C", "E", "Y"]
aa_basicos = ["A", "H", "K"]
pka_lista = [3.9, 8.37, 4.07, 10.5]
pkb_lista  = [1.52, 7.96, 3.46]
vpol = [0.048, 0.157, 0.076, 0.070, 0.069, 0.094, 0.103, None, 0.112, 0.126, 0.126, 0.128, 0.122, 0.149, 0.085, 0.047, 0.075, 0.186, 0.153, 0.099, 0.113]
aa_beads = ["C1", "C3", "P5", "P3", "C5", "C3", "P4", None, "P5","C1", "C1", "C3", "C5", "C4", "C2", "P1", "P1", "C5", "C5", "C1", "C1"]
nombres = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 'His', 'Ile','Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val', 'HC']
aa_BB1 = ["A","R","N","D","E","Q","G","H","K","P","S"]
aa_BB2 = ["C","I","L","M","F","T","W", "Y","V"]


pepanf = []

for i in PA_codigo:
    if i in aa:
       pepanf.append(i)

fullbeads = []
for i in pepanf:
    if i == "Z":
        fullbeads.append(i)
    elif i in aa_BB1:
        fullbeads.append("BB1")
    elif i in aa_BB2:
        fullbeads.append("BB2")
for i in pepanf:
    if i not in ("Z","G"):
        fullbeads.append(i)

long_valor = sum([beads[aa.index(i)] for i in pepanf])

###DEFINITIONS

#pka y pkb
pka = []
pkb = []
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

#vpol
vpol_valores = []
BB1_seen = False
BB2_seen = False
aa_seen = []

for i in fullbeads:
    if i == "BB1":
        if not BB1_seen:
            BB1_seen = True
            vpol_valores.append(0.113)
    elif i == "BB2":
        if not BB2_seen:
            BB2_seen = True
            vpol_valores.append(0.113)
    elif i not in aa_seen:
        aa_seen.append(i)
        vpol_valores.append(vpol[aa.index(i)])

#matrices diagonales
dimf = np.tril(np.ones((len(vpol_valores),len(vpol_valores)), dtype=int)) * 6
lseg = np.tril(np.ones((len(vpol_valores),len(vpol_valores)), dtype=int)) * 0.47
Utg = np.zeros((len(vpol_valores),2), dtype=int)

#carboxi y amino terminal, el 1ro indica amino, el 2do carboxi
terminal = re.findall(r'\d+', PA_codigo)
pkbterminal = None
pkaterminal = None
for i in terminal:
    if i == "00":
        pass
    elif i == "01":
        if pepanf[-1] == "Z":
            raise SystemExit("Z en el extremo derecho impide carboxilo libre")
        pkaterminal = 4.5
    elif i == "10":
        if pepanf[0] == "Z":
            raise SystemExit("Z en el extremo izquierdo impide amino libre")
        pkbterminal = 4.5
    elif i == "11":
        if pepanf[-1] == "Z":
            if pepanf[0] != "Z":
                raise SystemExit("Z en el extremo derecho impide carboxilo libre")
            if pepanf[0] == "Z":
                raise SystemExit("Z en los extremos impide amino libre y carboxilo libre")
        if pepanf[0] == "Z":
            raise SystemExit("Z en el extremo izquierdo impide amino libre")
        pkbterminal = 4.5
        pkaterminal = 4.5


###STRUCTURE

matriz_structure = np.zeros((long_valor,3))
bead_numero = []
bead_norepetido = {}
for j, k in enumerate(fullbeads):
    if not bead_norepetido:
        bead_norepetido[k] = 1
    elif k not in bead_norepetido:
        max_val = max(bead_norepetido.values())
        bead_norepetido[k] = max_val + 1
    bead_numero.append(bead_norepetido[k])

pka_vistos = []
k = 0
pkb_vistos = []
n = 0
for i in range(matriz_structure.shape[0]):
    matriz_structure[i,0] = bead_numero[i]
    j = fullbeads[i]
    if j in aa_acidos:
        if j not in pka_vistos:
            k = k + 1
            matriz_structure[i,1] = k
            pka_vistos.append(j)
        else:
            matriz_structure[i,1] = pka_vistos.index(j) + 1
    if j in aa_basicos:
        if j not in pkb_vistos:
            n = n + 1
            matriz_structure[i,2] = n
            pkb_vistos.append(j)
        else:
            matriz_structure[i,2] = pkb_vistos.index(j) + 1

BB_indices = [i for i, bead in enumerate(fullbeads) if bead in ("BB1", "BB2")]
primer_BB = BB_indices[0]
ultimo_BB = BB_indices[-1]
for i in terminal:
    if i == "00":
        pass
    elif i == "01":
        matriz_structure[ultimo_BB, 1] = len(pka) + 1
    elif i == "10":
        matriz_structure[primer_BB, 2] = len(pkb) + 1
    elif i == "11":
        matriz_structure[ultimo_BB, 1] = len(pka) + 1
        matriz_structure[primer_BB, 2] = len(pkb) + 1


###EPSILON

#matriz epsilon
beadslist = ["Qda","Qd","Qa","Q0","P5","P4","P3","P2","P1","Nda","Nd","Na","N0","C5","C4","C3","C2","C1","SC4","EO"]
n = len(vpol_valores) + 1
actual_path = os.path.dirname(os.path.abspath(__file__))
ruta_martini_tabla = os.path.join(actual_path, "table_martini.dat")
data=np.loadtxt(ruta_martini_tabla, skiprows=1, usecols=range(1,21))
epslist = data[-2,:]
sigmalist = data[-1,:]
epsilon=np.zeros((n,n))
epsilon_th=np.zeros((n,n))
interaction_index=np.zeros((n,n),dtype=np.int8)
sigma=np.zeros((n,n))
vol_mar=np.zeros((n,n))

bead_type_index = []
aa_repetidos = []
aa_to_beads = ["P4"]
vol_th = [0.03]
for i in fullbeads:
    if i not in aa_repetidos:
        if i == "BB1":
            aa_to_beads.append("P5")
            vol_th.append(0.113)
        elif i == "BB2":
            aa_to_beads.append("Nda")
            vol_th.append(0.113)
        else:
            vol_th.append(vpol[aa.index(i)])
            aa_to_beads.append(aa_beads[aa.index(i)])
        aa_repetidos.append(i)

for i in aa_to_beads:
        bead_type_index.append(beadslist.index(i))

for i in range(0,n):
   for j in range(0,n):
     interaction_index[i][j]=data[bead_type_index[i]][bead_type_index[j]]
     epsilon[i][j]=epslist[interaction_index[i][j]]
     sigma[i][j]=sigmalist[interaction_index[i][j]]
     vol_mar[i][j]=(sigma[i][j]/2)**3*4/3*np.pi


for i in range(0,n):
   for j in range(0,n):
      epsilon_th[i][j]=epsilon[i][j]/(vol_mar[i][i]*vol_mar[j][j])-epsilon[i][0]/(vol_mar[0][0]*vol_mar[j][j])-epsilon[0][j]/(vol_mar[i][i]*vol_mar[0][0])+epsilon[0][0]/(vol_mar[0][0]**2)
      epsilon_th[i][j]=epsilon_th[i][j]*vol_th[i]*vol_th[j]

factor = 0.36
epsilon_th *= factor
epsilon_th = epsilon_th[1:, 1:]

#texto epsilon

epsilon_nombres = []
epsilon_beads = []
epsilon_texto = []
for i in bead_norepetido.keys():
    if i == "BB1":
        epsilon_nombres.append("BB1")
        epsilon_beads.append("P5")
    elif i == "BB2":
        epsilon_nombres.append("BB2")
        epsilon_beads.append("Nda")
    else:
        epsilon_nombres.append(nombres[aa.index(i)])
        epsilon_beads.append(aa_beads[aa.index(i)])

for k, (i, j) in enumerate(zip(epsilon_nombres, epsilon_beads), start=1):
    epsilon_texto.append(f"{k}:   {i}   {j}")


#creo una carpeta con el nomnbre del codigo
dir_general = f"{PA_codigo}"
os.makedirs(dir_general, exist_ok=True)


#creo las carpetas curvatura
script_dir = os.path.dirname(os.path.abspath(__file__))
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

Npoorsv {len(vpol_valores)}
vpol {len(vpol_valores)}
{chr(10).join(str(j) for j in vpol_valores if j is not None)}

rsalt 0.3 0.3

dimf {len(vpol_valores)}
{chr(10).join(" ".join(str(j) for j in fila if j != 0) for fila in dimf)}

Nacids {len(pka) + (1 if pkaterminal is not None else 0)}
{chr(10).join(str(j) for j in pka)}
{pkaterminal if pkaterminal is not None else ''}
Nbasics {len(pkb) + (1 if pkbterminal is not None else 0)}
{chr(10).join(str(j) for j in pkb)}
{pkbterminal if pkbterminal is not None else ''}

PBCflag 1
infile 0
flagkai 0
flagonekais 1

npol {npol}

Xulimit 5
lseg {len(vpol_valores)}
{chr(10).join(" ".join(str(j) for j in fila if j != 0) for fila in lseg)}

lsegkai {len(vpol_valores)}
{chr(10).join(" ".join(str(j) for j in fila if j != 0) for fila in lseg)}

Utg {len(vpol_valores)}
{chr(10).join(" ".join(f"{j:.2f}" for j in fila) for fila in Utg)}

csalt 0.1
pHbulk 7
dielP 3.0

nbranches
{len(nbranches)}
{chr(10).join(str(j) for j in nbranches)}

saveflag 0
"""

    ruta_definitions = os.path.join(ruta_curvatura, "DEFINITIONS.txt")
    with open(ruta_definitions, 'w', encoding='utf-8') as f:
        f.write(definitions)

    structure = f"""{chr(10).join(chr(9).join(f"{int(i)}" for i in fila) for fila in matriz_structure)} """

    ruta_structure = os.path.join(ruta_curvatura, "structure.001.in")
    with open(ruta_structure, 'w', encoding='utf-8') as f:
        f.write(structure)


    epsilon = f"""
{chr(10).join(" ".join(f"{j}" for j in fila) for fila in epsilon_th)}

0:   W    P4
{chr(10).join(epsilon_texto)}
"""

    ruta_epsilon = os.path.join(ruta_curvatura, "epsilon.in")
    with open(ruta_epsilon, 'w', encoding='utf-8') as f:
        f.write(epsilon)


    archivo_kais = os.path.join(script_dir, f"kais{i}", "kais.001.001.in")
    shutil.copy(archivo_kais, ruta_curvatura)



