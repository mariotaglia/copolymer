import subprocess
import os
import re
import sys
import time

definitions = "DEFINITIONS.txt"

copolymer = sys.argv[1]

os.chdir("ramp")

#find maximum value for dimensions and replace by 1

with open(definitions, "r") as f:
    contenido = f.readlines()

with open(definitions, "w") as f:
    for line in contenido:
        if "dimensions" in line:
             line_split = line.split()
             max_ntot = int(line_split[3])
             line_split[3] = "1"
             line = " ".join(line_split) + "\n"
        f.write(line)

#correr copolymer
os.system(copolymer)
sys.stdout.flush()

#result = subprocess.run(copolymer, shell=True)
#if result.returncode != 0:
#    print("Error en copolymer")
#    exit()

#cambiar ultimo out por in.in
archivos_out = [f for f in os.listdir() if re.match(r"out\.001\.\d+\.dat", f)]
if not archivos_out:
    print("No se encontraron out.dat")
    sys.stdout.flush()
    exit()

out_nums = [int(re.findall(r'\d+', f)[-1]) for f in archivos_out]
ultimo_out = max(out_nums)
ultimo_archivo_out = f"out.001.{ultimo_out:03d}.dat"

subprocess.run(f"mv {ultimo_archivo_out} in.in", shell=True, check=True)
print(f"out movido a in.in")
sys.stdout.flush()
time.sleep(20)

#cambiar infilie 0 a 2
with open(definitions, "r") as f:
    contenido_definitions = f.read()

contenido_definitions = re.sub(r"(infile\s+)[^\n]+", r"\g<1>2", contenido_definitions)

with open(definitions, "w") as f:
    f.write(contenido_definitions)

#reemplazar npol del ultimo system y cambiar min y max a npol en DEFINITIONS
archivos_system = [f for f in os.listdir() if re.match(r"system\.001\.\d+\.dat", f)]
if not archivos_system:
    print("No se encontraron system.dat")
    sys.stdout.flush()
    exit()

system_nums = [int(re.findall(r'\d+', f)[-1]) for f in archivos_system]
ultimo_system = max(system_nums)
ultimo_archivo_system = f"system.001.{ultimo_system:03d}.dat"

with open(ultimo_archivo_system, "r") as f:
    linea_npol = f.readlines()[13].strip()
valor_npol = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", linea_npol)[0]

with open(definitions, "r") as f:
    contenido = f.read()

def reemplazar_npol(match):
    valores = match.group(1).split()
    valores[0] = valor_npol
    valores[1] = valor_npol
    valores[2] = valor_npol
    return "npol " + " ".join(valores)

contenido = re.sub(r"npol\s+([^\n]+)", reemplazar_npol, contenido)

with open(definitions, "w") as f:
    f.write(contenido)

print(f"Reemplazado npol por {valor_npol}.")
sys.stdout.flush()

for i in range(0,max_ntot-1):

    #sumarle 1 a dimensions
    with open(definitions, "r") as f:
        contenido = f.readlines()

    with open(definitions, "w") as f:
        for line in contenido:
            if "dimensions" in line:
             line_split = line.split()
             line_split[3] = str(int(line_split[3]) + 1)
             line = " ".join(line_split) + "\n"
            f.write(line)

    print(f"sumado 1 a dimensions")
    sys.stdout.flush()

    #correr copolymer
    os.system(copolymer)
#    result = subprocess.run(copolymer, shell=True)
#    print(result)
#    if result.returncode != 0:
#        print("Error en copolymer")
#        exit()
#    subprocess.run(f"mv out.001.001.dat in.in", shell=True, check=True)

    if not os.path.isfile("out.001.001.dat"):
        print("out.001.001.dat not found. exit")
        sys.stdout.flush()
        exit()
    os.system("mv out.001.001.dat in.in")
    print("out.001.001.dat movido a in.in")
    sys.stdout.flush()
    time.sleep(20)

# not go back one folder and scan densirty

os.chdir("..")
os.system("mv ./ramp/in.in .")
time.sleep(20)

# replace initial value

with open(definitions, "r") as f:
    contenido = f.read()

def reemplazar_npol(match):
    valores = match.group(1).split()
    valores[0] = valor_npol
    return "npol " + " ".join(valores)

contenido = re.sub(r"npol\s+([^\n]+)", reemplazar_npol, contenido)

with open(definitions, "w") as f:
    f.write(contenido)

#
# run one last time
#

os.system(copolymer)
sys.stdout.flush()


