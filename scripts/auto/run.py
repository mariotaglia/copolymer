import subprocess
import os
import re


definitions = "DEFINITIONS.txt"


#correr copolymer
result = subprocess.run("~/develop/copolymer/assembly", shell=True)
if result.returncode != 0:
    print("Error en copylimer")
    exit()



#cambiar ultimo out por in.in
archivos_out = [f for f in os.listdir() if re.match(r"out\.001\.\d+\.dat", f)]
if not archivos_out:
    print("No se encontraron out.dat")
    exit()

out_nums = [int(re.findall(r'\d+', f)[-1]) for f in archivos_out]
ultimo_out = max(out_nums)
ultimo_archivo_out = f"out.001.{ultimo_out:03d}.dat"

subprocess.run(f"mv {ultimo_archivo_out} in.in", shell=True, check=True)
print(f"out movido a in.in")



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




#sumarle 1 a dimensions
with open(definitions, "r") as f:
    contenido = f.read()
linea_dimensions = contenido.splitlines()[7]
valor_dimensions = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", linea_dimensions)[2]
valor_dimensions = int(valor_dimensions) + 1
 
def reemplazar_dimensions(match):
    valores = match.group(1).split()
    valores[2] = str(valor_dimensions)
    return "dimensions " + " ".join(valores)
 
contenido = re.sub(r"dimensions\s+([^\n]+)", reemplazar_dimensions, contenido)
 
with open(definitions, "w") as f:
    f.write(contenido)
 
print(f"Sumado una unidad a dimensions.")




for i in range(0,16):
    result = subprocess.run("~/develop/copolymer/assembly", shell=True)
    if result.returncode != 0:
        print("Error en copylimer")
        exit()

    with open(definitions, "r") as f:
        contenido = f.read()
    linea_dimensions = contenido.splitlines()[7]
    valor_dimensions = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", linea_dimensions)[2]
    valor_dimensions = int(valor_dimensions) + 1
 
    def reemplazar_dimensions(match):
        valores = match.group(1).split()
        valores[2] = str(valor_dimensions)
        return "dimensions " + " ".join(valores)
 
    contenido = re.sub(r"dimensions\s+([^\n]+)", reemplazar_dimensions, contenido)
 
    with open(definitions, "w") as f:
        f.write(contenido)
    print(f"sumado 1 a dimensions")
    subprocess.run(f"mv out.001.001.dat in.in", shell=True, check=True)
    print(f"out.001.001.dat movido a in.in")


