import subprocess
import os
import re
import sys
import time

definitions = "DEFINITIONS.txt"

copolymer = sys.argv[1]

#######################################################################
# go to ramp folder
os.chdir("ramp")

#######################################################################
# find maximum value for dimensions in DEFINITIONS.txt and replace by 1

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

#######################################################################
# run copolymer
os.system(copolymer)
sys.stdout.flush()


#######################################################################
# At this point copolymer crashes 
# We will attempt now ramping maxntot from 1 using the previos solutions

# create a folder and move solutions there


os.mkdir("first_scan")
os.system("mv system.* first_scan")
os.system("mv out.* first_scan")
os.system("rm *.dat")

# find number of out and system files
archivos_out = [f for f in os.listdir("./first_scan") if re.match(r"out\.001\.\d+\.dat", f)]
if not archivos_out:
    print("##### run.py: No se encontraron out.dat #####")
    sys.stdout.flush()
    exit()

out_nums = [int(re.findall(r'\d+', f)[-1]) for f in archivos_out]
ultimo_out = max(out_nums) 

archivos_system = [f for f in os.listdir("./first_scan") if re.match(r"system\.001\.\d+\.dat", f)]
if not archivos_system:
    print("##### run.py: No se encontraron system.dat #####")
    sys.stdout.flush()
    exit()

system_nums = [int(re.findall(r'\d+', f)[-1]) for f in archivos_system]
ultimo_system = max(system_nums) 

if(ultimo_system != ultimo_out):
    print("##### run.py: Number of out and system files differ. Exiting #####")
    sys.stdout.flush()
    exit()



#######################################################################
# loop until find solution or first solution was used

converged = 0 # flag
while (converged == 0):

#traer ultimo out

    ultimo_archivo_out = f"first_scan/out.001.{ultimo_out:03d}.dat"

    subprocess.run(f"mv {ultimo_archivo_out} out.001.001.dat", shell=True, check=True)
    print(f"##### run.py: out movido #####")
    sys.stdout.flush()
    time.sleep(20)

#cambiar infilie 0 a 2
    with open(definitions, "r") as f:
        contenido_definitions = f.read()

    contenido_definitions = re.sub(r"(infile\s+)[^\n]+", r"\g<1>2", contenido_definitions)

    with open(definitions, "w") as f:
        f.write(contenido_definitions)

#reemplazar npol del ultimo system y cambiar min y max a npol en DEFINITIONS
    ultimo_archivo_system = f"first_scan/system.001.{ultimo_system:03d}.dat"

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

    print(f"##### run.py: reemplazado npol por {valor_npol}. #####")
    sys.stdout.flush()

    for i in range(0,max_ntot-1):

    # mv out -> in for next iteration

        if not os.path.isfile("out.001.001.dat"):
            sys.stdout.flush()
            break
        os.system("mv out.001.001.dat in.in")
        print("##### run.py: out.001.001.dat movido a in.in #####")
        sys.stdout.flush()
        time.sleep(20)


    #dimensions = i+2
        with open(definitions, "r") as f:
            contenido = f.readlines()

        with open(definitions, "w") as f:
            for line in contenido:
                if "dimensions" in line:
                    line_split = line.split()
                    line_split[3] = str(i+2)
                    line = " ".join(line_split) + "\n"
                f.write(line)

        print(f"##### run.py: changed dimensions #####")
        sys.stdout.flush()

    #correr copolymer
        os.system(copolymer)

# here ends loop converged
#############################################################################################

# now go back one folder and scan density


    # check if converged or not

    if os.path.isfile("out.001.001.dat"):
        print("##### run.py: loop in maxntot converged #####")
        converged = 1
    else:
        print("##### run.py: loop in maxntot not converged #####")
        if(ultimo_system == 1):
            print("##### run.py: no more solutions to try -- exiting #####")
            sys.stdout.flush()
            exit()
        print("##### run.py: decreasing npol and trying again #####")
        ultimo_system = ultimo_system - 1
        ultimo_out = ultimo_out - 1
        os.system("rm *.dat")


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


