import os
import numpy as np
import pandas as pd

listfe = []
listindex = []
listnpol = []


for filename in os.listdir("."):
    if filename.endswith(".dat") and filename.startswith("system"):
        index = filename[7:14]
        listindex.append(index)
        with open(filename) as fp:
            lines = fp.readlines()
            for line in lines:
                if(line.startswith(" Free")):
                     listfe.append(float(line.split()[3]))
                if(line.startswith(" npol ")):
                     listnpol.append(line.split()[2])
                     
pos = listfe.index(min(listfe))
index = listindex[pos]
npol = listnpol[pos]

print(pos,index,npol,min(listfe))

os.system("mkdir equil")
os.chdir("./equil")
os.system("cp ../*.in .")
os.system("cp ../out."+index+".dat in.in")


with open("../DEFINITIONS.txt") as fp:
    lines = fp.readlines()
with open("DEFINITIONS.txt","w") as fp:
    for line in lines:
        if 'npol' in line:
            line = "npol "+npol+" "+npol+" "+npol+" 0.01"
        if 'saveflag' in line:
            line = "saveflag 1"
        fp.write(line)

