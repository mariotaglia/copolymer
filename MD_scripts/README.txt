1. Pasar el formato XX.mol a XX.pdb con AVOGADRO
2. pasar pdb a mol2: antechamber -i XX.pdb -fi pdb -o XX.mol2 -fo mol2 
3. da parametros:  parmchk2 -i XX.mol2 -f mol2 -o XX.prm -a 'Y' 
4. antechamber -i XX.mol2 -fi mol2 -o XX-fin.pdb -fo pdb
5. antechamber -i XX.mol2 -fi mol2 -o XX-fin.prepi -fo prepi
6. tleap script "tleap.in" 
     source leaprc.gaff2
     loadAmberParams XX.prm
     loadamberprep XX-fin.prepi
     c16 = loadpdb XX-fin.pdb
     saveamberparm c16 XX.prmtop XX.rst7
     quit
7.# Correr tleap
tleap -f tleap.in 
8. se necesita el archivo XX.rst7 y XX.prmtop y correr convert.py. Para poder hacer esto, se requiere:
	a. tener instalado ambertool
        b. tener instalado parmed (pip install parmed)
   se corre con python3.9

   esto genera el .crd y top.
   

9. Para armar el data para lammps: correr amber2lammps.py  que necesita XX.crd  XX.prmtop y XX.top
	modificar antes de los archivos:
          a. crd: 4 lineas vacias iniciales y en la 5 SOLO el numero de atomos
          b. .prmtop  : borrar las primeras 6 lineas que compienzan con  %. Y poner una linea con XX

     ejecutar con python2 (kMC environment)
