import parmed as pmd

#convert prmtop and inpcrd into top and gro
amber = pmd.load_file('c16k3_H.prmtop','c16k3_H.rst7')
#
amber.save('c16k3_H.top')
amber.save('c16k3_H.crd')
#
