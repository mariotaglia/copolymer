import numpy

beadslist=["Qda","Qd","Qa","Q0","P5","P4","P3","P2","P1","Nda","Nd","Na","N0","C5","C4","C3","C2","C1"]

print("------------------------------------------------------------------------------------------------")
print("Possible types of martini beads: Qda  Qd  Qa  Q0  P5  P4  P3  P2  P1  Nda  Nd  Na  N0  C5  C4  C3  C2  C1")
print("Please use exactly the same spelling, including lower and upper cases.")
print("------------------------------------------------------------------------------------------------")
print("")

n=0
n = input("Number of types of beads (including solvent): ")
print("")

name=""
bead_type_index=[0]*n
vol_th=[0]*n

for i in range(0,n):
   print("Name of bead #%s: " %i)
   name_tmp = raw_input()
   print("Martini type: ")
   type_tmp = raw_input()
   vol_th[i] = input("Molecular theory volume (nm): ")
   print("")
   bead_type_index[i]=beadslist.index(type_tmp)
   name= name+str(i)+": "+name_tmp+"	"

data=numpy.loadtxt("table_martini.dat", skiprows=1, usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18))
epslist=data[:][18]
sigmalist=data[:][19]

epsilon=numpy.zeros((n,n))
epsilon_th=numpy.zeros((n,n))
interaction_index=numpy.zeros((n,n),dtype=numpy.int8)
sigma=numpy.zeros((n,n))
vol_mar=numpy.zeros((n,n))

for i in range(0,n):
   for j in range(0,n):
     interaction_index[i][j]=data[bead_type_index[i]][bead_type_index[j]]
     epsilon[i][j]=epslist[interaction_index[i][j]]
     sigma[i][j]=sigmalist[interaction_index[i][j]]
     vol_mar[i][j]=(sigma[i][j]/2)**3*4/3*numpy.pi

print epsilon
print vol_mar

for i in range(0,n):
   for j in range(0,n):
      epsilon_th[i][j]=epsilon[i][j]/(vol_mar[i][i]*vol_mar[j][j])-epsilon[i][0]/(vol_mar[0][0]*vol_mar[j][j])-epsilon[0][j]/(vol_mar[i][i]*vol_mar[0][0])+epsilon[0][0]/(vol_mar[0][0]**2)
      epsilon_th[i][j]=epsilon_th[i][j]*vol_th[i]*vol_th[j]

#epsilon=epsilon/epsilon[n-1][n-1]*3.0

print epsilon_th

f = open("epsilon.in","w")

for i in range(1,n):
    for j in range(1,n):
       f.write("%s " % epsilon_th[i][j])
    f.write("\n")


f.write(name)




