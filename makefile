TARGET = assembly

#SRC = modules.f90 SPmain.f90 parser.f90 init.f90 allocation.f90 allocateell.f90 3D.f90 cadenas.f90 cadenasMK.f90 fe.f90  fkfun.f90  kai.f90  kinsol.f90  pxs.f90  savetodisk.f90 rands.f90 ellipsoid.f90 dielectric.f90 monomers.definitions-onck.f90 chains.definitions.f90 sphere.f90 kapfromfile.f90

SRC = modules.f90 read.f90 allocation.f90 main.f90 fkfun.f90 cadenas.f90 rands.f90 kinsol.f90 initmpi.f90 kai.f90 factorcurv.f90 free_energy.f90

HOST=$(shell hostname)
$(info HOST is ${HOST})


# some definitions
SHELL = /bin/bash
FFLAGS= -O3# -g -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow # -O3

ifeq ($(HOST),skay)
LFLAGS = -lm /usr/lib/x86_64-linux-gnu/librt.so  -L/usr/local/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/usr/local/lib
endif



ifeq ($(HOST),spinetta) 
LFLAGS = -lm /usr/lib/x86_64-linux-gnu/librt.so  -L/usr/local/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/usr/local/lib
endif



ifeq ($(HOST),master) 


LFLAGS = -L/shared/software/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath

endif
ifeq ($(HOST),mate.bme.northwestern.edu) 


LFLAGS = -L/home/mario/software/kinsol/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath

endif

ifeq ($(HOST),quser13) 

LFLAGS = -L/home/mta183/KINSOL2.7/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -lgfortranbegin -lgfortran -lm
endif


ifeq ($(HOST),quser12) 

LFLAGS = -L/home/mta183/KINSOL2.7/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -lgfortranbegin -lgfortran -lm
endif

ifeq ($(HOST),quser11)

LFLAGS = -L/home/mta183/KINSOL2.7/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -lgfortranbegin -lgfortran -lm
endif

ifeq ($(HOST),quser10)

LFLAGS = -L/home/mta183/KINSOL2.7/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -lgfortranbegin -lgfortran -lm
endif

GIT_VERSION := $(shell git describe --abbrev=6 --dirty --always --tags)
GFLAGS=-cpp -D_VERSION=\"$(GIT_VERSION)\"

FF = mpif77 #${F90}
VER = ~/bin/copolymer

all:	$(TARGET)

$(TARGET): $(SRC:.f90=.o)
	$(FF) -o $(TARGET) $(SRC:.f90=.o) $(LFLAGS) $(GFLAGS)
	cp $(TARGET) $(VER)

$(SRC:.f90=.o): $(SRC)
	${FF} -c ${FFLAGS}  $(SRC) $(LFLAGS) $(GFLAGS)

install: all
	cp $(TARGET) $(VER)

clean:	
	@rm -f $(SRC:.f90=.o) $(SRC:.f90=.d) $(TARGET) *~

realclean: clean
	@rm -f .depend

depend dep:
	@$(FF)  $(CFLAGS) -MM $(SRC) > .depend 

ifeq (.depend, $(wildcard .depend))
include .depend
endif
















































