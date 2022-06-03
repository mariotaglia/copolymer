TARGET = assembly

#SRC = modules.f90 SPmain.f90 parser.f90 init.f90 allocation.f90 allocateell.f90 3D.f90 cadenas.f90 cadenasMK.f90 fe.f90  fkfun.f90  kai.f90  kinsol.f90  pxs.f90  savetodisk.f90 rands.f90 ellipsoid.f90 dielectric.f90 monomers.definitions-onck.f90 chains.definitions.f90 sphere.f90 kapfromfile.f90

SRC = modules.f90 allocation.f90 main.f90 fkfun.f90 cadenas.f90 rands.f90 initmpi.f90 kai.f90 factorcurv.f90 free_energy.f90  dielectric.f90 parser.f90 PBC.f90 creador.f90 init.f90 save2disk.f90 kinsol.f90

HOST=$(shell hostname)
$(info HOST is ${HOST})

USER=$(shell id -u -n)
$(info USER is ${USER})


# some definitions
SHELL = /bin/bash

# FFLAGS= -Wall 
FLAGS= -O3 -fbounds-check
# FFLAGS=-g -fbacktrace -fbounds-check # -ffpe-trap=zero,overflow,underflow 

LFLAGS = -lm /usr/lib/x86_64-linux-gnu/librt.so  -L/usr/local/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/usr/local/lib


ifeq ($(HOST),login.tusker.hcc.unl.edu)
LFLAGS = -L/home/conda/gzaldivar/bin/kinsol/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/home/conda/gzaldivar/bin/kinsol/lib
endif


ifeq ($(HOST),mdq)
LFLAGS = -lm /usr/lib/x86_64-linux-gnu/librt.so  -L/usr/local/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/usr/local/lib
endif


ifeq ($(HOST),login.crane.hcc.unl.edu)
ifeq ($(USER),gzaldivar)
LFLAGS = -L/home/conda/gzaldivar/bin/kinsol/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/home/conda/gzaldivar/bin/kinsol/lib
endif
endif

ifeq ($(HOST),service0)
LFLAGS = -lm  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial
endif

ifeq ($(HOST),skay)
LFLAGS = -lm /usr/lib/x86_64-linux-gnu/librt.so  -L/usr/local/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/usr/local/lib
endif
ifeq ($(HOST),carapa)
LFLAGS = -lm /usr/lib/x86_64-linux-gnu/librt.so  -L/usr/local/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/usr/local/lib
endif


ifeq ($(HOST),odin)
LFLAGS = -lm /usr/lib64/librt.so -L/home/gzaldivar/software/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/home/gzaldivar/software/lib
endif

ifeq ($(HOST),logan.q1.fcen.uba.ar)
LFLAGS = -lm /usr/lib64/librt.so -L/home/gerva/bin/kinsol/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/home/gerva/bin/kinsol/lib
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

FF = mpif90 #${F90}
VER = ~/bin/copolymer

all:	$(TARGET)

$(TARGET): $(SRC:.f90=.o)
	$(FF) -o $(TARGET) $(SRC:.f90=.o) $(LFLAGS) 
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
















































