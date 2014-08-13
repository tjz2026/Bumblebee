.SUFFIXES: .f90

include ./make.sys

default: all

modm = global-mpi_mod.o SCMFT_constants.o
modn = utility_mod.o 
mpiu = control_mod.o mod_mpi.o 
dmft = cosgmma_mod.o initialize_mod.o SCMFT_dump_mod.o
core = scmft_mod.o MCstatistics_mod.o order_para_mod.o
main = SCMFTbrush-slab.o

mkllib=/root/intel_nfs/composer_xe_2013.3.163/mkl/lib/intel64
mklinc=/root/intel_nfs/composer_xe_2013.3.163/mkl/include
MKLFCCFLAG= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -liomp5 -lpthread
objects = $(modm) $(modn) $(mpiu) $(dmft) $(core) $(main)

all: SCMFTbrush-slab

SCMFTbrush-slab: $(objects)
	$(LINKER) $(objects) -o SCMFTbrush-slab $(LFLAGS) $(LIBS) -I$(mklinc) -L$(mkllib) $(MKLFCCFLAG)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f SCMFTbrush-slab

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
