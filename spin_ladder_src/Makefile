#For Cineca's Marconi cluster:
#For GNU: module load gnu; module load profile/advanced; module load blas/3.6.0--gnu--6.1.0 lapack/3.6.1--gnu--6.1.0
#For Intel: module load intel/pe-xe-2018--binary; module load blas/3.6.0--intel--pe-xe-2016--binary lapack/3.8.0--intel--pe-xe-2018--binary RUNTIME REQUIRES INTEL MODULE
#For Cineca's Galileo cluster:
#For Intel: module load intel/pe-xe-2018--binary; module load blas/3.8.0--intel--pe-xe-2016--binary lapack/3.8.0--intel--pe-xe-2018--binary RUNTIME REQUIRES INTEL MODULE
#FC=gfortran
FC=ifort

OPT_DEF = -O2 -Wall
OPT_HIGH = -O3  -unroll -Wall
OPT_HIGH_INTEL = -O3  -unroll -warn all
OPT_HIGH_KNL = -O3  -unroll -xMIC-AVX512
OPT_DEBUG = -g  -O0 -Wall -fbounds-check -frange-check

FFLAGS =  $(OPT_HIGH_INTEL)


LIBS = -llapack -lblas
#For Cineca's Marconi cluster:
#LIBS = /cineca/prod/opt/libraries/lapack/3.6.1/gnu--6.1.0/lib/liblapack.so /cineca/prod/opt/libraries/blas/3.6.0/gnu--6.1.0/lib/libblas.so
#LIBS = /cineca/prod/opt/libraries/lapack/3.8.0/intel--pe-xe-2018--binary/lib/liblapack.so /cineca/prod/opt/libraries/blas/3.8.0/intel--pe-xe-2018--binary/lib/libblas.so
#For Cineca's Galileo cluster:
LIBS = /cineca/prod/opt/libraries/lapack/3.8.0/intel--pe-xe-2018--binary/lib/liblapack.so /cineca/prod/opt/libraries/blas/3.8.0/intel--pe-xe-2018--binary/lib/libblas.so

##########################################################

PREFIX=$(PWD)
MODDIR=$(PREFIX)/modules
SRCDIR=$(PREFIX)/src
OBJDIR=$(PREFIX)/objects
BINDIR=$(PREFIX)/bin

vpath %.x $(BINDIR)
vpath %.f90 $(SRCDIR)
vpath %.o $(OBJDIR)

#INC=-I$(MODDIR) -J$(MODDIR) 
INC=-I$(MODDIR)

##########################################################

SANDBOX = spin_ladder.x
MODULES = parameters.o input.o print_tools.o CG_tools.o J_SVD.o spin_h.o force.o


all: $(SANDBOX)

.SUFFIXES:
.SUFFIXES: .f90 .x .o

%.x : %.o $(OBJECTS)
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(INC) $(OBJDIR)/$< $(addprefix $(OBJDIR)/,$(MODULES)) $(LIBS)

%.o : %.f90
	$(FC) $(FFLAGS) $(INC) -c $(SRCDIR)/$*.f90 -o $(OBJDIR)/$*.o

#%.o : %.f90
#	$(FC) $(FFLAGS) -c $(INC) $< -o $(OBJDIR)/$@


clean:
	rm -f $(BINDIR)/$(SANDBOX) $(OBJDIR)/*.o $(MODDIR)/*.mod *~

spin_ladder.o: spin_ladder.f90 $(MODULES)
