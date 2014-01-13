# Makefile for PERFECT

# The following 2 lines are required to find the HDF5 libraries and include files on loki:
#MATTS_INCLUDE = -I/usr/local/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -OPT:Olimit=0
#MATTS_INCLUDE = -O 0
MATTS_INCLUDE =
MATTS_LIBS =
#MATTS_LIBS = /usr/local/lib/libhdf5_fortran.a /usr/local/lib/libhdf5.a /usr/local/lib/libhdf5hl_fortran.a /usr/local/lib/libhdf5_hl.a -lz -lsz
#MATTS_LIBS = -lhdf5_fortran -lhdf5 -lhdf5hl_fortran -lhdf5_hl
#MATTS_LIBS = -L/usr/local/lib -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran
# -L/usr/local/lib/libhdf_hl.a -L/usr/local/lib/libhdf_fortran.a

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
#FLINKER = ftn

TARGET = perfect

OBJS = $(TARGET).o xGrid.o geometry.o writeHDF5Output.o globalVariables.o scan.o readInput.o \
	polynomialDiffMatrices.o uniformDiffMatrices.o \
	uniformInterpolationMatrix.o printToStdout.o ChebyshevGrid.o \
	profiles.o sparsify.o \
	polynomialInterpolationMatrix.o solveDKE.o quadpack/quadpack.a 

# Dependencies:
geometry.o: globalVariables.o
writeHDF5Output.o: globalVariables.o scan.o
printToStdout.o: globalVariables.o
readInput.o: globalVariables.o
solveDKE.o: globalVariables.o xGrid.o geometry.o polynomialDiffMatrices.o printToStdout.o ChebyshevGrid.o profiles.o sparsify.o
scan.o: globalVariables.o
$(TARGET).o: globalVariables.o geometry.o writeHDF5Output.o printToStdout.o solveDKE.o readInput.o scan.o

%.o: %.F90
	$(FLINKER) $(MATTS_INCLUDE) -c $<

#erf.o: erf.c
#	$(CC) -o erf.o -c erf.c -I$(GSL_DIR)/include

quadpack/quadpack.a:
	cd quadpack; make all

all: $(OBJS)
#	$(FLINKER) -o $(TARGET) $(OBJS) $(MATTS_LIBS)
#	$(FLINKER) -o $(TARGET) $(OBJS) $(MATTS_LIBS) -L$(GSL_DIR)/lib -lgsl
#	$(FLINKER) -o $(TARGET) $(OBJS) -lm -lc  $(PETSC_LIB) $(MATTS_LIBS) -L$(GSL_DIR)/lib -lgsl
	$(FLINKER) -o $(TARGET) $(OBJS) -lm -lc $(MATTS_LIBS) $(PETSC_LIB)

cleaner: clean
#clean: 
	rm -f *.o *.mod *.a *.dat *~ $(TARGET) *.m *.h5
	cd quadpack; rm -f *.o *.a
