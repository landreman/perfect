# Makefile for PERFECT

MATTS_INCLUDE =
MATTS_LIBS =

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
#FLINKER = ftn

TARGET = perfect

OBJS = $(TARGET).o xGrid.o geometry.o writeHDF5Output.o globalVariables.o scan.o readInput.o \
	polynomialDiffMatrices.o uniformDiffMatrices.o \
	uniformInterpolationMatrix.o printToStdout.o \
	profiles.o sparsify.o \
	polynomialInterpolationMatrix.o solveDKE.o quadpack/quadpack.a 

# Dependencies:
geometry.o: globalVariables.o
writeHDF5Output.o: globalVariables.o scan.o
printToStdout.o: globalVariables.o
readInput.o: globalVariables.o
solveDKE.o: globalVariables.o xGrid.o geometry.o polynomialDiffMatrices.o printToStdout.o profiles.o sparsify.o
scan.o: globalVariables.o
$(TARGET).o: globalVariables.o geometry.o writeHDF5Output.o printToStdout.o solveDKE.o readInput.o scan.o

%.o: %.F90
	$(FLINKER) $(MATTS_INCLUDE) -c $<

#erf.o: erf.c
#	$(CC) -o erf.o -c erf.c -I$(GSL_DIR)/include

quadpack/quadpack.a:
	cd quadpack; make all

all: $(OBJS)
	$(FLINKER) -o $(TARGET) $(OBJS) $(MATTS_LIBS) $(PETSC_LIB)

cleaner:
	rm -f *.o *.mod *.a *.dat *~ $(TARGET) *.m *.h5
	cd quadpack; rm -f *.o *.a

