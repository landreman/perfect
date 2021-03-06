# Common commands
# -------------------------------
# To build PERFECT:
#   make
# or, equivalently,
#   make all
#
# To run tests:
#   make test
#
# To run tests on perfectOutput.h5 for all examples that would be run for this system with "make test", without actually re-running these examples:
#   make retest
#
# To delete .mod and .o files:
#   make clean

# First, make sure the environment variable PERFECT_SYSTEM is set appropriately:

ifndef PERFECT_SYSTEM
$(error You must set the environment variable PERFECT_SYSTEM to select which makefile to use. For example, to use makefiles/makefile.edison, set PERFECT_SYSTEM=edison. To set this environment variable, use the export or setenv command, depending on your shell)
endif

SYSTEM_MAKEFILE = makefiles/makefile.$(PERFECT_SYSTEM)
ifeq ($(wildcard $(SYSTEM_MAKEFILE)),)
$(error Based on the environment variable PERFECT_SYSTEM, make is trying to include the file $(SYSTEM_MAKEFILE), but this file does not exist)
endif

TARGET = perfect

.PHONY: all clean test shorttest retest perfectVersion.h

all: $(TARGET)

include $(SYSTEM_MAKEFILE)

include makefile.depend

# Export variables that will be used in the quadpack subdirectory or for "make test":
export

%.o: %.F90
	$(FC) $(EXTRA_COMPILE_FLAGS) -c $<

quadpack/quadpack.a:
	$(MAKE) -C quadpack

perfectVersion.h:
	echo "#ifndef _PERFECTVERSION_H" > perfectVersion.h
	echo "#define _PERFECTVERSION_H" >> perfectVersion.h
	echo "#define GIT_COMMIT \"$$(git describe --always --dirty --abbrev=40)\"" >> perfectVersion.h
	echo "#endif" >> perfectVersion.h

$(TARGET): $(OBJ_FILES)
	$(FLINKER) -o $(TARGET) $(OBJ_FILES) $(EXTRA_LINK_FLAGS) $(PETSC_LIB)

# Since the sometimes-included PETSc makefile uses a double-colon to define clean, we must do the same here:
clean::
	rm -f *.o *.mod *.a *.dat *~ $(TARGET) perfectMatrices.m perfectOutput.h5
	cd quadpack; rm -f *.o *.a

test_make:
	@echo PERFECT_SYSTEM is $(PERFECT_SYSTEM)
	@echo FC is $(FC)
	@echo FLINKER is $(FLINKER)
	@echo EXTRA_COMPILE_FLAGS is $(EXTRA_COMPILE_FLAGS)
	@echo EXTRA_LINK_FLAGS is $(EXTRA_LINK_FLAGS)
	@echo PETSC_LIB is $(PETSC_LIB)
	@echo TARGET is $(TARGET)
	@echo OBJ_FILES is $(OBJ_FILES)

test: $(TARGET)
	@echo "Beginning functional tests." && cd examples && export PERFECT_RETEST=no && export PERFECT_RUN_LARGE_EXAMPLES=yes && ./runExamples.py

shorttest: $(TARGET)
	@echo "Beginning short functional tests." && cd examples && export PERFECT_RETEST=no && ./runExamples.py

retest: $(TARGET)
	@echo "Testing perfectOutput.h5 for examples without re-running then." && cd examples && export PERFECT_RETEST=yes && ./runExamples.py
