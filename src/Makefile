# GNU Makefile
#
.SUFFIXES: .o .f90
FC=gfortran
FCFLAGS= -Ofast -march=native -funsafe-math-optimizations
#FCFLAGS= -O3
#FCFLAGS= -O -g -fbounds-check -Wuninitialized -Wsurprising -Wall -Wextra
#FCFLAGS += -I/usr/include
FOMP= -fopenmp

PROGRAM = ../Ndyn

## "make" builds all
#all: $(PROGRAM)

# list all source files
OBJ = mods.o domtype.o read_input.o potentials.o nd.o

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
${PROGRAM}:${OBJ}
	$(FC) -o ${PROGRAM} ${FOMP} ${OBJ} ${FCFLAGS}

# to clean up stuff
clean:
	rm -f *.o *.mod *.MOD
veryclean: clean
	rm -f *~ $(PROGRAMS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
nd.o:nd.f90
	$(FC) $(FCFLAGS) $(FOMP) -c nd.f90
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<


