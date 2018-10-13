#  © 2002 – 2014 NorthWest Research Associates, Inc. All Rights Reserved
#    Author: Thomas S. Lund, lund@cora.nwra.com

include compiler.spec

include dependencies.make

#-------------------------------------------------------------------------------
#   Defile lists of mpi and serial programs

MPI_PROGRAMS    = sam             \
                  tst_input       \
                  tst_transpose   \

SERIAL_PROGRAMS = check           \
                  compare_vel     \
                  planes2vtk      \
                  remesh          \
                  tst_pointer     \


#-------------------------------------------------------------------------------
#   Define lists of mpi and serial object files

MPI_OBJECTS    := $(shell grep -l 'mpif\.h' *.f | sed 's:\.f:\.o:')
SERIAL_OBJECTS := $(shell grep -L 'mpif\.h' *.f | sed 's:\.f:\.o:')


#-------------------------------------------------------------------------------
#   Rule for syncing the header file with set_labels.f

header: 
	syncHeader

cgcam.h: set_labels.f
	syncHeader


#-------------------------------------------------------------------------------
#   Rule for generating dependency information

depend: 
	get_dependencies $(MPI_PROGRAMS) $(SERIAL_PROGRAMS)


#-------------------------------------------------------------------------------
#   Rule for doing some gentle cleaning

clean: 
	rm *.o


#-------------------------------------------------------------------------------
#   Rule for doing some more serious cleaning

clean_all: 
	rm *.o $(SERIAL_PROGRAMS) $(MPI_PROGRAMS)


#-------------------------------------------------------------------------------
#   Rules for building both mpi and serial programs

$(MPI_PROGRAMS): %:
	$(MPIFC) $(FFLAGS) $^ -o $@

$(SERIAL_PROGRAMS): %:
	$(FC)    $(FFLAGS) $^ -o $@


#-------------------------------------------------------------------------------
#   Rules for compiling both mpi and serial object files

$(MPI_OBJECTS): %.o: %.f
	$(MPIFC) $(FFLAGS) -c $<

$(SERIAL_OBJECTS): %.o: %.f
	$(FC) $(FFLAGS) -c $<
