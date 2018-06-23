
         Spectral Atmospheric Model  (SAM)
           
Fourier-spectral Navier stokes solver using pressure projection.

To build the code do the following:
 a) Move into the path_to_sam/sam/src directory.
 b) Type configure to detect your compiler and to write the compiler.[spec,
    opt,debug] files.  If the configure command fails, make sure you have 
    perl installed and do in fact have a mpi/fortran environment set up.  
    If you use a compiler other than intel, gnu, pgi, or cray, you will
    need to edit the compiler.[spec,opt,debug] files to suit your compiler.
    Note that the makefile reads the compiler information from the file
    compiler.spec.  The files compiler.opt and compiler.debug are used to
    keep track of optimized and debug compiler flags.  If you need to debug
    the code then do 'make clean', 'cp compiler.debug compiler.spec',
    'make cgcam'.  Once you are done debugging do 'make clean', 
    'cp compiler.opt compiler.spec', 'make cgcam'.
 c) make sam
 d) make planes2vtk

There are a couple of test cases in the examples directory.  Look at the 
Readme files contained therein for instructions on how to run those cases.
