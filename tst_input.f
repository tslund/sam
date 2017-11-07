      program tst_input

      include 'sam.h'
      include 'mpif.h'

      real work(500)

      character(12) labels(L_params)
      real          values(L_params)
      logical        fixed(L_params), write_params

      write_params = .true.

c   *** Initialize MPI, get myid, numprocs, and test if on root process

      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myid,ierr)
      call mpi_comm_size(mpi_comm_world,numprocs,ierr)

      i_root = 0
      l_root = .false.
      if( myid .eq. i_root ) l_root = .true.

c   *** Write a copyright message

      if(l_root) then
         write(6,6)
6        format(/,'© 2017 NorthWest Research Associates, '
     &            'Inc.  All Rights Reserved',/,
     &            'Author: Thomas S. Lund, lund@cora.nwra.com',/)
      end if

      call input_p( 'input.dat', labels, values, fixed, write_params )

c   *** Determine index ranges for each process.

      call set_range( Nx, Nz )

      call read_header( 'header.in', labels, values, fixed,
     &                  work(1), work(n_params+1) )

      if(l_root) then
         call write_header( 'header.out', labels, values )
      end if

      call mpi_finalize( ierr )

      stop
      end
