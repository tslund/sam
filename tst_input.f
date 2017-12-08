      program tst_input

      include 'sam.h'
      include 'mpif.h'

      real work(500)

      character(12) labels(L_params)
      real          values(L_params)
      logical        fixed(L_params), required(L_params), write_params

      ierr = 0
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

      if(l_root) then
      call input( 'input.dat', labels, values, required, fixed,
     &             write_params, ierr )
      end if
      call stop_on_error( ierr, 1 )
      call assign_static( values )

c   *** Determine index ranges for each process.

      call set_range( Nx, Nz )

      if(l_root) then
         call read_header( 'header.in', labels, values, fixed,
     &                     work(1), work(n_params+1), ierr )
      end if
      call stop_on_error( ierr, 1 )
      call assign_run_time( work )

      if(l_root) then
         call write_header( 'header.out', labels, values )
      end if

      call mpi_finalize( ierr )

      stop
      end
