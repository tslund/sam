      subroutine read_header_p( fname, values_h )

      include 'sam.h'
      include 'mpif.h'

      character(*) fname
      real         values_h(L_params)

c   *** Only root reads the parameter values

      if(l_root) then
         call read_header( fname, values_h, ierr )
      end if
      call stop_on_error( ierr, 1 )

c   *** Root broadcasts the run_time parameters to everyone else

      n_run_time = n_dynpar_r + n_dynpar_i

      call mpi_bcast( values_h(n_inputs+1), n_run_time,
     &     mpi_double_precision, i_root, mpi_comm_world, ierr )

c   *** Assign values to the /run_time/ common block.  Root already did
c   *** this in read_header.

      if( .not. l_root ) then
         call assign_run_time( values_h )
      end if

      return
      end
