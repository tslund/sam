      subroutine assign_run_time( values_h )

      include 'sam.h'
      include 'mpif.h'

      real values_h(L_params)

c   *** Root broadcasts the run_time parameters to everyone else

      n_run_time = n_dynpar_r + n_dynpar_i

      call mpi_bcast( values_h(n_inputs+1), n_run_time,
     &     mpi_double_precision, i_root, mpi_comm_world, ierr )

c   *** Assign values to the /run_time/ common block

      dt_input = dt

      call assign_params( time, nt, values_h(n_inputs+1), 
     &                    n_dynpar_r, n_dynpar_i )

      t_start  = time
      nt_start = nt
      dt = dt_input
      n_frame_p = 0

      return
      end
