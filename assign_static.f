      subroutine assign_static( values )

      include 'sam.h'
      include 'mpif.h'

      real    values(L_params)
      integer isend(20)

c   *** Root broadcasts the parameter values determined in check_inputs
c   *** to everyone else

      if(l_root) then
         isend( 1) = n_xy_planes
         isend( 2) = n_xz_planes
         isend( 3) = n_yz_planes
         isend( 4) = i_gw_type
      end if

      call mpi_bcast( isend, 4, mpi_integer, i_root,
     &                mpi_comm_world, ierr )

      n_xy_planes = isend( 1)
      n_xz_planes = isend( 2)
      n_yz_planes = isend( 3)
      i_gw_type   = isend( 4)

c   *** Root broadcasts the input parameter values to everyone else

      call mpi_bcast( values, n_inputs, mpi_double_precision, i_root,
     &                mpi_comm_world, ierr )

c   *** Assign the input data to the /params/ common block

      call assign_params( xL, Nx, values, n_inputs_r, n_inputs_i )

      return
      end
