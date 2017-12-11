      subroutine input_p( fname, write_params )

      include 'sam.h'
      include 'mpif.h'

      character(*) fname
      integer      isend(20)
      logical      write_params

c   *** Only root reads the input file.  Everyone else just sets the
c   *** labels.

      if(l_root) then
         call input( fname, write_params, ierr )
      else
         call set_labels( labels, values, required, fixed, L_params, 
     &                    n_inputs_r, n_inputs_i, n_inputs,
     &                    n_dynpar_r, n_dynpar_i, n_dynpar, n_params,
     &                    ierr )
      end if
      call stop_on_error( ierr, 1 )

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

c   *** Assign the input data to the /params/ common block (root has
c   *** already done this in input).

      if( .not. l_root ) then
         call assign_params( loc(xL), loc(Nx), values,
     &                       n_inputs_r, n_inputs_i )
      end if

      return
      end
