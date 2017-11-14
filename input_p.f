      subroutine input_p( fname, labels, values, fixed, write_params )

      include 'sam.h'
      include 'mpif.h'

      character( *) fname
      character( 8) default_flag
      character(12) labels(*)
      real          values(*)
      logical        fixed(*), write_params
      logical     required(L_params), found(L_params)
      pointer(ipr, rdat)
      pointer(ipi, idat)
      real    rdat(L_params)
      integer idat(L_params), isend(5)

      call set_labels( labels, values, required, fixed )

c   *** Only root reads the input parameter values

      if( l_root ) then
         call read_params(fname, labels, values, found, n_inputs, ierr)
         values(n_inputs+1) = ierr
      end if   

c   *** Root broadcasts the parameter vaues to everyone else, including
c   *** the error status from the call to read_params (which is sent in
c   *** values(n_inputs+1))

      call mpi_bcast( values, n_inputs+1, mpi_double_precision,
     &                i_root, mpi_comm_world, ierr )

c   *** Halt if an error was encountered in the read_params routine

      ierr = values(n_inputs+1)
      if( ierr .ne. 0 ) goto 99

c   *** Assign the input data to the /params/ common block

      ipr = loc(xL)
      ipi = loc(Nx)

      n_inputs_i = n_inputs - n_inputs_r

      rdat(1:n_inputs_r) = values(1:n_inputs_r)
      idat(1:n_inputs_i) = values(n_inputs_r+1:n_inputs)

c   *** Root checks for errors in the input values

      if( l_root ) then
         call check_inputs( labels, values, found, required, ierr )
         isend(1) = ierr
         isend(2) = n_xy_planes
         isend(3) = n_xz_planes
         isend(4) = n_yz_planes
         isend(5) = i_gw_type
      end if

      call mpi_bcast( isend, 5, mpi_integer, i_root, mpi_comm_world,
     &                ierr_mpi )

      ierr        = isend(1)
      n_xy_planes = isend(2)
      n_xz_planes = isend(3)
      n_yz_planes = isend(4)
      i_gw_type   = isend(5)

      if( ierr .ne. 0 ) goto 99

c   *** Root writes the input parameters to file params.out

      if( l_root .and. write_params ) then

         open(newunit=i_unit,file='params.out')

         do i=1, n_inputs
            default_flag = ' '
            if( .not. found(i) ) default_flag = 'DEFAULT'
            write(i_unit,*) labels(i), values(i), default_flag
         end do

         close(i_unit)

      end if

      return

99    call mpi_finalize( ierr )
      stop

      end
