      subroutine input( fname, write_params, ierr )

      include 'sam.h'

      character( *) fname
      character( 8) default_flag
      logical       write_params

      call set_labels( labels, values, required, fixed, L_params, 
     &                 n_inputs_r, n_inputs_i, n_inputs,
     &                 n_dynpar_r, n_dynpar_i, n_dynpar, n_params )

      if( n_params .gt. L_params ) then
         print *, 'ERROR: The number of input parameters exceeds ',
     &             the value of L_params'
         ierr = 1
         return
      end if

      call read_params(fname, labels, values, found, n_inputs, ierr)
      if( ierr .ne. 0 ) return

c   *** Assign input data to the /params/ common block

      call assign_params( loc(xL), loc(Nx), values,
     &                    n_inputs_r, n_inputs_i )

      call check_inputs( ierr )
      if( ierr .ne. 0 ) return

c   *** Write the input parameters to file params.out

      if( write_params ) then

         open(newunit=i_unit,file='params.out')

         do i=1, n_inputs
            default_flag = ' '
            if( .not. found(i) ) default_flag = 'DEFAULT'
            write(i_unit,*) labels(i), values(i), default_flag
         end do

         close(i_unit)

      end if

      return
      end
