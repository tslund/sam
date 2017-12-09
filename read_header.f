      subroutine read_header( fname, values_h, ierr )

      include 'sam.h'

      character( *) fname
      real          values_h(L_params)

      ierr = 0

      squak = .true.
      i = index_param( 'i_restart', labels, n_inputs )
      squak(i) = .false.
      i = index_param( 'nt_restart', labels, n_inputs )
      squak(i) = .false.

c   *** Read the parameter values

      call read_params( fname, labels, values_h, found, n_params,
     &                  ierr )
      if( ierr .ne. 0 ) return

c   *** Check for inconsistencies between the input and header
c   *** file values.

      do i=1, n_inputs
         if( abs(values_h(i)-values(i)) .gt. 1.0e-12 ) then
            if( fixed(i) ) then
               print *, 'ERROR: parameter ', trim(labels(i)),
     &         ' differs between the header and input file'
               print *, values_h(i), values(i), values_h(i)-values(i)
               ierr = 1
            else if( squak(i) ) then
               print *, 'WARNING: parameter ', trim(labels(i)),
     &         ' differs between the header and input file'
               print *, values_h(i), values(i)
            end if
         end if
      end do

c   *** Assign values to the /run_time/ common block

      call assign_run_time( values_h )

      return
      end
