      subroutine read_header( fname, labels, values, fixed, 
     &                        values_h, found, ierr )

      include 'sam.h'

      character( *) fname
      character(12) labels(L_params)
      real          values(L_params), values_h(L_params)
      logical        fixed(L_params),    found(L_params),
     &               squak(L_params)

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

      dt_input = dt

      call assign_params( time, nt, values_h(n_inputs+1), 
     &                    n_dynpar_r, n_dynpar_i )

      t_start  = time
      nt_start = nt
      dt = dt_input

      return
      end
