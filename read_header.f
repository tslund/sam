      subroutine read_header( fname, labels, values, fixed, 
     &                        values_h, found )

      include 'sam.h'
      include 'mpif.h'

      character( *) fname
      character(12) labels(L_params)
      real          values(L_params), values_h(L_params)
      logical        fixed(L_params),    found(L_params)
      pointer(ipr, rdat)
      pointer(ipi, idat)
      real    rdat(L_params)
      integer idat(L_params)

c   *** Only root reads the parameter values

      if( l_root ) then
         call read_params( fname, labels, values_h, found, n_params,
     &                     ierr )
         values_h(n_params+1) = ierr
      end if

c   *** Root broadcasts the parameter vaues to everyone else, including
c   *** the error status from the call to read_params (which is sent in
c   *** values(n_params+1))

      call mpi_bcast( values, n_params+1, mpi_double_precision,
     &                i_root, mpi_comm_world, ierr )

c   *** Halt if an error was encountered in the read_params routine

      ierr = values(n_params+1)
      if( ierr .ne. 0 ) goto 99

c   *** Root checks for inconsistencies between the input and header
c   *** file values.

      if( l_root ) then
         do i=1, n_inputs
            if( values_h(i) .ne. values(i) ) then
               if( fixed(i) ) then
                  print *, 'ERROR: parameter ', trim(labels(i)),
     &            ' differs between the header and input file'
                  print *, values_h(i), values(i)
                  ierr = 1
               else
                  print *, 'WARNING: parameter ', trim(labels(i)),
     &            ' differs between the header and input file'
                  print *, values_h(i), values(i)
               end if
            end if
         end do
      end if

c   *** Halt if an error was encountered while checking the parameter values

      call mpi_bcast( ierr, 1, mpi_integer, i_root, mpi_comm_world,
     &                ierr_mpi )

      if( ierr .ne. 0 ) goto 99

c   *** Assign values to the /run_time/ common block

      dt_input = dt

      ipr = loc(time)
      ipi = loc(nt)

      i_end = n_inputs
      i_beg=i_end+1;  i_dat=n_dynpar_r;  i_end=i_beg+i_dat-1
         rdat(1:i_dat) = values_h(i_beg:i_end)
      i_beg=i_end+1;  i_dat=n_dynpar_i;  i_end=i_beg+i_dat-1
         idat(1:i_dat) = values_h(i_beg:i_end)

      t_start  = time
      nt_start = nt
      dt = dt_input

      return

99    call mpi_finalize( ierr )
      stop

      end
