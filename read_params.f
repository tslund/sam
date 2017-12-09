      subroutine read_params( fname, labels, values, found, n_vars,
     &                        ierr )

      real          values(n_vars)
      character(* ) fname
      character(12) labels(n_vars), label_lower, label
      logical        found(n_vars), done

      found = .false.
      done  = .false.

      open(newunit=i_unit,file=fname)

      do while( .not. done )
         call read_line( i_unit, label, value, done, ierr )
         if( ierr .ne. 0 ) return
         do i=1, n_vars
            call to_lower( labels(i), label_lower, 12 )
            if( label .eq. label_lower ) then
               found( i) = .true.
               values(i) = value
               goto 20
            end if
         end do
         print *, 'Warning: extraneous variable ', trim(label), 
     &            ' found in file ', trim(fname)
20       continue
      end do

      close(i_unit)

      return
      end 
