      subroutine grab_input( i_unit, var_name, value, ivalue, ierr )

      character(*) var_name
      character(12) label
      logical done

      ierr = 0

      done = .false.

      rewind(i_unit)

      do while( .not. done )
         call read_line( i_unit, label, value, ierr, done )
         if( trim(label) .eq. trim(var_name) ) then
            ivalue = nint(value)
            return
         end if
      end do

      print *, 'Error in grab_input: unable to find ', var_name
      ierr = 1

      return
      end
