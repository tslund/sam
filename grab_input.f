      subroutine grab_input( i_unit, var_name, value, ivalue, ierr )

      character(*) var_name
      character(12) label, var_name1
      logical done

      ierr = 0
      done = .false.

c   *** Make sure that var_name is all lower case.  We need to make a
c   *** copy first since the input value may not have a memory address

      var_name1 = var_name
      call to_lower( var_name1 )

      rewind(i_unit)

      do while( .not. done )
         call read_line( i_unit, label, value, ierr, done )
         if( trim(label) .eq. trim(var_name1) ) then
            ivalue = nint(value)
            return
         end if
      end do

      print *, 'Error in grab_input: unable to find ', var_name
      ierr = 1

      return
      end
