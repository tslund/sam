      subroutine grab_input( i_unit, var_name, value, ivalue, ierr )

      character(*)  var_name
      character(len(var_name)) var_name_lower
      character(12) label
      logical done

      ierr = 0
      done = .false.

      rewind(i_unit)

c   *** read_line already converts labels to lower case.  We convert
c   *** var_name to lower case here in line in order to compare strictly
c   *** lower case strings.

      call to_lower( var_name, var_name_lower, len(var_name) )

      do while( .not. done )
         call read_line( i_unit, label, value, done, ierr )
         if( label .eq. var_name_lower ) then
            ivalue = nint(value)
            return
         end if
      end do

      print *, 'Error in grab_input: unable to find ', var_name
      ierr = 1

      return
      end
