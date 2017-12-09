      subroutine to_lower( str_in, str_out, len_str )

      character(len_str) str_in, str_out, lower

c      str_out = lower(str_in)  ! this works instead of the lines blow

      str_out = str_in

      do i=1, len_str
         ic = ichar(str_in(i:i))                                             
         if( ic .ge. 65 .and. ic .le. 90 ) then
            str_out(i:i) = achar(ic+32)
         end if
      end do

      return
      end


      character(*) function lower( string )

      character(*) string

      len_string = len(string)

      lower(1:len_string) = string(1:len_string)

      do i=1, len_string
         ic = ichar(string(i:i))
         if( ic .ge. 65 .and. ic .le. 90 ) then
            lower(i:i) = achar(ic+32)
         end if
      end do

      return
      end
