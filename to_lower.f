      subroutine to_lower( string_in )

      character(*) string_in

      n = len(string_in)

      do i=1, n
         ic = ichar(string_in(i:i))
         if( ic .ge. 65 .and. ic .le. 90 ) then
            string_in(i:i) = achar(ic+32)
         end if
      end do

      return
      end
