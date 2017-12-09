      subroutine read_line( i_unit, label, value, done, ierr )

c  © 2002 – 2014 NorthWest Research Associates, Inc. All Rights Reserved
c    Author: Thomas S. Lund, lund@cora.nwra.com

c   *** Error codes:
c          0 - normal exit
c          2 - non-numeric character in numeric field
c          4 - invalic character in label
c          5 - missing numeric entry

      character*12 label
      character*80 line
      logical done

      ierr = 0
      done = .false.

10    continue

c   *** Read a line and extract the label.  Comments are flagged with #
c   *** Either a blank line or the string 'eof' signifies the end of file

      read(i_unit,'(a80)',end=50) line

      if( line(1:1) .eq. '#' ) goto 10

      do i=1, 80
         if( line(i:i) .ne. ' ' ) goto 20
      end do
c   *** If flow gets here we found a blank line, which implies end of file
      goto 50

20    i_lab_beg = i

      i_lab_end = index(line(i_lab_beg+1:80),' ')-1 + i_lab_beg
      len_label = i_lab_end-i_lab_beg+1
      label = line(i_lab_beg:i_lab_end)

c   *** Check for non-alpha characters and convert to lower case

      do i=1, len_label
         ic = ichar(label(i:i))
         if(  ic .lt. 48                   .or.
     &       (ic .gt. 57 .and. ic .lt. 65) .or. 
     &       (ic .gt. 90 .and. ic .lt. 95) .or.
     &        ic .eq. 96                   .or.
     &        ic .gt. 122                       ) then
            print *, 'ERROR: Invalid character in label = ', label
            ierr = 4
            return
         end if
         if( ic .ge. 65 .and. ic .le. 90 ) then
            label(i:i) = achar(ic+32)
         end if
      end do

      if( label(1:len_label) .eq. 'eof' ) goto 50

c   *** Isolate the numeric input and convert it from character to
c   *** floating point value.

      i_num_beg = scan( line(i_lab_end+1:80), '+-.0123456789' ) +
     &            i_lab_end
      if( i_num_beg .eq. i_lab_end ) then
         print *, 'ERROR: No numeric entry found'
         ierr = 5
         return
      end if

      i_num_end = index(line(i_num_beg+1:80),' ')-1 + i_num_beg

c   *** Look for 10^2, 10^3, and 10^6 multipliers appended to
c   *** numerical input

      mult = 1
      select case(line(i_num_end:i_num_end))
      case('c')
         mult = 100
      case('k')
         mult = 1000
      case('m')
         mult = 1000000
      end select
      if( mult .ne. 1 ) i_num_end = i_num_end - 1

      read( line(i_num_beg:i_num_end), *, err=40 ) value

      if( mult .ne. 1 ) value = value*float(mult)

      return

40    print *, 'ERROR: Malformed numeric entry'
      ierr = 2
      return

50    done = .true.
      return

      end
