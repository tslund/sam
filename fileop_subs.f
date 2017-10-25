      subroutine get_fsize( fname, word_size, fsize, n_rec )

c  © 2002 – 2014 NorthWest Research Associates, Inc. All Rights Reserved
c    Author: Thomas S. Lund, lund@cora.nwra.com

c   *** Determine the size of a file in bytes (up to 1.0 Tb max size).

c   *** One should be able to use the inquire intrinsic to determine the
c   *** file size.  I tried this, but could not get it to work properly
c   *** on all platforms where I needed it.  This routine is somewhat
c   *** complicated, but at least it is portable.

c   *** Earlier versions of this routine made use of the error code for a
c   *** read past the end of file.  Unfortunately this code is variable
c   *** across different compilers and/or across different platforms and
c   *** thus an attempt was made to discover the error code by reading 
c   *** something like 1 TB into the file (assuming that the file was not 
c   *** this large and thus would give an end of file error).
c   *** This approach does not work (on the intel compiler at least) since
c   *** different error codes are returned for end of file less that 2 GB
c   *** in size (code 36) and greater than 2 GB in size (error code 25).
c   *** The work around at the moment is to assume that any error code other 
c   *** than 0 signals a read past the end of file.

c   *** Taking the word_size to be > 8 can work, but may lead to
c   *** erroneous results.  Not sure why, but iostat errors do not arise
c   *** consistently when a read is attempted past the end of file when
c   *** the record size is large.  It is best to have the word_size match
c   *** the type of data being read.  This appears to work reliably in all
c   *** cases.

      implicit none
      character(*) fname
      integer iunit, ios, word_size, eof_code, n
      byte, allocatable, dimension(:) :: t
      integer(8) max_fsize, max_rec, start, end, i, fsize, n_rec

c      if( word_size .gt. 8 ) then
c         print *, 'WARNING: In get_fize word_size = ', word_size,
c     &            ' which is greater than 8'
c         print *, 'This may be ok, but check to see if this is what ',
c     &            'was intended'
c      end if

      allocate( t(word_size) )

      max_fsize = 1.0e+12   ! 1.0 Tb
      max_rec = max_fsize/word_size

      open( newunit=iunit, file=fname, status='old', access='direct',
     &      form='unformatted', recl=word_size )

c   *** Make sure that we can read the file.  Exit with an error if not.

      read( iunit, rec=1,      iostat=ios      ) t
      if( ios .ne. 0 ) then
         fsize = -1
         n_rec = -1
         goto 99
      end if

c   *** Attempt to discover the end of file exit code.
c   *** We do not use this any longer since the error code returned can depend
c   *** on the file size.

c      read( iunit, rec=max_rec, iostat=eof_code ) t
c
cc      print *, 'eof_code = ', eof_code
c
c      if( eof_code .eq. 0 ) then
c         fsize = -99
c         n_rec = -99
c         goto 99
c      end if

c   *** Check to see if the file size is equal to the input value of fsize.
c   *** If so, we jump to the end thereby bypassing the binary search.

      if( fsize .gt. 0 ) then
         n_rec = fsize/word_size
         read( iunit, rec=n_rec,   iostat=ios ) t
         if( ios .ne. 0 ) goto 10
         read( iunit, rec=n_rec+1, iostat=ios ) t
         if( ios .ne. 0 ) goto 99
      end if

10    continue
         
c   *** Binary search for file end.  if we get iostat=0, we're before the 
c   *** end of file, if we get an iostat!=0, we're past the end of file.

      start = 0
      end = max_rec

      do n=1, 100
      
         i = ( start + end )/2

         read(iunit, rec=i, iostat=ios ) t

c         print *, n, i, ios, i*word_size

         if( ios .eq. 0 ) then
            start = i
c         else if( ios .eq. eof_code ) then
         else 
            end = i
c         else
c            fsize = -1
c            n_rec = -1
c            goto 99
         end if

c   *** The search is done if the interval is reduced to 1

         if( end-start .le. 1 ) then
            fsize = start*word_size
            n_rec = start
            goto 99
         end if

      end do

c   *** If control gets here the binary search failed

      fsize = -2
      n_rec = -2

99    close( iunit )

      deallocate( t )

      return
      end


      subroutine delete_file( fname )

      implicit none
      integer iunit
      character(*) fname
      logical file_exists

      inquire(file=fname,exist=file_exists)

      if( file_exists ) then
         open(newunit=iunit,file=fname)
         close(iunit,status='delete')
         print *, 'deleted file ', trim(fname)
      end if

      return
      end


      subroutine move_file( fname1, fname2 )

      implicit none
      character(*) fname1, fname2

      call copy_file( fname1, fname2, 1 )

      return
      end


      subroutine copy_file( fname1, fname2, i_delete )

      implicit none
      character(*)       :: fname1, fname2
      integer            :: i_delete
      integer, parameter :: ibs=131072  ! blocksize = 256*512 bytes
      integer            :: i, iunit1, iunit2, word_size, nbks,
     &                      i_left_over
      integer(8)         :: fsize, n_rec, n_bks
      character(1)       :: data(ibs)
      logical            :: file_exists

      inquire(file=fname1,exist=file_exists)

      if( .not. file_exists ) then
         print *, 'ERROR: can not copy non-existant file ', trim(fname1)
         stop
      end if

      word_size = 1
      call get_fsize( fname1, word_size, fsize, n_rec )
      n_bks = fsize/int(ibs,kind=8)
      i_left_over = fsize - int(n_bks*ibs,kind=8)

      call delete_file( fname2 )

      open(newunit=iunit1,file=fname1,form='unformatted',
     &     access='stream',action= 'read')
      open(newunit=iunit2,file=fname2,form='unformatted',
     &     access='stream',action='write')

      do i=1, n_bks
         read( iunit1) data
         write(iunit2) data
      end do

      read( iunit1) data(1:i_left_over)
      write(iunit2) data(1:i_left_over)

      print *, 'copied file ', trim(fname1), ' to ', trim(fname2)

      if( i_delete .eq. 1 ) then
         close(iunit1,status='delete')
         print *, 'deleted file ', trim(fname1)
      else
         close(iunit1)
      end if

      close(iunit2)

      return
      end
