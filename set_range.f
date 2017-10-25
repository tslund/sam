      subroutine set_range( N1, N2 )

c  © 2002 – 2014 NorthWest Research Associates, Inc. All Rights Reserved
c    Author: Thomas S. Lund, lund@cora.nwra.com

c   *** Get ranges for x and z variables.

      include 'sam.h'

      logical write_ranges
      character(3)  num

      write_ranges = .false.
c      write_ranges = .true.

      do i=0, numprocs-1
         call range( 1, N1, numprocs, i, iz_s(i), iz_e(i) )
         call range( 1, N2, numprocs, i, ix_s(i), ix_e(i) )
         nx_p( i) = ix_e(i) - ix_s(i) + 1
         nz_p( i) = iz_e(i) - iz_s(i) + 1
      end do

      ixs  = ix_s(myid)
      ixe  = ix_e(myid)
      nxp  = nx_p(myid)
      izs  = iz_s(myid)
      ize  = iz_e(myid)
      nzp  = nz_p(myid)

      if( write_ranges ) then
         write(num,'(i3.3)') myid
         open(unit=80,file='ranges.'//num)
         write(80,5) ixs, ixe, izs, ize, nxp, nzp
5        format(6i6)
         close(80)
      end if

      return
      end

      subroutine range( n1, n2, nprocs, irank, ista, iend )

c  © 2002 – 2014 NorthWest Research Associates, Inc. All Rights Reserved
c    Author: Thomas S. Lund, lund@cora.nwra.com

c   *** The j. tuccillo range getter to balance load.

      iwork1 = (n2 - n1 + 1)/nprocs
      iwork2 = mod(n2 - n1 +1, nprocs)
      ista = irank*iwork1 + n1 + min(irank,iwork2)
      iend = ista + iwork1 - 1
      if(iwork2 .gt. irank) iend = iend + 1

      return
      end
