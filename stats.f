      subroutine stats( u, uiuj, strain, vt, k )

      include 'speciso.h'

      real u(Lx,Ly,Lz,3), uiuj(Lx,Ly,6), strain(Lx,Ly,6),
     &     vt(Lx,Ly)
c      real p(Lx,Ly,Lz)

      do n=1, Ls
         stat_tmp(n) = 0.0
      end do

      do i=1, Nx
         do j=1, Ny

            stat_tmp(1 ) = stat_tmp(1 ) + u(i,j,k,1)
            stat_tmp(2 ) = stat_tmp(2 ) + u(i,j,k,2)
            stat_tmp(3 ) = stat_tmp(3 ) + u(i,j,k,3)

c            stat_tmp(4 ) = stat_tmp(4 ) + p(i,j,k  )

            stat_tmp(5 ) = stat_tmp(5 ) + uiuj(i,j,1)
            stat_tmp(6 ) = stat_tmp(6 ) + uiuj(i,j,2)
            stat_tmp(7 ) = stat_tmp(7 ) + uiuj(i,j,3)
            stat_tmp(8 ) = stat_tmp(8 ) + uiuj(i,j,4)
            stat_tmp(9 ) = stat_tmp(9 ) + uiuj(i,j,5)
            stat_tmp(10) = stat_tmp(10) + uiuj(i,j,6)

c            stat_tmp(11) = stat_tmp(11) + p(i,j,k)**2

c            stat_tmp(12) = stat_tmp(12) + u(i,j,k,1)*p(i,j,k)
c            stat_tmp(13) = stat_tmp(13) + u(i,j,k,2)*p(i,j,k)
c            stat_tmp(14) = stat_tmp(14) + u(i,j,k,3)*p(i,j,k)

            stat_tmp(15) = stat_tmp(15) + strain(i,j,1)
            stat_tmp(16) = stat_tmp(16) + strain(i,j,2)
            stat_tmp(17) = stat_tmp(17) + strain(i,j,3)
            stat_tmp(18) = stat_tmp(18) + strain(i,j,4)
            stat_tmp(19) = stat_tmp(19) + strain(i,j,5)
            stat_tmp(20) = stat_tmp(20) + strain(i,j,6)

            stat_tmp(21) = stat_tmp(21) + vt(i,j)

            stat_tmp(22) = stat_tmp(22) + 2*vt(i,j)*strain(i,j,1)
            stat_tmp(23) = stat_tmp(23) + 2*vt(i,j)*strain(i,j,2)
            stat_tmp(24) = stat_tmp(24) + 2*vt(i,j)*strain(i,j,3)
            stat_tmp(25) = stat_tmp(25) + 2*vt(i,j)*strain(i,j,4)
            stat_tmp(26) = stat_tmp(26) + 2*vt(i,j)*strain(i,j,5)
            stat_tmp(27) = stat_tmp(27) + 2*vt(i,j)*strain(i,j,6)

            stat_tmp(28) = stat_tmp(28) + 2*vis*strain(i,j,1)
            stat_tmp(29) = stat_tmp(29) + 2*vis*strain(i,j,2)
            stat_tmp(30) = stat_tmp(30) + 2*vis*strain(i,j,3)
            stat_tmp(31) = stat_tmp(31) + 2*vis*strain(i,j,4)
            stat_tmp(32) = stat_tmp(32) + 2*vis*strain(i,j,5)
            stat_tmp(33) = stat_tmp(33) + 2*vis*strain(i,j,6)

            sum = 0.0
            do n=1, 6
               sum = sum + ifact(n)*strain(i,j,n)**2
            end do
            stat_tmp(34) = stat_tmp(34) + sum
            stat_tmp(35) = stat_tmp(35) + 2*vt(i,j)*sum

         end do
      end do

c      *** Add statistics for the current time step to the running time average

      if( k .eq. Nz ) then
         wt = dt*NxNyNz_inv
         do n=1, Ls
            stat(n) = stat(n) + wt*stat_tmp(n)
         end do
      end if

      return
      end
