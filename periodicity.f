      subroutine periodicity( f, n )

      include 'speciso.h'

      real f(Lx,Ly,Lz,n)

      do j=2, Ny+1
         do k=2, Nz+1
            do m=1, n
               f(Nx+2,j,k,m) = f(2   ,j,k,m)
               f(1   ,j,k,m) = f(Nx+1,j,k,m)
            end do
         end do
      end do

      do i=1, Nx+2
         do k=2, Nz+1
            do m=1, n
               f(i,Ny+2,k,m) = f(i,2   ,k,m)
               f(i,1   ,k,m) = f(i,Ny+1,k,m)
            end do
         end do
      end do

      do i=1, Nx+2
         do j=1, Ny+2
            do m=1, n
               f(i,j,Nz+2,m) = f(i,j,2   ,m)
               f(i,j,1   ,m) = f(i,j,Nz+1,m)
            end do
         end do
      end do

      return
      end
