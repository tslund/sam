      subroutine periodicity2( f, n_var )

      include 'speciso.h'

      real f(Lx,Ly,n_var)

      do j=2, Ny+1
         do n=1, n_var
            f(Nx+2,j,n) = f(2   ,j,n)
            f(1   ,j,n) = f(Nx+1,j,n)
         end do
      end do

      do i=1, Nx+2
         do n=1, n_var
            f(i,Ny+2,n) = f(i,2   ,n)
            f(i,1   ,n) = f(i,Ny+1,n)
         end do
      end do

      return
      end
