      subroutine deriv1( u, du, Nx, Ny, trigx, k_x, kx_max, wave_x,
     &                   work )

      real u(Nx,Ny), du(Nx,Ny), trigx(*), work(*), Nx_inv
      integer k_x(Nx)

      Nx_inv = 1.0/float(Nx)

      do j=1, Ny
         do i=1, Nx
            du(i,j) = u(i,j)*Nx_inv
         end do
         call rfftf( Nx, du(1,j), trigx )
         du(1,j) = 0.0
         do i=2, kx_max+1
            ir = 2*(i-1)
            ii = ir+1
            rkx = k_x(i)*wave_x
            tmp = du(ir,j)
            du(ir,j) = -rkx*du(ii,j)
            du(ii,j) =  rkx*tmp
         end do
         do i=2*kx_max+2, Nx
            du(i,j) = 0.0
         end do
         call rfftb( Nx, du(1,j), trigx )
      end do

      return
      end
