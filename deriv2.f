      subroutine deriv2( u, du, Nx, Ny, trigy, k_y, ky_max, wave_y,
     &                   work )

      real u(Nx,Ny), du(Nx,Ny), trigy(*), work(Ny), Ny_inv
      integer k_y(Ny)

      Ny_inv = 1.0/float(Ny)

      do i=1, Nx
         do j=1, Ny
            work(j) = u(i,j)*Ny_inv
         end do
         call rfftf( Ny, work, trigy )
         work(1) = 0.0
         do j=2, ky_max+1
            jr = 2*(j-1)
            ji = jr+1
            rky = k_y(j)*wave_y
            tmp = work(jr)
            work(jr) = -rky*work(ji)
            work(ji) =  rky*tmp
         end do
         do j=2*ky_max+2, Ny
            work(j) = 0.0
         end do
         call rfftb( Ny, work, trigy )
         do j=1, Ny
            du(i,j) = work(j)
         end do
      end do

      return
      end
