      subroutine filter( f_bar, f_hat, n_var )

      include 'speciso.h'

      real f_bar(Lx,Ly,n_var,3), f_hat(Lx,Ly,n_var)

      w_inv = 1.0 / ( 1.0 + 4.0 + 1.0 )**3

      do n=1, n_var
         do i=1, Nx+2
            do j=1, Ny+2
               f_hat(i,j,n) = f_bar(i,j,n,1) + 4*f_bar(i,j,n,2) +
     &                        f_bar(i,j,n,3)
            end do
         end do
      end do

      do n=1, n_var
         do i=2, Nx+1
            do j=1, Ny+2
               w_r(j) = f_hat(i,j,n)
            end do
            do j=2, Ny+1
               f_hat(i,j,n) = w_r(j-1) + 4*w_r(j) + w_r(j+1)
            end do
         end do
      end do

      call periodicity2( f_hat, n_var )

      do n=1, n_var
         do j=2, Ny+1
            do i=1, Nx+2
               w_r(i) = f_hat(i,j,n)
            end do
            do i=2, Nx+1
               f_hat(i,j,n) = w_inv*( w_r(i-1) + 4*w_r(i) + w_r(i+1) )
            end do
         end do
      end do

      call periodicity2( f_hat, n_var )

      return
      end
