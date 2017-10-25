      subroutine force( u, ek_0, work, ek_g )

      include 'sam.h'
      include 'mpif.h'

      complex u(Nze,Ny_min,nxp,Lu)
      real    ek_0(Lm), work(Lm), ek_g(Lm)
      real    k_mag

c   *** Compute energy in spherical shells up to the forcing radius

      work = 0.0
      do i=1, nxp
         ii = ixs + i-1
         kx2 = k_x(ii)**2
         do j=1, Ny_min
            k_sq2 = kx2 + k_y(j)**2
            do k=1, Nz_min
               factor = 2.0
               if( ii .eq. 1 .or. ii .eq. Nx2+1 ) factor = 1.0
               k_sq  = k_sq2 + k_z(k)**2
               k_mag = sqrt(float(k_sq))
               is = 1 + int( k_mag + 0.5 )
               if( is .le. k_force+1 ) then
                  u_sq = u(k,j,i,1)*conjg(u(k,j,i,1)) +
     &                   u(k,j,i,2)*conjg(u(k,j,i,2)) +
     &                   u(k,j,i,3)*conjg(u(k,j,i,3))
                  work(is) = work(is) + factor*u_sq
               end if
            end do
         end do
      end do

      call mpi_allreduce( work, ek_g, k_force+1, mpi_double_precision,
     &                    mpi_sum, mpi_comm_world, ierr )

c   *** Multiply the fourier coefficients by a correction factor to
c   *** keep the spectrum fixed

      work(1) = 1.0
      do is=2, k_force+1
         work(is) = sqrt( ek_0(is)/ek_g(is) )
      end do

      do i=1, nxp
         ii = ixs + i-1
         kx2 = k_x(ii)**2
         do j=1, Ny_min
            k_sq2 = kx2 + k_y(j)**2
            do k=1, Nz_min
               factor = 2.0
               if( ii .eq. 1 .or. ii .eq. Nx2+1 ) factor = 1.0
               k_sq  = k_sq2 + k_z(k)**2
               k_mag = sqrt(float(k_sq))
               is = 1 + int( k_mag + 0.5 )
               if( is .le. k_force+1 ) then
                  u(k,j,i,1) = work(is)*u(k,j,i,1)
                  u(k,j,i,2) = work(is)*u(k,j,i,2)
                  u(k,j,i,3) = work(is)*u(k,j,i,3)
               end if
            end do
         end do 
      end do 

      return
      end
