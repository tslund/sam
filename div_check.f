      subroutine div_check( u, div_max, div_rms, n_div )

      include 'sam.h'
      include 'mpif.h'

      complex u(Nze,Ny_min,nxp,Lu), div
      real sum1(4), sum1_g(4)

      sum1 = 0.0

      do i=1, nxp
         ii = ixs + i-1
         rkx = wave_x*float(k_x(ii))
         do j=1, Ny_min
            rky = wave_y*float(k_y(j))
            do k=1, Nz_min
               rkz = wave_z*float(k_z(k))
               div = rkx*u(k,j,i,1) +
     &               rky*u(k,j,i,2) +
     &               rkz*u(k,j,i,3)
               div_sq = div*conjg(div)
               sum1(1) = sum1(1) + sqrt(div_sq)
               sum1(2) = sum1(2) +      div_sq
               sum1(3) = max( sum1(3), div_sq )
               if( div_sq .gt. 1.e-20 ) sum1(4) = sum1(4) + 1
            end do
         end do
      end do

      call mpi_allreduce( sum1, sum1_g, 4, mpi_double_precision,
     &                    mpi_sum, mpi_comm_world, ierr )

      div_avg = sum1_g(1)
      div_rms = sqrt(sum1_g(2))
      div_max = sqrt(sum1_g(3))
      n_div   = sum1_g(4)

      return
      end
