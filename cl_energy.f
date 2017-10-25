      subroutine cl_energy( u )

      include 'sam.h'
      include 'mpif.h'

      complex u(Nze,Ny_min,nxp,Lu)

      cl_start_sq = (float(k_cl)-0.5)**2
      cl_stop_sq  = (float(k_cl)+0.5)**2

      sum_e = 0.0
      do i=1, nxp
         ii = ixs + i-1
         factor = 2.0
         if( ii .eq. 1 .or. ii .eq. Nx/2+1 ) factor = 1.0
         kx2 = k_x(ii)**2
         do j=1, Ny_min
            k_sq2 = kx2 + k_y(j)**2
            if( k_sq2 .gt. cl_start_sq ) then
               k_start = 1
            else
               k_start = int(sqrt(cl_start_sq-k_sq2)) + 2
            end if
            if( k_sq2 .gt. cl_stop_sq  ) then
               k_stop  = 0
            else
               k_stop  = int(sqrt(cl_stop_sq -k_sq2)) + 1
            end if
            do k=k_start, k_stop
               kk = Nz_min+2 - k
               u_sq1 = u( k,j,i,1)*conjg(u( k,j,i,1)) +
     &                 u( k,j,i,2)*conjg(u( k,j,i,2)) +
     &                 u( k,j,i,3)*conjg(u( k,j,i,3))
               u_sq2 = u(kk,j,i,1)*conjg(u(kk,j,i,1)) +
     &                 u(kk,j,i,2)*conjg(u(kk,j,i,2)) +
     &                 u(kk,j,i,3)*conjg(u(kk,j,i,3))
               sum_e = sum_e + factor*( u_sq1 + u_sq2 )
            end do
         end do
      end do

      call mpi_allreduce( sum_e, E_cl, 1, mpi_double_precision,
     &                    mpi_sum, mpi_comm_world, ierr )

      return
      end
