      subroutine vel_max( u, work )

      include 'sam.h'
      include 'mpif.h'
      parameter( nw=3 )

      real u(Nx,Ny,nzp,Lu), work(Lu,2*nw)

      work(1:Lu,1:nw) = 0.0

      do n=1, Lu
         do k=1, nzp
            do j=1, Ny
               do i=1, Nx
                  work(n,1) = work(n,1) + u(i,j,k,n)
                  work(n,2) = work(n,2) + u(i,j,k,n)**2
                  work(n,3) = max( work(n,3), abs(u(i,j,k,n)) )
               end do
            end do
         end do
      end do

      call mpi_allreduce( work(1,1), work(1,nw+1), Lu*2,
     &     mpi_double_precision, mpi_sum, mpi_comm_world, ierr )
      call mpi_allreduce( work(1,3), work(1,nw+3), Lu,
     &     mpi_double_precision, mpi_max, mpi_comm_world, ierr )

      do n=1, Lu
         u_avg(n) = work(n,nw+1)*NxNyNz_inv
         u_var(n) = work(n,nw+2)*NxNyNz_inv - u_avg(n)**2
         u_max(n) = work(n,nw+3)
      end do

      energy = 0.5*( u_var(1) + u_var(2) + u_var(3) )

      cfl_x_by_dt = u_max(1)*kx_max*wave_x
      cfl_y_by_dt = u_max(2)*ky_max*wave_y
      cfl_z_by_dt = u_max(3)*kz_max*wave_z

      if( i_cfl .eq. 1 ) then
c         cfl_by_dt = sqrt( cfl_x_by_dt**2 + cfl_y_by_dt**2 +
c     &                     cfl_z_by_dt**2 )
         cfl_by_dt = max( cfl_x_by_dt, cfl_y_by_dt, cfl_z_by_dt )
         dt = cfl0*adv_cfl_limit/cfl_by_dt
      end if

      cfl_x = cfl_x_by_dt*dt
      cfl_y = cfl_y_by_dt*dt
      cfl_z = cfl_z_by_dt*dt

      return
      end
