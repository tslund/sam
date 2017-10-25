      subroutine z2x_decomp( u1, u2, n_var, u3 )

c ------ Transpose u1(nzp,Ny_min,Nx_min,n_var) into u2(Nze,Ny_min,nxp,n_var)
c ------ u3() is workspace that must not overlap u1() or u2().
c ------ u1() and u2() can be the same array.

c ------ On input the u1 array is arranged in blocks like
c ------ [nzp,Ny_min,nxp(0),n_var], [nzp,Ny_min,nxp(1),n_var], ...
c ------ After the all_to_all call, the data is arranged in blocks like
c ------ [nzp(0),Ny_min,nxp,n_var], [nzp(1),Ny_min,nxp,n_var], ...
c ------ We then reorder to get [Nze,Ny_min,nxp,n_var].

      include 'sam.h'
      include 'mpif.h'

      complex u1(nzp*Ny_min*Nx_min*n_var), u2(Nze*Ny_min*nxp*n_var),
     &        u3(Nze*Ny_min*nxp*n_var)
      integer n_send(0:numprocs-1), i_dspl_send(0:numprocs-1),
     &        n_recv(0:numprocs-1), i_dspl_recv(0:numprocs-1)

c ----- Determine send and receive counts as well as their displacements

      n_block_send = Ny_min*nzp*n_var
      n_block_recv = Ny_min*nxp*n_var
      n_send(0) = n_block_send*nx_p(0)
      n_recv(0) = n_block_recv*nz_p(0)
      i_dspl_send(0) = 0
      i_dspl_recv(0) = 0
      do i=1, numprocs-1
         n_send(i) = n_block_send*nx_p(i)
         n_recv(i) = n_block_recv*nz_p(i)
         i_dspl_send(i) = i_dspl_send(i-1) + n_send(i-1)
         i_dspl_recv(i) = i_dspl_recv(i-1) + n_recv(i-1)
      end do

c ----- Swap data across all processors

      call mpi_alltoallv( u1, n_send, i_dspl_send, mpi_double_complex,
     &                    u3, n_recv, i_dspl_recv, mpi_double_complex,
     &                    mpi_comm_world, ierr )

c ----- Reorder the transpose data

      NzeNy_min = Nze*Ny_min
      NzeNy_minNxp = NzeNy_min*Nxp

      do m=0, numprocs-1
         NzpNy_min = nz_p(m)*Ny_min
         NzpNy_minNxp = NzpNy_min*nxp
         ind0 = i_dspl_recv(m)
         do n=1, n_var
            ind21 =        NzeNy_minNxp*(n-1)
            ind31 = ind0 + NzpNy_minNxp*(n-1)
            do i=1, nxp
               ind22 = ind21 + (i-1)*NzeNy_min
               ind32 = ind31 + (i-1)*NzpNy_min
               do j=1, Ny_min
                  ind23 = ind22 + (j-1)*Nze
                  ind33 = ind32 + (j-1)*nz_p(m)
                  do k=iz_s(m), iz_e(m)
                     ind2 = ind23 + k
                     ind3 = ind33 + k-iz_s(m)+1
                     u2(ind2) = u3(ind3)
                  end do
               end do
            end do
         end do
      end do

      return
      end
