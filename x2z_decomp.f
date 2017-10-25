      subroutine x2z_decomp( u1, u2, n_var, u3 )

c ------ Transpose u1(Ny_min,nxp,Nze,n_var) into u2(Ny,Nx/2,nzp,n_var)
c ------ u3() is workspace that must not overlap u1().  It may be
c ------ the same array as u2().

c ------ On input the u1 array is arranged in blocks like
c ------ [Ny_min,nxp,nzp(0),n_var], [Ny_min,nxp,nzp(1),n_var], ...
c ------ After the all_to_all call, the data is arranged in blocks like
c ------ [Ny_min,nxp(0),nzp,n_var], [Ny_min,nxp(1),nzp,n_var], ...
c ------ We then reorder to get [Ny_min,Nx/2,nzp,n_var].


      include 'sam.h'
      include 'mpif.h'

      complex u1(Ny_min*nxp*Nze*n_var), u2(Ny*Nx2*nzp*n_var),
     &                                  u3(Ny*Nx2*nzp*n_var)
      integer n_send(0:maxp), i_dspl_send(0:maxp),
     &        n_recv(0:maxp), i_dspl_recv(0:maxp)

c ----- Determine send and receive counts as well as their displacements

      n_block_send = Ny_min*nxp*n_var
      n_block_recv = Ny_min*nzp*n_var
      n_send(0) = n_block_send*nz_p(0)
      n_recv(0) = n_block_recv*nx_p(0)
      i_dspl_send(0) = 0
      i_dspl_recv(0) = 0
      do i=1, numprocs-1
         n_send(i) = n_block_send*nz_p(i)
         n_recv(i) = n_block_recv*nx_p(i)
         i_dspl_send(i) = i_dspl_send(i-1) + n_send(i-1)
         i_dspl_recv(i) = i_dspl_recv(i-1) + n_recv(i-1)
      end do

c ----- Swap the data across all processors

      call mpi_alltoallv( u1, n_send, i_dspl_send, mpi_double_complex,
     &                    u3, n_recv, i_dspl_recv, mpi_double_complex,
     &                    mpi_comm_world, ierr )

c ----- Reorder the transpose data, inserting zeros for the truncated
c ----- modes in y.  This is tricky since we want to reorder the data
c ----- in place (u2 and u3 the same array).  Since the output array u2
c ----- is larger than the input u1, we can perform a portion of the
c ----- reordering in place by running the loops backwards over the
c ----- appropriate portion of the data.  The remainder of the data 
c ----- can not be reordered in place and thus it is copied in and out
c ----- of a temporary array.  u1 is used for this purpose, and thus its
c ----- values are altered upon exit.  
c ----- In reality a portion of the data copied to the temporary could
c ----- be reordered in place.  The workable data is located in several
c ----- small chunks, however, which would require complex index 
c ----- bookkeeping.  Reading and writing these small segments is also 
c ----- likely to be inefficient from a cache perspective.  It is 
c ----- probably not worth the effort to try to take advantage of this.
      
      js3 = Ny-Ny_min
      je1 = ky_max+1
      jb2 = je1+1
      je2 = je1+js3
      jb3 = je2+1

c ----- Calculate the fraction of the data that can be reordered in place.
c ----- Here n_stop is the position in the combined index nzp*n_var.

      m = numprocs-1
      n_stop = ( i_dspl_recv(m) + (Nx2-ix_e(m))*Ny )/
     &         ( Ny*Nx2 - Ny_min*nx_p(m) )

      Ny_minNx_min = Ny_min*Nx_min
      NyNx2 = Ny*Nx2

c ----- Load the portion of the data that can not be reordered in place
c ----- into a temporary array.  Use u1() as the temporary.

      do n=1, n_stop
         nm1 = (n-1)/nzp
         km1 = n-1 - nm1*nzp
         ind22 = (n-1)*Ny_minNx_min
         do m=0, numprocs-1
            ind0 = i_dspl_recv(m)
            Ny_minNxp   = Ny_min*nx_p(m)
            ind32 = ind0 + ( nm1*Nzp + km1 )*Ny_minNxp
            do i=ix_s(m), ix_e(m)
               ind23 = ind22 + (i-1      )*Ny_min
               ind33 = ind32 + (i-ix_s(m))*Ny_min
               do j=1, Ny_min
                  ind2 = ind23 + j
                  ind3 = ind33 + j
                  u1(ind2) = u3(ind3)
               end do
            end do
         end do
      end do

c ----- Reorder the remainder of the data in place.

      do n=nzp*n_var, n_stop+1, -1
         nm1 = (n-1)/nzp
         km1 = n-1 - nm1*nzp
         ind22 = (n-1)*NyNx2
         do m=numprocs-1, 0, -1
            ind0 = i_dspl_recv(m)
            Ny_minNxp   = Ny_min*nx_p(m)
            ind32 = ind0 + ( nm1*Nzp + km1 )*Ny_minNxp
            do i=ix_e(m), ix_s(m), -1
               ind23 = ind22 + (i-1      )*Ny
               ind33 = ind32 + (i-ix_s(m))*Ny_min
               do j=Ny, jb3, -1
                  ind2 = ind23 + j
                  ind3 = ind33 + j - js3
                  u2(ind2) = u3(ind3)
               end do
               do j=je2, jb2, -1
                  ind2 = ind23 + j
                  u2(ind2) = cmplx(0.0,0.0)
               end do
               do j=je1, 1, -1
                  ind2 = ind23 + j
                  ind3 = ind33 + j
                  u2(ind2) = u3(ind3)
               end do
            end do
         end do
      end do

c ----- Now copy the data back from the temporary array u1.

      do n=1, n_stop
         ind11 = (n-1)*Ny_minNx_min
         ind21 = (n-1)*NyNx2
         do i=1, Nx_min
            ind12 = ind11 + (i-1)*Ny_min
            ind22 = ind21 + (i-1)*Ny
            do j=1, je1
               ind1 = ind12 + j
               ind2 = ind22 + j
               u2(ind2) = u1(ind1)
            end do
            do j=jb2, je2
               ind2 = ind22 + j
               u2(ind2) = cmplx(0.0,0.0)
            end do
            do j=jb3, Ny
               ind1 = ind12 + j - js3
               ind2 = ind22 + j
               u2(ind2) = u1(ind1)
            end do
         end do
      end do

      return
      end
