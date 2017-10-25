      subroutine xy_trans_f( u1, u2, n_var, work )

      include 'sam.h'

      real    u1(Nx,Ny  ,nzp,n_var)
      complex u2(nzp*Ny_min*Nx_min*n_var), work(Ny,Nx2)

c -------- Transform in x and then y and reorder the data to that it is
c -------- ready to go to the z2x_decomp transpose routine.  The output
c -------- array u2 is arranged in the blocks
c -------- u2([nzp*Ny_min*nxp(0)*n_var],[nzp*Ny_min*nxp(1)*n_var]...)

      NzpNy_min = nzp*Ny_min
      NzpNy_minN_var = NzpNy_min*N_var

      js2 = Ny - Ny_min
      je1 = ky_max+1
      jb2 = je1+1

      fn = NxNy_inv

      do n=1, n_var
         do k=1, nzp
            ind1 = k

            do j=1, Ny
               call rfftf( Nx, u1(1,j,k,n), trigx )
               work(j,1) = cmplx( u1(1,j,k,n)*fn, 0.0 )
               do i=2, min( kx_max+1, Nx2 )
                  ir = 2*(i-1)
                  ii = ir + 1
                  work(j,i) = cmplx( u1(ir,j,k,n), u1(ii,j,k,n) )*fn
               end do
               if( kx_max+1 .eq. Nx2+1 ) then
                  work(j,Nx2+1) = cmplx( u1(Nx,j,k,n)*fn, 0.0 )
               end if
            end do

            do m=0, numprocs-1
               ind2 = ind1 + NzpNy_min*nx_p(m)*(n-1) +
     &                NzpNy_minN_var*(ix_s(m)-1)
               do i=ix_s(m), ix_e(m)
                  ind3 = ind2 + (i-ix_s(m))*NzpNy_min
                  call cfftf( Ny, work(1,i), trigy )
                  do j=1, je1
                     ind = ind3 + (j-1)*nzp
                     u2(ind) = work(j,i)
                  end do
                  do j=jb2, Ny_min
                     ind = ind3 + (j-1)*nzp
                     u2(ind) = work(j+js2,i)
                  end do
               end do
            end do

         end do
      end do

      return
      end
