      subroutine z_trans_b( u1, u2, n_var, isymm, u_c, u_r, u_i )

      include 'sam.h'

      complex u1(Nze,Ny_min,nxp,n_var), u2(Ny_min*nxp*Nze*n_var),
     &        u_c(Nze)
      real    u_r(Nze), u_i(Nze)
      integer isymm(n_var)

c -------- Inverse transform in z and reorder the data so that it is
c -------- ready to go to the x2z_decomp transpose routine.  The output
c -------- array u2 is arranged in the blocks
c -------- u2([Ny_min*nxp*nzp(0)*n_var],[Ny_min*nxp*nzp(1)*n_var]...)

      ke1 = kz_max+1
      kb2 = ke1+1
      ks3 = Nze-Nz_min
      ke2 = ke1 + ks3
      kb3 = ke2+1

      Nzm = Nz-1

      Ny_minNxp   = Ny_min*nxp
      Ny_minNxpN_var = Ny_minNxp*n_var

      do n=1, n_var
         do i=1, nxp
            ind1 = Ny_min*(i-1)

            select case( isymm(n) )
            case( 0 )
               do j=1, Ny_min
                  ind2 = ind1 + j
                  do k=1, ke1
                     u_c(k) = u1(k,j,i,n)
                  end do
                  do k=kb2, ke2
                     u_c(k) = cmplx(0.0,0.0)
                  end do
                  do k=kb3, Nz
                     u_c(k) = u1(k-ks3,j,i,n)
                  end do
                  call cfftb( Nz, u_c, trigz )
                  do m=0, numprocs-1
                     ind3 = ind2 + Ny_minNxp*nz_p(m)*(n-1) +
     &                      Ny_minNxpN_var*(iz_s(m)-1)       
                     do k=iz_s(m), iz_e(m)
                        ind = ind3 + Ny_minNxp*(k-iz_s(m))
                        u2(ind) = u_c(k)
                     end do
                  end do
               end do
            case( -1 )
               do j=1, Ny_min
                  ind2 = ind1 + j
                  u_r( 1 ) = 0.0
                  u_i( 1 ) = 0.0
                  do k=2, k1e
                     u_r(k) = real( u1(k,j,i,n))
                     u_i(k) = aimag(u1(k,j,i,n))
                  end do
                  do k=k2b, Nze
                     u_r(k) = 0.0
                     u_i(k) = 0.0
                  end do
                  call sint( Nzm, u_r(2), trigzs )
                  call sint( Nzm, u_i(2), trigzs )
                  do m=0, numprocs-1
                     ind3 = ind2 + Ny_minNxp*nz_p(m)*(n-1) +
     &                      Ny_minNxpN_var*(iz_s(m)-1)
                     do k=iz_s(m), iz_e(m)
                        ind = ind3 + Ny_minNxp*(k-iz_s(m))
                        u2(ind) = cmplx( u_r(k), u_i(k) )
                     end do
                  end do
               end do
            case( 1 )
               do j=1, Ny_min
                  ind2 = ind1 + j
                  do k=1, k1e
                     u_r(k) = real( u1(k,j,i,n))
                     u_i(k) = aimag(u1(k,j,i,n))
                  end do
                  do k=k2b, Nze
                     u_r(k) = 0.0
                     u_i(k) = 0.0
                  end do
                  call cost( Nze, u_r, trigzc )
                  call cost( Nze, u_i, trigzc )
                  do m=0, numprocs-1
                     ind3 = ind2 + Ny_minNxp*nz_p(m)*(n-1) +
     &                      Ny_minNxpN_var*(iz_s(m)-1)
                     do k=iz_s(m), iz_e(m)
                        ind = ind3 + Ny_minNxp*(k-iz_s(m))
                        u2(ind) = cmplx( u_r(k), u_i(k) )
                     end do
                  end do
               end do
            end select

         end do
      end do

      return
      end
