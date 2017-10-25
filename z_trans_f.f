      subroutine z_trans_f( u, n_var, isymm, u_r, u_i )

      include 'sam.h'

      complex u(Nze,Ny_min,nxp,n_var)
      real    u_r(Nze), u_i(Nze)
      integer isymm(n_var)

c -------- Forward transform in z, truncate and store coefficients on
c -------- the minimal grid in z.

      ke1 = kz_max+1
      kb2 = ke1+1
      ks2 = Nze-Nz_min

      Nzm = Nz-1

      do n=1, n_var
         do i=1, nxp

            select case( isymm(n) )
            case( 0 )
               do j=1, Ny_min
                  call cfftf( Nz, u(1,j,i,n), trigz )
                  do k=1, ke1
                     u(k,j,i,n) = u(k    ,j,i,n)*Nz_inv
                  end do
                  do k=kb2, Nz_min
                     u(k,j,i,n) = u(k+ks2,j,i,n)*Nz_inv
                  end do
               end do
            case( -1 )
               do j=1, Ny_min
                  u_r( 1 ) = 0.0
                  u_i( 1 ) = 0.0
                  do k=2, Nz
                     u_r(k) = real( u(k,j,i,n))
                     u_i(k) = aimag(u(k,j,i,n))
                  end do
                  u_r(Nze) = 0.0
                  u_i(Nze) = 0.0
                  call sint( Nzm, u_r(2), trigzs )
                  call sint( Nzm, u_i(2), trigzs )
                  do k=1, ke1
                     u(k,j,i,n) = cmplx( u_r(k), u_i(k) )*Nz_inv
                  end do
                  do k=kb2, Nze
                     u(k,j,i,n) = cmplx( 0.0, 0.0 )
                  end do
               end do
            case( 1 )
               do j=1, Ny_min
                  do k=1, Nze
                     u_r(k) = real( u(k,j,i,n))
                     u_i(k) = aimag(u(k,j,i,n))
                  end do
                  call cost( Nze, u_r, trigzc )
                  call cost( Nze, u_i, trigzc )
                  do k=1, ke1
                     u(k,j,i,n) = cmplx( u_r(k), u_i(k) )*Nz_inv
                  end do
                  do k=kb2, Nze
                     u(k,j,i,n) = cmplx( 0.0, 0.0 )
                  end do
               end do
            end select

         end do
      end do

      return
      end
