      subroutine rhs( u, uu, r, work )

      include 'sam.h'

      real u(Nx,Ny,nzp,Lu), uu(Nx,Ny,nzp,Lr), work(Lw)
      complex r(Nze,Ny_min,nxp,Lr+ipad_r), p
      integer i_symm2(10)

      Nr = Lr

c ----- The non-linear terms for the upper half of the matrix are computed
c ----- in column order.  Set up the symmetry factors accordingly.

      L = 0
      do n=1, Lu
         do m=1, n
            L = L + 1
            i_symm2(L) = i_symm(m)*i_symm(n)
         end do
      end do

c   *** Compute velocity maxima, CFL numbers, and timestep

      if( nrk .eq. 1 ) then
         call vel_max( u, work )
      end if

c   *** Compute the non-linear terms

      if( i_strat .eq. 0 ) then
         do k=1, nzp
            do j=1, Ny
               do i=1, Nx
                  uu(i,j,k,6) = u(i,j,k,3)*u(i,j,k,3)
                  uu(i,j,k,5) = u(i,j,k,2)*u(i,j,k,3)
                  uu(i,j,k,4) = u(i,j,k,1)*u(i,j,k,3)
                  uu(i,j,k,3) = u(i,j,k,2)*u(i,j,k,2)
                  uu(i,j,k,2) = u(i,j,k,1)*u(i,j,k,2)
                  uu(i,j,k,1) = u(i,j,k,1)*u(i,j,k,1)
               end do
            end do
         end do
      else
         do k=1, nzp
            do j=1, Ny
               do i=1, Nx
                  uu(i,j,k,9) = u(i,j,k,3)*u(i,j,k,4)
                  uu(i,j,k,8) = u(i,j,k,2)*u(i,j,k,4)
                  uu(i,j,k,7) = u(i,j,k,1)*u(i,j,k,4)
                  uu(i,j,k,6) = u(i,j,k,3)*u(i,j,k,3)
                  uu(i,j,k,5) = u(i,j,k,2)*u(i,j,k,3)
                  uu(i,j,k,4) = u(i,j,k,1)*u(i,j,k,3)
                  uu(i,j,k,3) = u(i,j,k,2)*u(i,j,k,2)
                  uu(i,j,k,2) = u(i,j,k,1)*u(i,j,k,2)
                  uu(i,j,k,1) = u(i,j,k,1)*u(i,j,k,1)
               end do
            end do
         end do
      end if

c   *** Transform the non-linear terms to wave space

      call xy_trans_f( uu, r, Nr, work )

      call z2x_decomp(  r, r, Nr,   uu )

      call z_trans_f( r, Nr, i_symm2, work(1), work(Nze+1) )

c   *** Compute and project the stress divergence.  Save the pressure
c   *** in r(:,:,:,Lu+1).

      do i=1, nxp
         ii = ixs + i-1
         rkx = wave_x*float(k_x(ii))
         rkx2 = rkx**2
         do j=1, Ny_min
            rky = wave_y*float(k_y(j))
            rk_sq2 = rkx2 + rky**2
            do k=1, Nz_min
               rkz = wave_z*float(k_z(k))
               rk_sq = rk_sq2 + rkz**2
               rk_sq = max( rk_sq, 1.0e-10 )
               r(k,j,i,1) = -iunit*( rkx*r(k,j,i,1) +
     &                               rky*r(k,j,i,2) +
     &                               rkz*r(k,j,i,4)   )
               r(k,j,i,2) = -iunit*( rkx*r(k,j,i,2) +
     &                               rky*r(k,j,i,3) +
     &                               rkz*r(k,j,i,5)   )
               r(k,j,i,3) = -iunit*( rkx*r(k,j,i,4) +
     &                               rky*r(k,j,i,5) +
     &                               rkz*r(k,j,i,6)   )
               p = -iunit*( rkx*r(k,j,i,1) + rky*r(k,j,i,2) +
     &                      rkz*r(k,j,i,3) )/rk_sq
               r(k,j,i,1) = r(k,j,i,1) - iunit*rkx*p
               r(k,j,i,2) = r(k,j,i,2) - iunit*rky*p
               r(k,j,i,3) = r(k,j,i,3) - iunit*rkz*p
               r(k,j,i,Lu+1) = p
            end do
         end do
      end do

      if( i_strat .ne. 0 ) then
         do i=1, nxp
            ii = ixs + i-1
            rkx = wave_x*float(k_x(ii))
            do j=1, Ny_min
               rky = wave_y*float(k_y(j))
               do k=1, Nz_min
                  rkz = wave_z*float(k_z(k))
                  r(k,j,i,4) = -iunit*( rkx*r(k,j,i,7) +
     &                                  rky*r(k,j,i,8) +
     &                                  rkz*r(k,j,i,9)   )
               end do
            end do
         end do
      end if

      return
      end
