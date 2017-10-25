      subroutine dm( u, u_bar, s_bar, s_mag_s, u_bar_u_bar,
     &               u_hat, s_hat, L, M, vt )

c***********************************************************************
c                                                                      *
c Dynamic Smagorinsky model computation.  In order to save memory,     *
c this routine process the field in z=const planes.  In order to       *
c perform the reqired differentiations and filterings, 4 adjacent      *
c planes of velocity components and 3 adjacent planes of most other    *
c quantities are needed.  The work arrays for multi-plane data are     *
c dimensioned according to (Lx,Ly,n_var,n_plane), where n_var is the   *
c number of variables (3 for velocity, 6 for symmetric tensor) and     *
c n_plane is the number of adjacent planes (usually 3 or 4).  Plane 2  *
c is the current working plane; planes 1 and 3 are to the left and     *
c right ( k-1, k+1) respectively.  In the main loop, velocity data is  *
c added to plane 4 and the strain rate is computed on plane 3.  All    *
c quantities were computed on planes 1 and 2 on previous passes through*
c the loop and thus it is possible to compute the tensors L and M at   *
c plane 2.  This strategy requires several initialization passes,      *
c but those occur automatically when k<2 in the main loop.             *
c   Note that the velocity components are averaged to the cell center  *
c at the start of the computation.  This step makes the procedure much *
c simpler and has virtually no effect on the volume-averaged value of  *
c the Smagorinsky constant.                                            *
c                                                                      *
c INPUT (ARGUMENT LIST):                                               *
c   u(Lx,Ly,Lz,3) - Velocity field.                                    *
c   u_bar-M       - Work space arrays dimensioned as shown below       *
c                                                                      *
c OUTPUT:                                                              *
c   vt(Lx,Ly,Lz)  - Eddy viscosity.                                    *
c                                                                      *
c***********************************************************************

      include 'speciso.h'

      real u(Lx,Ly,Lz,3), u_bar(Lx,Ly,3,4), s_bar(Lx,Ly,6,3),
     &     s_mag_s(Lx,Ly,6,3), u_bar_u_bar(Lx,Ly,6,3),
     &     u_hat(Lx,Ly,3), s_hat(Lx,Ly,6), 
     &     L(Lx,Ly,6), M(Lx,Ly,6), vt(Lx,Ly,Lz)
      real num

c   *** Initialize stuff.

      num = 0.0
      den = 0.0

      do i=1, Lx*Ly*12
         u_bar(i,1,1,1) = 0.0
      end do
      do i=1, Lx*Ly*18
         s_bar(      i,1,1,1) = 0.0
         s_mag_s(    i,1,1,1) = 0.0
         u_bar_u_bar(i,1,1,1) = 0.0
      end do

c   *** Process the field in z=const planes.
c   *** Initializing for k=-3,...,1.  Model computation starts for k=2

      do k=-3, Nz+1

         kk = k + 2
         if( kk .lt.    2 ) kk = kk + Nz
         if( kk .gt. Nz+2 ) kk = kk - Nz
         do i=2, Nx+1
            do j=2, Ny+1
               u_bar(i,j,1,4) = 0.5*( u(i,j,kk,1) + u(i-1,j  ,kk  ,1) )
               u_bar(i,j,2,4) = 0.5*( u(i,j,kk,2) + u(i  ,j-1,kk  ,2) )
               u_bar(i,j,3,4) = 0.5*( u(i,j,kk,3) + u(i  ,j  ,kk-1,3) )
            end do
         end do

         call periodicity2( u_bar(1,1,1,4), 3 )

         if( k .ge. 0 ) then

            do i=2, Nx+1
               do j=2, Ny+1
               s_bar(i,j,1,3) = 
     &           0.5*(u_bar(i+1,j  ,1,3)-u_bar(i-1,j  ,1,3))*dx_inv
               s_bar(i,j,2,3) = 
     &           0.5*(u_bar(i  ,j+1,2,3)-u_bar(i  ,j-1,2,3))*dy_inv
               s_bar(i,j,3,3) = 
     &           0.5*(u_bar(i  ,j  ,3,4)-u_bar(i  ,j  ,3,2))*dz_inv
               s_bar(i,j,4,3) = 
     &         .25*( (u_bar(i  ,j+1,1,3)-u_bar(i  ,j-1,1,3))*dy_inv +
     &               (u_bar(i+1,j  ,2,3)-u_bar(i-1,j  ,2,3))*dx_inv   )
               s_bar(i,j,5,3) = 
     &         .25*( (u_bar(i  ,j  ,2,4)-u_bar(i  ,j  ,2,2))*dz_inv +
     &               (u_bar(i  ,j+1,3,3)-u_bar(i  ,j-1,3,3))*dy_inv   )
               s_bar(i,j,6,3) = 
     &         .25*( (u_bar(i+1,j  ,3,3)-u_bar(i-1,j  ,3,3))*dx_inv +
     &               (u_bar(i  ,j  ,1,4)-u_bar(i  ,j  ,1,2))*dz_inv   )
               end do
            end do

            call periodicity2( s_bar(1,1,1,3), 6 )

            do i=2, Nx+1
               do j=2, Ny+1
                  sum=0.0
                  do n=1, 6
                     sum = sum + ifact(n)*s_bar(i,j,n,3)**2
                  end do
                  s_mag = sqrt( 2.0*sum )
                  if( k .ge. 1 ) vt(i,j,k+1) = s_mag
                  do n=1, 6
                     s_mag_s(i,j,n,3) = s_mag*s_bar(i,j,n,3)
                     u_bar_u_bar(i,j,n,3) = u_bar(i,j,ip(n),3) *
     &                                      u_bar(i,j,jp(n),3)
                  end do
               end do
            end do

            call periodicity2( s_mag_s(    1,1,1,3), 6 )
            call periodicity2( u_bar_u_bar(1,1,1,3), 6 )

         end if

         if( k .ge. 2 ) then

            call filter( u_bar(1,1,1,1), u_hat, 3 )
            call filter( s_bar(1,1,1,1), s_hat, 6 )

            call filter( s_mag_s(    1,1,1,1), M, 6 )
            call filter( u_bar_u_bar(1,1,1,1), L, 6 )

            do i=2, Nx+1
               do j=2, Ny+1
                  sum=0.0
                  do n=1, 6
                     sum = sum + ifact(n)*s_hat(i,j,n)**2
                  end do
                  s_mag = sqrt( 2.0*sum )
                  do n=1, 6
                     M(i,j,n) = -M(i,j,n) + alpha_sq*s_mag*s_hat(i,j,n)
                     L(i,j,n) =  L(i,j,n) - u_hat(i,j,ip(n)) *
     &                                      u_hat(i,j,jp(n))
                  end do
               end do
            end do

            do i=2, Nx+1
               do j=2, Ny+1
                  do n=1, 6
                     num = num + ifact(n)*L(i,j,n)*M(i,j,n)
                     den = den + ifact(n)*M(i,j,n)*M(i,j,n)
                  end do
               end do
            end do

         end if

         do i=1, Nx+2
            do j=1, Ny+2
               do n=1, 3
                  u_bar(i,j,n,1) = u_bar(i,j,n,2)
                  u_bar(i,j,n,2) = u_bar(i,j,n,3)
                  u_bar(i,j,n,3) = u_bar(i,j,n,4)
               end do
               do n=1, 6
                  s_bar(      i,j,n,1) = s_bar(      i,j,n,2)
                  s_bar(      i,j,n,2) = s_bar(      i,j,n,3)
                  s_mag_s(    i,j,n,1) = s_mag_s(    i,j,n,2)
                  s_mag_s(    i,j,n,2) = s_mag_s(    i,j,n,3)
                  u_bar_u_bar(i,j,n,1) = u_bar_u_bar(i,j,n,2)
                  u_bar_u_bar(i,j,n,2) = u_bar_u_bar(i,j,n,3)
               end do
            end do
         end do

      end do

c   *** Compute volume-averaged model coefficient

      c_delta_sq = -0.5*num / den

      write(8,10) nt, t, c_delta_sq*delta_sq_inv
10    format(i5, 1p2e12.4)

c   *** Compute eddy viscosity (strain mag was loaded previously)

      do i=2, Nx+1
         do j=2, Ny+1
            do k=2, Nz+1
               vt(i,j,k) = c_delta_sq*vt(i,j,k)
            end do
         end do
      end do

      call periodicity( vt, 1 )

      return
      end
