      subroutine sm( u, u_bar, s_bar, vt )

c***********************************************************************
c                                                                      *
c Conventional Smagorinsky model computation.  In order to save memory,*
c this routine process the field in z=const planes.  In order to       *
c perform the reqired differentiations, 3 adjacent planes of velocity  *
c components are needed.  The work array for velocity is dimensioned   *
c according to (Lx,Ly,3,3), where the first 3 is for the 3 velocity    *
c components and the second 3 is for the 3 planes.  Plane 2 is the     *
c current working plane; planes 1 and 3 are to the left and right      *
c ( k-1, k+1) respectively.  In the main loop, velocity data is added  *
c to plane 3 and the strain rate is computed on plane 2.  The model    *
c is then computed.                                                    *
c                                                                      *
c INPUT (ARGUMENT LIST):                                               *
c   u(Lx,Ly,Lz,3) - Velocity field.                                    *
c   u_bar, s_bar  - Work space arrays dimensioned as shown below       *
c                                                                      *
c OUTPUT:                                                              *
c   vt(Lx,Ly,Lz)  - Eddy viscosity.                                    *
c                                                                      *
c***********************************************************************

      include 'speciso.h'

      real u(Lx,Ly,Lz,3), u_bar(Lx,Ly,3,3), s_bar(Lx,Ly,6),
     &     vt(Lx,Ly,Lz)

c   *** Initialize stuff.

      c_delta_sq = c_smag*delta_sq

      do i=1, Lx*Ly*9
         u_bar(i,1,1,1) = 0.0
      end do

c   *** Process the field in z=const planes.
c   *** Initializing for k=0, 1.  Model computation starts for k=2

      do k=0, Nz+1

         kk = k + 1
         if( kk .lt. 2 ) kk = kk + Nz
         do i=2, Nx+1
            do j=2, Ny+1
               u_bar(i,j,1,3) = 0.5*( u(i,j,kk,1) + u(i-1,j  ,kk  ,1) )
               u_bar(i,j,2,3) = 0.5*( u(i,j,kk,2) + u(i  ,j-1,kk  ,2) )
               u_bar(i,j,3,3) = 0.5*( u(i,j,kk,3) + u(i  ,j  ,kk-1,3) )
            end do
         end do

         call periodicity2( u_bar(1,1,1,3), 3 )

         if( k .ge. 2 ) then

            do i=2, Nx+1
               do j=2, Ny+1
               s_bar(i,j,1) = 
     &           0.5*(u_bar(i+1,j  ,1,2)-u_bar(i-1,j  ,1,2))*dx_inv
               s_bar(i,j,2) = 
     &           0.5*(u_bar(i  ,j+1,2,2)-u_bar(i  ,j-1,2,2))*dy_inv
               s_bar(i,j,3) = 
     &           0.5*(u_bar(i  ,j  ,3,3)-u_bar(i  ,j  ,3,1))*dz_inv
               s_bar(i,j,4) = 
     &         .25*( (u_bar(i  ,j+1,1,2)-u_bar(i  ,j-1,1,2))*dy_inv +
     &               (u_bar(i+1,j  ,2,2)-u_bar(i-1,j  ,2,2))*dx_inv   )
               s_bar(i,j,5) = 
     &         .25*( (u_bar(i  ,j  ,2,3)-u_bar(i  ,j  ,2,1))*dz_inv +
     &               (u_bar(i  ,j+1,3,2)-u_bar(i  ,j-1,3,2))*dy_inv   )
               s_bar(i,j,6) = 
     &         .25*( (u_bar(i+1,j  ,3,2)-u_bar(i-1,j  ,3,2))*dx_inv +
     &               (u_bar(i  ,j  ,1,3)-u_bar(i  ,j  ,1,1))*dz_inv   )
               end do
            end do

            do i=2, Nx+1
               do j=2, Ny+1
                  sum=0.0
                  do n=1, 6
                     sum = sum + ifact(n)*s_bar(i,j,n)**2
                  end do
                  s_mag = sqrt( 2.0*sum )
                  vt(i,j,k) = c_delta_sq*s_mag
               end do
            end do

         end if

         do i=1, Nx+2
            do j=1, Ny+2
               do n=1, 3
                  u_bar(i,j,n,1) = u_bar(i,j,n,2)
                  u_bar(i,j,n,2) = u_bar(i,j,n,3)
               end do
            end do
         end do

      end do

      call periodicity( vt, 1 )

      return
      end
