      subroutine initial_field( u, random, ierr )

c***********************************************************************
c                                                                      *
c    Initialize a velocity field so that it has prescribed spectrum,   *
c zero divergence and random phase. The velocity field can also be     *
c conditioned by a spherical truncation if desired (see k_truncate in  *
c the argument list description below).                                *
c                                                                      *
c INPUT (ARGUMENT LIST):                                               *
c   i_prob            - flag to declare initial spectrum shape:        *
c                       0-white noise, 1-k^(-5/3), 2-Comte-Bellot      *
c PASSED IN COMMON                                                     *
c   k_truncate        - maximum wavenumber passed by spherical         *
c                       truncation.  Declared as real.  Set > Nx for   *
c                       no truncation.                                 *
c                                                                      *
c OUTPUT                                                               *
c   u(Lx/2,Ly,Lz,3) - velocity field (dimensioned as complex here).  *
c                                                                      *
c***********************************************************************

      include 'sam.h'

      complex u(Nze,Ny_min,nxp,Lu)
      complex  alpha, beta1, div
      real random(Nz_min,Ny_min,3)
      real k_mag, k_mag1
      real E(250,3), S(250)
      integer iseed(2)

      pi = acos(-1.0)
      two_pi = 2.0*pi
      n_div = 0

      time = 0.0
      nt = 0
      t_start  = time
      nt_start = nt

      E = 0.0
      S = 0.0

      factor = 1.0/two_pi
      do i=1, nxp
         ii = ixs + i-1
         fac = 2.0
         if( ii .eq. 1 .or. ii .eq. Nx/2+1 ) fac=1.0
         kx2 = k_x(ii)**2
         rkx = wave_x*float(k_x(ii))
         rkx2 = rkx**2
         iseed(1) = 1234 + ii
         iseed(2) = 6789 + ii
         call random_seed( put=iseed )
         call random_number( random )
         do j=1, Ny_min
            k_sq2 = kx2 + k_y(j)**2
            rky = wave_y*float(k_y(j))
            rk_sq2 = rkx2 + rky**2
            rk12_mag = sqrt(rk_sq2)
            do k=1, Nz_min
               k_sq = k_sq2 + k_z(k)**2
               k_mag = sqrt(float(k_sq))
               k_mag1 = max(k_mag,1.0e-8)
               is = int( k_mag + 0.5 ) + 1
               rkz = wave_z*float(k_z(k))
               rk_sq = rk_sq2 + rkz**2
               rk_mag = sqrt(rk_sq)

c            *** compute energy density from curve fit

               f = sqrt(e_k( k_mag, i_prob )*factor)/k_mag1

c           *** Construct a divergence-free velocity field with random phase

               theta1 = two_pi*random(k,j,1)
               theta2 = two_pi*random(k,j,2)
               phi    = two_pi*random(k,j,3)

               alpha = f*cexp(iunit*theta1)*cos(phi)
               beta1 = f*cexp(iunit*theta2)*sin(phi)

               if( rk12_mag .eq. 0.0 ) then
                  u(k,j,i,1) = alpha
                  u(k,j,i,2) = beta1
                  u(k,j,i,3) = 0.0
               else
                  denom = 1.0/(rk_mag*rk12_mag)
                  u(k,j,i,1) = ( alpha*rk_mag*rky + 
     &                           beta1*   rkx*rkz   )*denom
                  u(k,j,i,2) = ( beta1*   rky*rkz - 
     &                           alpha*rk_mag*rkx   )*denom
                  u(k,j,i,3) = -beta1*rk12_mag / rk_mag
               end if
               div = rkx*u(k,j,i,1) +
     &               rky*u(k,j,i,2) +
     &               rkz*u(k,j,i,3)
               if( cabs(div) .gt. 1.e-8 ) n_div = n_div + 1

               E(is,1) = E(is,1) + f*u(k,j,i,1)*conjg(u(k,j,i,1))
               E(is,2) = E(is,2) + f*u(k,j,i,2)*conjg(u(k,j,i,2))
               E(is,3) = E(is,3) + f*u(k,j,i,3)*conjg(u(k,j,i,3))
               S(is) = S(is) + fac
            end do
         end do
      end do

c      do i=1, 40
c         write(20,20) i-1, (E(i,n),n=1,3), S(i)
c20       format(i5,1p,4e12.4)
c      end do

c   *** Enforce Conjugate Symmetry on plane kx=0

      do i=1, nxp
         ii = ixs + i-1
         if( ii .eq. 1 ) then 
            call enforce_symm( u, i )
         end if
      end do

c   *** Rescale the velocity so that the velocity fluctuation is unity

      if( i_prob .lt. 2 ) call rescale( u )

c   *** Check initial divergence

      call div_check( u, div_max, div_rms, n_div )

      ierr = 0
      if( div_max .gt. 1.e-7 ) then
         if(l_root) print *, 'ERROR in initial_field max divergence = ',
     &              div_max
         ierr = 1
      end if

c   *** Add a mean velocty

c      if( izs .eq. 1 ) then
c         Uo = 1.0e+3
c         u(1,1,1,1) = cmplx(Uo,0.0)
c      end if

      return
      end


      function e_k( k, i_prob )

      real k
      real l, a42(0:6), a98(0:6), a171(0:6)

c   *** Initialize constants for curve fit of Comte-Bellot and Corrsin spectrum

      pi = acos(-1.0)
      l     = 55.0
      u_0   = 27.1893

      const1 = 2.0 * pi / l
      const2 = const1 / ( u_0**2 )

      a42(0) =  0.56102E+01   
      a42(1) = -0.11236E+01
      a42(2) = -0.30961E+00
      a42(3) =  0.33172E+00
      a42(4) = -0.10959E+00
      a42(5) = -0.22320E-01
      a42(6) =  0.66575E-02

      a98(0) =  0.43649E+01
      a98(1) = -0.11793E+01
      a98(2) =  0.67320E-01
      a98(3) =  0.88554E-01
      a98(4) = -0.13372E+00
      a98(5) =  0.28385E-01
      a98(6) = -0.99708E-02

      a171(0) =  0.36567E+01
      a171(1) = -0.11641E+01
      a171(2) = -0.51571E-02
      a171(3) =  0.38064E-02
      a171(4) = -0.10484E+00
      a171(5) =  0.12676E-01
      a171(6) = -0.54562E-02


      if( i_prob .eq. 0 ) then

         e_k = k**2

      else if( i_prob .eq. 1 ) then

         if( k .eq. 0.0 ) then
            e_k = 0.0
         else
            e_k = k**(-5.0/3.0)
         end if

      else if( i_prob .eq. 2 ) then

         if( k .eq. 0.0 ) then
            e_k = 0.0
         else
            sum = 0.0
            do n=0, 6
               sum = sum + a42(n)*alog( const1*k )**n
            end do
            e_k = const2*exp( sum )
         end if

      else if( i_prob .eq. 3 ) then

         if( k .eq. 0.0 ) then
            e_k = 0.0
         else
            sum = 0.0
            do n=0, 6
               sum = sum + a98(n)*alog( const1*k )**n
            end do
            e_k = const2*exp( sum )
         end if

      else if( i_prob .eq. 4 ) then

         if( k .eq. 0.0 ) then
            e_k = 0.0
         else
            sum = 0.0
            do n=0, 6
               sum = sum + a171(n)*alog( const1*k )**n
            end do
            e_k = const2*exp( sum )
         end if

      else

         print *, 'ERROR bad i_prob in e_k'
         stop

      end if

      return
      end


      subroutine enforce_symm( u, i )

      include 'sam.h'

      complex u(Nze,Ny_min,nxp,Lu)

      j=1; k=1
      do n=1, Lu
         u(k,j,i,n) = 0.5*( u(k,j,i,n) + conjg(u(k,j,i,n)) )
      end do

      k=1
      do n=1, Lu
         do j=ky_max+2, Ny_min
            jj = Ny_min+2 - j
            u(k,j,i,n) = conjg(u(k,jj,i,n))
         end do
      end do

      j=1
      do n=1, Lu
         do k=kz_max+2, Nz_min
            kk = Nz_min+2 - k
            u(k,j,i,n) = conjg(u(kk,j,i,n))
         end do
      end do

      do n=1, Lu
         do j=2, Ny_min
            jj = Ny_min+2 - j
            do k=kz_max+2, Nz_min
               kk = Nz_min+2 - k
               u(k,j,i,n) = conjg(u(kk,jj,i,n))
            end do
         end do
      end do

      return
      end 


      subroutine rescale( u )

      include 'sam.h'
      include 'mpif.h'

      complex u(Nze,Ny_min,nxp,Lu)

      sum1 = 0.0
      do i=1, nxp
         do j=1, Ny_min
            do k=1, Nz_min
            sum1 = sum1 + u(k,j,i,1)*conjg(u(k,j,i,1)) +
     &                    u(k,j,i,2)*conjg(u(k,j,i,2)) +
     &                    u(k,j,i,3)*conjg(u(k,j,i,3))
            end do
         end do
      end do

      call mpi_allreduce( sum1, sum1g, 1, mpi_double_precision,
     &                    mpi_sum, mpi_comm_world, ierr )

      u_rms_inv = sqrt(3.0)/sqrt( sum1g  )

      do i=1, nxp
         do j=1, Ny_min
            do k=1, Nz_min
               u(k,j,i,1) = u(k,j,i,1)*u_rms_inv
               u(k,j,i,2) = u(k,j,i,2)*u_rms_inv
               u(k,j,i,3) = u(k,j,i,3)*u_rms_inv
            end do
         end do
      end do

      return
      end
