      subroutine initial_field( u, random, ierr )

c Initialize a velocity field so that it has prescribed spectrum,
c zero divergence and random phase. The velocity field can also be
c conditioned by a spherical truncation if desired (see k_truncate in
c the argument list description below).
c
c INPUT (ARGUMENT LIST):
c   i_prob            - flag to control the flow to be initialized:
c                       0-white noise
c                       1-k^(-5/3)
c                       2-Comte-Bellot
c                       3-Pao spectrum
c                       4-gravity wave in a tilted box
c                      -4-gravity wave in a non-tilted box
c                       5-shear profile in the y-z plane
c OUTPUT
c   u(Nze,Ny_min,nxp,Lu) - velocity field

      include 'sam.h'

      complex u(Nze,Ny_min,nxp,Lu)
      complex  alpha, beta1, div
      real random(Nze,Ny_min,3)
      real k_mag, k_mag1, lambda
      real E(250,3), S(250)
      integer, allocatable, dimension (:) ::  iseed

      pi = acos(-1.0)
      two_pi = 2.0*pi
      n_div = 0

      time = 0.0
      nt = 0
      t_start  = time
      nt_start = nt
      lapse = lapse0
      buoy_fac_x = 0.0
      buoy_fac_z = 0.0
      n_frame_p = 0

      select case( i_prob )
         case(0,4,-4,5); i_type = 0	! white noise
         case(1  )     ; i_type = 1	! k^-5/3
         case(2  )     ; i_type = 2	! cbc t=42
         case(3  )     ; i_type = 5	! Pao
         case default  ; i_type = 0
      end select

      call random_seed( size=L_seed)
      allocate( iseed(L_seed) )

      if( flct_u .eq. 0.0 .and. i_prob .ne. 2 .and. i_prob .ne. 3 ) then
         u(1:Nze,1:Ny_min,1:nxp,1:3) = cmplx(0.0,0.0)
         goto 50
      end if

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
         call set_seed( iseed, L_seed, ii )
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

               f = sqrt(e_k( k_mag, i_type )*factor)/k_mag1

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
            call enforce_symm( u, 3, i )
         end if
      end do

c   *** Rescale the velocity so that the fluctuation is flct_u

      if( i_prob .ne. 2 ) call rescale( u, 3, flct_u )

c   *** Check initial divergence

      call div_check( u, div_max, div_rms, n_div )

      ierr = 0
      if( div_max .gt. 1.e-7 ) then
         if(l_root) print *, 'ERROR in initial_field max divergence = ',
     &              div_max
         ierr = 1
      end if

50    continue

c   *** Initialize the temperature fluctuations if required.

      if( i_strat .eq. 0 ) goto 60

      buoy_fac = To_inv*grav
      buoy_fac_x = 0.0
      buoy_fac_z = buoy_fac

      if( flct_T .eq. 0.0 ) then
         u(1:Nze,1:Ny_min,1:nxp,4) = cmplx(0.0,0.0)
         goto 60
      end if

      do i=1, nxp
         ii = ixs + i-1
         kx2 = k_x(ii)**2
         call set_seed( iseed, L_seed, -ii )
         call random_seed( put=iseed )
         call random_number( random )
         do j=1, Ny_min
            k_sq2 = kx2 + k_y(j)**2
            do k=1, Nz_min
               k_sq = k_sq2 + k_z(k)**2
               k_mag = sqrt(float(k_sq))
               k_mag1 = max(k_mag,1.0e-8)
               f = sqrt(e_k( k_mag, 0 )*factor)/k_mag1
               theta = two_pi*random(k,j,1)
               u(k,j,i,4) = f*cexp(iunit*theta)
            end do
         end do
      end do

c   *** Enforce Conjugate Symmetry on plane kx=0

      do i=1, nxp
         ii = ixs + i-1
         if( ii .eq. 1 ) then
            call enforce_symm( u(1,1,1,4), 1, i )
         end if
      end do

c   *** Rescale the temperature so that the fluctuation is flct_t

      if( i_prob .ne. 2 ) call rescale( u(1,1,1,4), 1, flct_t )

60    continue

c   *** Add a mean velocty

      if( ixs .eq. 1 ) then
         u(1,1,1,1) = cmplx(Uo,0.0)
      end if

c   *** Add a fine-scale structure

      if( i_fs .ne. 0 ) then

         iz_fs1 = nint(zL/lambda_z_fs) + 1
         iz_fs2 = Nz_min+2 - iz_fs1
         amp_fs = lambda_z_fs/(2.0*pi)*sqrt(N_sq/Ri_fs)

         if( ixs .eq. 1 ) then
            u(iz_fs1,1,1,1) = u(iz_fs1,1,1,1) + amp_fs*cmplx(0.0, 0.5)
            u(iz_fs2,1,1,1) = u(iz_fs2,1,1,1) + amp_fs*cmplx(0.0,-0.5)
         end if

      end if

c   *** Initialize a gravity wave

      if( i_prob .eq. -4 ) then

         ix_w = nint(xL/lambda_x) + 1
         iz_w = nint(zL/lambda_z) + 1

         c_i = omega_i/k_w
         amp_u = amplitude*c_i
         amp_T = amplitude*lapse0/m_w
         r1 = k_w/m_w

         do i=1, nxp
            ii = ixs + i-1
            if( ii .eq. ix_w ) then
               u(iz_w,1,i,1) = u(iz_w,1,i,1) +    amp_u*cmplx( 0.0,-0.5)
               u(iz_w,1,i,3) = u(iz_w,1,i,3) - r1*amp_u*cmplx( 0.0,-0.5)
               u(iz_w,1,i,4) = u(iz_w,1,i,4) +    amp_T*cmplx(-0.5, 0.0)
            end if
         end do

         if( ixs .eq. 1 ) then
            u(1,1,1,1) = cmplx(Uo,0.0)
         end if

      end if

      if( i_prob .eq. 4 ) then

         denom = 1.0/sqrt( lambda_x**2 + lambda_z**2 )
         cos_theta = lambda_x*denom
         sin_theta = lambda_z*denom
         lambda = lambda_z*cos_theta
         r1 = 1.0/cos_theta

         iz_w  = nint(zL/lambda) + 1
         iz_w1 = Nz_min+2 - iz_w

         c_i = omega_i/k_w
         amp_u = amplitude*r1*c_i
         amp_T = amplitude*lapse0/m_w

         if( l_root ) then
            u(iz_w ,1,1,1) = u(iz_w,1,1,1) + amp_u*cmplx(0.0,0.5)
            u(iz_w ,1,1,4) = u(iz_w,1,1,4) + amp_T*cmplx(0.5,0.0)
            u(iz_w1,1,1,1) = conjg(u(iz_w,1,1,1))
            u(iz_w1,1,1,4) = conjg(u(iz_w,1,1,4))
         end if

         buoy_fac_x = -sin_theta*buoy_fac
         buoy_fac_z =  cos_theta*buoy_fac

         if( ixs .eq. 1 ) then
            u(1,1,1,1) = cmplx(Uo*cos_theta,0.0)
            u(1,1,1,3) = cmplx(Uo*sin_theta,0.0)
         end if

      end if

      if( i_prob .eq. 5 ) then

         a = 1.0 + shear_ratio
         tanh_a = tanh(a)
         do i=1, 50
            a_old = a
            a = (shear_ratio+1.0)*tanh_a/
     &          (1.0+shear_ratio*(1.0-tanh_a**2))
            tanh_a = tanh(a)
            if( abs(a-a_old) .lt. 1.0e-8 ) goto 20
         end do

20       continue

         scale = (shear*0.5*zL)*1.0/(a - tanh_a)

         dzz = 2.0/float(Nz)

         do k=1, Nz
            zz = -1.0 + float(k-1)*dzz
            random(k,1,1) = scale*( tanh(a*zz) - tanh_a*zz )
         end do

         call rfftf( Nz, random, trigzr )

         u(1,1,1,2) = Nz_inv*cmplx(random(1,1,1),0.0)
         do k=2, kz_max+1
            kr = 2*(k-1)
            ki = kr+1
            kk = Nz_min+2 - k
            u( k,1,1,2) = Nz_inv*cmplx(random(kr,1,1), random(ki,1,1))
            u(kk,1,1,2) = Nz_inv*cmplx(random(kr,1,1),-random(ki,1,1))
         end do

      end if

      deallocate( iseed )

      return
      end


      function e_k( k, i_type )

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


      if( i_type .eq. 0 ) then

         e_k = k**2

      else if( i_type .eq. 1 ) then

         if( k .eq. 0.0 ) then
            e_k = 0.0
         else
            e_k = k**(-5.0/3.0)
         end if

      else if( i_type .eq. 2 ) then

         if( k .eq. 0.0 ) then
            e_k = 0.0
         else
            sum = 0.0
            do n=0, 6
               sum = sum + a42(n)*alog( const1*k )**n
            end do
            e_k = const2*exp( sum )
         end if

      else if( i_type .eq. 3 ) then

         if( k .eq. 0.0 ) then
            e_k = 0.0
         else
            sum = 0.0
            do n=0, 6
               sum = sum + a98(n)*alog( const1*k )**n
            end do
            e_k = const2*exp( sum )
         end if

      else if( i_type .eq. 4 ) then

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

         print *, 'ERROR bad i_type in e_k'
         stop

      end if

      return
      end


      subroutine enforce_symm( u, n_var, i )

      include 'sam.h'

      complex u(Nze,Ny_min,nxp,n_var)

      j=1; k=1
      do n=1, n_var
         u(k,j,i,n) = 0.5*( u(k,j,i,n) + conjg(u(k,j,i,n)) )
      end do

      k=1
      do n=1, n_var
         do j=ky_max+2, Ny_min
            jj = Ny_min+2 - j
            u(k,j,i,n) = conjg(u(k,jj,i,n))
         end do
      end do

      j=1
      do n=1, n_var
         do k=kz_max+2, Nz_min
            kk = Nz_min+2 - k
            u(k,j,i,n) = conjg(u(kk,j,i,n))
         end do
      end do

      do n=1, n_var
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


      subroutine rescale( u, n_var, flct )

      include 'sam.h'
      include 'mpif.h'

      complex u(Nze,Ny_min,nxp,n_var)

      sum1 = 0.0
      do n=1, n_var
         do i=1, nxp
            do j=1, Ny_min
               do k=1, Nz_min
                  sum1 = sum1 + u(k,j,i,n)*conjg(u(k,j,i,n))
               end do
            end do
         end do
      end do

      call mpi_allreduce( sum1, sum1g, 1, mpi_double_precision,
     &                    mpi_sum, mpi_comm_world, ierr )

      scale = flct*sqrt(float(n_var))/sqrt( sum1g  )

      do n=1, n_var
         do i=1, nxp
            do j=1, Ny_min
               do k=1, Nz_min
                  u(k,j,i,n) = u(k,j,i,n)*scale
               end do
            end do
         end do
      end do

      return
      end

      subroutine set_seed( iseed, L_seed, my_seed )

      integer iseed(L_seed)

      iseed(1) = 987654321 + my_seed

      if( L_seed .gt. 1 ) then
         iseed(2) = 123456789
      end if

      i_flip = 1
      do i=3, L_seed
         i_flip = -i_flip
         iseed(i) = iseed(i-1) + i_flip*iseed(i-2)
      end do

      return
      end
