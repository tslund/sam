      subroutine init( write_params, ierr )

      include 'sam.h'

      real lambda
      logical write_params

      pi = acos(-1.0)
      two_pi = 2.0*pi
      iunit = cmplx( 0.0, 1.0 )

      ierr = 0

      NxNy = Nx*Ny
      NxNyNz = NxNy*Nz
      NxNyNzp = NxNy*Nzp
      Nx2 = Nx/2
      Nz_inv = 1.0/float(Nz)
      NxNy_inv = 1.0/float(NxNy)
      NxNyNz_inv = 1.0/float(NxNyNz)

      dx = xL/float(Nx)
      dy = yL/float(Ny)
      dz = zL/float(Nz)
      dx_inv = 1.0/dx
      dy_inv = 1.0/dy
      dz_inv = 1.0/dz
      delta_sq = ( dx*dy*dz )**(2.0/3.0)
      delta_sq_inv = 1.0 / delta_sq

      xL_inv = 1.0/xL
      yL_inv = 1.0/yL
      zL_inv = 1.0/zL

      if( i_strat .eq. 0 ) then
         Lu = 3
         Lr = 6
      else
         Lu = 4
         Lr = 9
      end if

c------ Define the maximum active wavenumbers for wave space data.  This 
c------ depends on the dealiasing scheme.  The Nyquist modes are excluded
c------ when n_dealias=0 since we can not differentiate them.

      Nz2 = Nz
      if( iubc_z .ne. 0 ) Nz2 = 2*Nz
      select case( n_dealias )
      case( 0 )
         kx_max = Nx /2 - 1
         ky_max = Ny /2 - 1
         kz_max = Nz2/2 - 1
      case( 1 )
         kx_max = int((float(Nx )-0.1)/3.0)
         ky_max = int((float(Ny )-0.1)/3.0)
         kz_max = int((float(Nz2)-0.1)/3.0)
      case( 2 )
         kx_max = Nx/2 - 1
         ky_max = int((float(Ny )-0.1)/3.0)
         kz_max = int((float(Nz2)-0.1)/3.0)
      end select

c------ Define the array dimensions for complex wave space data.
c------ Due to conjgate symmetry, only the positive wavenumber
c------ coefficients are stored for x.  Both positive and negative
c------ wavenumber coefficients are needed in the other directions if
c------ they are Fourier.  Only positive wavenumber coefficients are
c------ stored if sines and cosines are used in z.

      Nx_min =   kx_max+1
      Ny_min = 2*ky_max+1
      Nz_min = 2*kz_max+1
      if( iubc_z .ne. 0 ) Nz_min = kz_max+1

c------ Set z index ranges

      if( iubc_z .eq. 0 ) then
         Nze = Nz
         Nz_spec = Nz_min/2 + 1
      else
         Nze = Nz + 1
         Nz_spec = Nz_min + 1
      end if

c------ Set up for the Chollet-Lesieur model

      if( i_les .eq. 1 ) then
         a_cl = 0.145
         b_cl = 5.01
         c_cl = 3.03
         k_cl = int( Nx/3.0 )
      end if

c------ Set up for stratified cases.

      To_inv = 1.0/To
      pr_inv = 1.0/Pr

      scale_h = To/lapse0
      N_sq = grav*lapse0*To_inv

      buoy_fac_x = 0.0
      buoy_fac_z = 0.0
      if( i_strat .eq. 1 ) then
         buoy_fac_z = To_inv*grav
      end if

      time = 0.0
      dt = dt0

c------ Set up for the Runge Kutta time stepping.  Technically both
c ----- RK1 and RK2 schemes are asymptotically unstable for pure advection.
c ----- We set advective cfl limits for these at 1.0, but instabilities
c ----- may still result.

      gamma = 0.0
      zeta  = 0.0
      if( nrk_max .eq. 1 ) then
         gamma(1) = 1.0
         zeta(1)  = 0.0
         vis_cfl = 0.0
         adv_cfl_limit = 1.0
         vis_cfl_limit = 2.0
      end if
      if( nrk_max .eq. 2 ) then
         gamma(1) = 1.0
         gamma(2) = 0.5
         zeta(1)  = 0.0
         zeta(2)  =-0.5
         vis_cfl = 2.0
         adv_cfl_limit = 1.0
         vis_cfl_limit = 2.0
      end if
      if( nrk_max .eq. 3 ) then
         gamma(1) = 8.0/15.0
         gamma(2) = 5.0/12.0
         gamma(3) = 3.0/4.0
         zeta(1)  = 0.0
         zeta(2)  =-17.0/60.0
         zeta(3)  =-5.0/12.0
         vis_cfl = 4.0
         adv_cfl_limit = sqrt(3.0)
         vis_cfl_limit = 2.51
      end if
      do n=1, 4
         beta(n) = 0.5*( zeta(n) + gamma(n) )
      end do

c------ Set up for the Comte Bellot and Corrsin initial condition

      l_cbc = 55.0
      m_cbc = 5.08
      u_inf_cbc = 1000.0
      u_0_cbc = 27.1893
      t0_cbc = 42.0

      tfact = ( l_cbc/(2*pi) / m_cbc ) * ( u_inf_cbc / u_0_cbc )

c   *** Compute a complete set of gravity wave parameters according to the
c   *** problem specification type.  In all cases we specify the horizontal
c   *** wavelength.  For iprob=20 we also specify Uo and omega and compute
c   *** Gam.  For iprob=21 we specify omega and Gam and compute Uo.  For
c   *** iprob=22 we specify Uo and Gam and compute omega.  For
c   *** upward-propagating momentum and energy flux, (wave specified at the
c   *** lower boundary) the vertical wavenumber component (m_w) must be
c   *** negative.  We assume an upward-propagating wave here.  In case it
c   *** is ever needed, a wave with downward energy flux (initiated
c   *** at the upper boundary) would have positive m_w.  The horizontal
c   *** wavenumber is specified via k_w = -2*pi/char_L.  The minus sign
c   *** results in negative horizontal phase speed when char_L is positive.
c   *** This convention is useful for phase-locked simulations (omega=0),
c   *** such as terrain-generated waves where a mean wind in the positive x
c   *** direction results in an equal and opposite (and hence negative) phase
c   *** velocity.  Thus this situation is conveniently specified by taking
c   *** both Uo and char_L to be positive.  For a similar situation with
c   *** wind in the negative x direction, one would specify both Uo and
c   *** char_L as negative.  Note that Gam = lambda_x/lambda_z = abs(m_w/k_w),
c   *** which is a positive definite quantity.

c   *** Compute the GW parameters using the state at the first solution
c   *** point (k=2).

      if( abs(i_prob) .eq. 4 ) then

         k_w = -two_pi/lambda_x

         if( i_gw_type .eq. 3 ) then
            Gam = sqrt( N_sq/(omega-k_w*Uo)**2 - 1.0 )
            i = index_param( 'gam', labels, n_inputs )
            values(i) = Gam
         end if

         m_w = -Gam*abs(k_w)    ! negative for upward propagating energy flux
         omega_i = sqrt( N_sq/( 1.0 + Gam**2 ) )

         if( i_gw_type .eq. 1 ) then
            Uo = ( omega - omega_i )/k_w
            i = index_param( 'uo', labels, n_inputs )
            values(i) = Uo
         end if

         if( i_gw_type .eq. 2 ) then
            omega = omega_i + k_w*Uo
            i = index_param( 'omega', labels, n_inputs )
            values(i) = omega
         end if

         lambda_z = two_pi/abs(m_w)
         char_u = sqrt(N_sq/m_w**2)

         if( l_root .and. write_params ) then
            write(6,16) k_w, m_w, sqrt(N_sq), omega_i, omega,
     &                  Gam, Uo, char_u
16          format('k_w, m_w = ', 1p,2e16.8,/,
     &             'N        = ', 1p, e16.8,/,
     &             'omega_i  = ', 1p, e16.8,/,
     &             'omega    = ', 1p, e16.8,/,
     &             'Gam      = ', 1p, e16.8,/,
     &             'Uo       = ', 1p, e16.8,/,
     &             'char_u   = ', 1p, e16.8,/ )
         end if

         if( i_prob .eq. -4 ) then
            xL_old = xL
            xL = nint(xL/lambda_x)*lambda_x
            zL_old = zL
            zL = nint(zL/lambda_z)*lambda_z
            if( l_root .and. write_params .and.
     &          abs(xL-xL_old) .gt. 1.0e-12 ) then
               print *, 'WARNING: The computational box size xL ',
     &                  'was recomputed to be an exact integer'
               print *, 'multiple of the input GW parameter lambda_x'
               print *, 'The old and new xL are ', xL_old, xL
            end if
            if( l_root .and. write_params .and. 
     &          abs(zL-zL_old) .gt. 1.0e-12 ) then
               print *, 'INFO: The computational box size zL ',
     &                  'was recomputed to be an exact integer'
               print *, 'multiple of the computed GW parameter lambda_z'
               print *, 'The old and new zL are ', zL_old, zL
            end if
         end if

         lambda = lambda_z*lambda_x/sqrt(lambda_x**2+lambda_z**2)
         if( i_prob .eq. 4 ) then
            zL_old = zL
            zL = nint(zL/lambda)*lambda
            if( l_root .and. write_params .and. 
     &          abs(zL-zL_old) .gt. 1.0e-12 ) then
               print *, 'INFO: The computational box size zL ',
     &                  'was recomputed to be an exact integer'
               print *, 'multiple of the computed GW parameter lambda'
               print *, 'The old and new zL are ', zL_old, zL
            end if
         end if

c   *** Compute non-dimensional parameters.

         Re = char_u*lambda_z/(vis+1.0e-20)
         Fr = char_u/sqrt(grav*lambda_z)

         if( l_root .and. write_params ) then
            write(6,15) lambda_z/scale_h, Re, Pr, Fr
15          format(/,'lambda_z/scale_h = ',1p,e15.5,/,
     &               'Reynolds number  = ',1p,e15.5,/,
     &               'Prandtl  number  = ',1p,e15.5,/,
     &               'Froude   number  = ',1p,e15.5,/ )
         end if

      end if

      return
      end
