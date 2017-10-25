      subroutine init( )

      include 'sam.h'

      pi = acos(-1.0)
      two_pi = 2.0*pi
      iunit = cmplx( 0.0, 1.0 )

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

      if(l_root) print *, 'kx_max, ky_max, kz_max = ', 
     &                     kx_max, ky_max, kz_max
      if(l_root) print *, 'Nx_min, Ny_min, Nz_min, Nze = ',
     &                     Nx_min, Ny_min, Nz_min, Nze

c------ Set up for the Chollet-Lesieur model

      if( i_les .eq. 1 ) then
         a_cl = 0.145
         b_cl = 5.01
         c_cl = 3.03
         k_cl = int( Nx/3.0 )
      end if

      time = 0.0
      dt = dt0

c------ Set up for the Runge Kutta time stepping

      gamma = 0.0
      zeta  = 0.0
      if( nrk_max .eq. 1 ) then
         gamma(1) = 1.0
         zeta(1)  = 0.0
         vis_cfl = 0.0
      end if
      if( nrk_max .eq. 2 ) then
         gamma(1) = 1.0
         gamma(2) = 0.5
         zeta(1)  = 0.0
         zeta(2)  =-0.5
         vis_cfl = 2.0
      end if
      if( nrk_max .eq. 3 ) then
         gamma(1) = 8.0/15.0
         gamma(2) = 5.0/12.0
         gamma(3) = 3.0/4.0
         zeta(1)  = 0.0
         zeta(2)  =-17.0/60.0
         zeta(3)  =-5.0/12.0
         vis_cfl = 4.0
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

      return
      end
