      subroutine fft_init( )

      include 'sam.h'

      logical dump_wavenumbers

      pi = acos(-1.0)

      dump_wavenumbers = .false.

      ncx = Nx/2 + 1
      ncy = Ny/2 + 1 
      ncz = Nz/2 + 1
      
      call rffti( Nx, trigx )
      call cffti( Ny, trigy )
      if( iubc_z .eq. 0 ) then
         call cffti( Nz  , trigz  )
      else 
         call sinti( Nz-1, trigzs )
         call costi( Nz+1, trigzc )
      end if
         
      i_symm = 0
      if( iubc_z .ne. 0 ) then
         i_symm(1) =  1
         i_symm(2) =  1
         i_symm(3) = -1
         i_symm(4) = -1
      end if

      if( n_dealias .eq. 2 ) then
         k_truncate = float(Nx)/3.0
      end if
      if( n_dealias .eq. 3 ) then
         k_truncate = float(Nx)*sqrt(2.0)/3.0
      end if

      k_truncate_sq = k_truncate**2

      Nx2 = Nx/2
      Nz_inv = 1.0/float(Nz)
      NxNy = Nx*Ny
      NxNy_inv = 1.0/float(NxNy)
      NxNyNz = Nx*Ny*Nz
      NxNyNz_inv = 1.0/float(NxNyNz)

      wave_x = 2.0*pi/xL
      do i=1, kx_max+1
         k_x(i) = i-1
      end do

      wave_y = 2.0*pi/yL
      do j=1, ky_max+1
         k_y(j) = j-1
      end do
      do j=ky_max+2, Ny_min
         k_y(j) = j-1-Ny_min
      end do

      if( iubc_z .eq. 0 ) then
         wave_z = 2.0*pi/zL
         do k=1, kz_max+1
            k_z(k) = k-1
         end do
         do k=kz_max+2, Nz_min
            k_z(k) = k-1-Nz_min
         end do
      else
         wave_z = pi/zL
         do k=1, kz_max+1
            k_z(k) = k-1
         end do
      end if

      if( l_root .and. dump_wavenumbers ) then
         open(unit=11,file='kx.out')
         open(unit=12,file='ky.out')
         open(unit=13,file='kz.out')
         do i=1, Nx_min
            write(11,11) i, k_x(i)
11          format(2i5)
         end do
         do j=1, Ny_min
            write(12,11) j, k_y(j)
         end do
         do k=1, Nz_min
            write(13,11) k, k_z(k)
         end do
         close(11)
         close(12)
         close(13)
      end if

      return
      end
