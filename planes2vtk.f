      program planes2vtk

c  © 2002 – 2014 NorthWest Research Associates, Inc. All Rights Reserved
c    Author: Thomas S. Lund, lund@cora.nwra.com

      include 'sam.h'

      real,    allocatable, dimension (:,:,:) :: u, du
      real,    allocatable, dimension (:)     :: work
      real(4), allocatable, dimension (:,:,:) :: xz_rot
      real(4), allocatable, dimension (:,:)   :: u_out
      real(4), allocatable, dimension (:)     :: x, y, z
      real          values(L_params)
      logical        fixed(L_params), write_params
      character(12) labels(L_params), fname
      character( 1) newline, ext1
      character( 4) ext
      character( 8) Nx_str, Ny_str, Nz_str, Nxy_str, Nxz_str, Nyz_str
      character(11) time_string, label(8)

      newline = achar(10)
      l_root = .true.
      write_params = .false.

c   *** Read input parameters and initialize stuff

      call input( 'input.dat', labels, values, fixed, write_params )

      call init( write_params, ierr )
      if( ierr .ne. 0 ) stop

      call fft_init( )

      print *, 'enter the time step number start, stop, and skip factor'
      read(5,*) n_start, n_stop, n_skip
      n_skip = n_skip*n_skip_p

      i_rotate = 0
      if( i_prob .eq. 4 ) then
         print *, 'rotate the box into earth-fixed coordinates? [0,1]'
         read(5,*) i_rotate
      end if

c   *** Allocate arrays

      Nxp1 = Nx + 1
      Nyp1 = Ny + 1
      Nzp1 = Nz + 1

      allocate( u(Nx,Ny,2*Lu), du(Nx,Ny,6), u_out(Nxp1,Nyp1), x(Nxp1),
     &          y(Nyp1), z(Nzp1), work(max(Nx,Ny,Nz)) )

      if( i_rotate .eq. 1 ) allocate( xz_rot(3,Nxp1,Nzp1) )

c   *** Define x and y mesh points

      do i=1, Nx+1
         x(i) = -0.5*xL + float(i-1)*dx
      end do
      do j=1, Ny+1
         y(j) = -0.5*yL + float(j-1)*dy
      end do
      do k=1, Nz+1
         z(k) = -0.5*zL + float(k-1)*dz
      end do

      if( i_rotate .eq. 1 ) then
         denom = 1.0/sqrt(lambda_x**2+lambda_z**2)
         cos_theta = lambda_x*denom
         sin_theta = lambda_z*denom
         do k=1, Nzp1
            do i=1, Nxp1
               xz_rot(1,i,k) =  x(i)*cos_theta + z(k)*sin_theta
               xz_rot(2,i,k) =  0.0
               xz_rot(3,i,k) = -x(i)*sin_theta + z(k)*cos_theta
            end do
         end do
      end if

c   *** Define mesh point number character strings.

      write(Nx_str ,'(i8)') Nx+1
      write(Ny_str ,'(i8)') Ny+1
      write(Nz_str ,'(i8)') Nz+1
      write(Nxy_str,'(i8)') Nxp1*Nyp1
      write(Nxz_str,'(i8)') Nxp1*Nzp1
      write(Nyz_str,'(i8)') Nyp1*Nzp1

      label(1) = "u'"
      label(2) = "v'"
      label(3) = "w'"
      if( i_strat .eq. 1 ) then
         label(4) = "T'"
      end if
      label(Lu+1) = "p'"

      open(unit=3,file='time.out')

c   *** Loop over the xy frames, reading and writing the data

      u_out = 0.0

      read(3,*,end=99)

      do m=1, n_xy_planes

         write(ext1,'(i1)') m

         open(unit=10,file='xy'//ext1//'.out',form='unformatted',
     &        access='direct',recl=Nx*Ny*2*Lu*8,action='read')

         do nf=1, 10000

            read(3,*,end=20) nt, time

            if( nt .ge. n_start .and. mod(nt,n_skip) .eq. 0 ) then

               n_frame = nt/n_skip_p
               write(ext,'(i4.4)') n_frame
               write(time_string,'(1p,e11.4)') time

               read(10,rec=n_frame+1) u

               call deriv1( u(1,1,1), du(1,1,1), Nx, Ny, trigx ,
     &                      k_x, kx_max, wave_x, work )
               call deriv2( u(1,1,1), du(1,1,2), Nx, Ny, trigyr,
     &                      k_y, ky_max, wave_y, work )
               call deriv1( u(1,1,2), du(1,1,3), Nx, Ny, trigx ,
     &                      k_x, kx_max, wave_x, work )
               call deriv2( u(1,1,2), du(1,1,4), Nx, Ny, trigyr,
     &                      k_y, ky_max, wave_y, work )
               call deriv1( u(1,1,3), du(1,1,5), Nx, Ny, trigx ,
     &                      k_x, kx_max, wave_x, work )
               call deriv2( u(1,1,3), du(1,1,6), Nx, Ny, trigyr,
     &                      k_y, ky_max, wave_y, work )

               fname = 'xy'//ext1//'_'//ext//'.vtk'
               open(unit=16,file=fname,form='binary',
     &              convert='big_endian',access='stream',action='write')

               write(16) "# vtk DataFile Version 2.0"//newline//
     &                   "written by planes2vtk t ="//time_string//
     &                   newline//"BINARY"//newline
               write(16) "DATASET RECTILINEAR_GRID"//newline//
     &                   "DIMENSIONS"//Nx_str//Ny_str//" 1 "//newline
               write(16) "X_COORDINATES "//Nx_str//" float"//newline,
     &                   x, newline
               write(16) "Y_COORDINATES "//Ny_str//" float"//newline,
     &                   y, newline
               write(16) "Z_COORDINATES 1 float"//newline,
     &                   z(k_xy_plane(m)), newline
               write(16) "POINT_DATA"//Nxy_str//newline

               do n=1, Lu+1
                  do j=1, Ny
                     do i=1, Nx
                        u_out(i,j) = u(i,j,n)
                     end do
                  end do
                  call periodicity4( u_out, Nx, Ny )
                  write(16) "Scalars "//trim(label(n))//" float"//
     &                      newline//"LOOKUP_TABLE default"//newline,
     &                      u_out, newline
               end do

               do j=1, Ny
                  do i=1, Nx
                     u_out(i,j) = du(i,j,3) - du(i,j,2)
                  end do
               end do
               call periodicity4( u_out, Nx, Ny )
               write(16) "Scalars vort_z float"//
     &                   newline//"LOOKUP_TABLE default"//newline,
     &                   u_out, newline

c   wx =  (dwdy - dvdz)
c   wy = -(dwdx - dudz)
c   wz =  (dvdx - dudy)

               do j=1, Ny
                  do i=1, Nx
                     w_x =  (du(i,j,6) -  u(i,j,Lu+3))
                     w_y = -(du(i,j,5) -  u(i,j,Lu+2))
                     w_z =  (du(i,j,3) - du(i,j,   2))
                     u_out(i,j) = sqrt( w_x**2 + w_y**2 + w_z**2 )
                  end do
               end do
               call periodicity4( u_out, Nx, Ny )
               write(16) "Scalars vort_mag float"//
     &                   newline//"LOOKUP_TABLE default"//newline,
     &                   u_out, newline

            end if

            close(16)

         end do

         close(10)

      end do

20    continue

c   *** Loop over the xz frames, reading and writing the data

      deallocate( u, du, u_out )
      allocate( u(Nx,Nz,2*Lu), du(Nx,Nz,6), u_out(Nxp1,Nzp1) )

      rewind(3)
      read(3,*,end=99)

      do m=1, n_xz_planes

         write(ext1,'(i1)') m

         open(unit=10,file='xz'//ext1//'.out',form='unformatted',
     &        access='direct',recl=Nx*Nz*2*Lu*8,action='read')

         do nf=1, 10000

            read(3,*,end=30) nt, time

            if( nt .ge. n_start .and. mod(nt,n_skip) .eq. 0 ) then

               n_frame = nt/n_skip_p
               write(ext,'(i4.4)') n_frame
               write(time_string,'(1p,e11.4)') time

               read(10,rec=n_frame+1) u

               call deriv1( u(1,1,1), du(1,1,1), Nx, Nz, trigx ,
     &                      k_x, kx_max, wave_x, work )
               call deriv2( u(1,1,1), du(1,1,2), Nx, Nz, trigzr,
     &                      k_z, kz_max, wave_z, work )
               call deriv1( u(1,1,2), du(1,1,3), Nx, Nz, trigx ,
     &                      k_x, kx_max, wave_x, work )
               call deriv2( u(1,1,2), du(1,1,4), Nx, Nz, trigzr,
     &                      k_z, kz_max, wave_z, work )
               call deriv1( u(1,1,3), du(1,1,5), Nx, Nz, trigx ,
     &                      k_x, kx_max, wave_x, work )
               call deriv2( u(1,1,3), du(1,1,6), Nx, Nz, trigzr,
     &                      k_z, kz_max, wave_z, work )

               fname = 'xz'//ext1//'_'//ext//'.vtk'
               open(unit=16,file=fname,form='binary',
     &              convert='big_endian',access='stream',action='write')

               write(16) "# vtk DataFile Version 2.0"//newline//
     &                   "written by planes2vtk t ="//time_string//
     &                   newline//"BINARY"//newline
               if( i_rotate .eq. 1 ) then
                  write(16) "DATASET STRUCTURED_GRID"//newline//
     &                   "DIMENSIONS"//Nx_str//" 1 "//Nz_str//newline//
     &                   "POINTS"//Nxz_str//" float"//newline
                  write(16) xz_rot, newline
               else
                  write(16) "DATASET RECTILINEAR_GRID"//newline//
     &                      "DIMENSIONS"//Nx_str//" 1 "//Nz_str//newline
                  write(16) "X_COORDINATES "//Nx_str//" float"//newline,
     &                      x, newline
                  write(16) "Y_COORDINATES 1 float"//newline,
     &                      y(j_xz_plane(m)), newline
                  write(16) "Z_COORDINATES "//Nz_str//" float"//newline,
     &                      z, newline
               end if
               write(16) "POINT_DATA"//Nxz_str//newline

               do n=1, Lu+1
                  do k=1, Nz
                     do i=1, Nx
                        u_out(i,k) = u(i,k,n)
c      if( i .eq. 60 ) print *, k, u(i,k,n)
                     end do
                  end do
                  call periodicity4( u_out, Nx, Nz )
                  write(16) "Scalars "//trim(label(n))//" float"//
     &                      newline//"LOOKUP_TABLE default"//newline,
     &                      u_out, newline
               end do

               do k=1, Nz
                  do i=1, Nx
                     u_out(i,k) = du(i,k,2) - du(i,k,5)
                  end do
               end do
               call periodicity4( u_out, Nx, Nz )
               write(16) "Scalars vort_y float"//
     &                   newline//"LOOKUP_TABLE default"//newline,
     &                   u_out, newline

c   wx =  (dwdy - dvdz)
c   wy = -(dwdx - dudz)
c   wz =  (dvdx - dudy)

               do k=1, Nz
                  do i=1, Nx
                     w_x =  ( u(i,k,Lu+3) - du(i,k,   4))
                     w_y = -(du(i,k,   5) - du(i,k,   2))
                     w_z =  (du(i,k,   3) -  u(i,k,Lu+2))
                     u_out(i,k) = sqrt( w_x**2 + w_y**2 + w_z**2 )
                  end do
               end do
               call periodicity4( u_out, Nx, Nz )
               write(16) "Scalars vort_mag float"//
     &                   newline//"LOOKUP_TABLE default"//newline,
     &                   u_out, newline

            end if

            close(16)

         end do

         close(10)

      end do

30    continue

c   *** Loop over the yz frames, reading and writing the data

      deallocate( u, du, u_out )
      allocate( u(Ny,Nz,2*Lu), du(Ny,Nz,6), u_out(Nyp1,Nzp1) )

      rewind(3)
      read(3,*,end=99)

      do m=1, n_yz_planes

         write(ext1,'(i1)') m

         open(unit=10,file='yz'//ext1//'.out',form='unformatted',
     &        access='direct',recl=Ny*Nz*2*Lu*8,action='read')

         do nf=1, 10000

            read(3,*,end=40) nt, time

            if( nt .ge. n_start .and. mod(nt,n_skip) .eq. 0 ) then

               n_frame = nt/n_skip_p
               write(ext,'(i4.4)') n_frame
               write(time_string,'(1p,e11.4)') time

               read(10,rec=n_frame+1) u

               call deriv1( u(1,1,1), du(1,1,1), Ny, Nz, trigyr,
     &                      k_y, ky_max, wave_y, work )
               call deriv2( u(1,1,1), du(1,1,2), Ny, Nz, trigzr,
     &                      k_z, kz_max, wave_z, work )
               call deriv1( u(1,1,2), du(1,1,3), Ny, Nz, trigyr,
     &                      k_y, ky_max, wave_y, work )
               call deriv2( u(1,1,2), du(1,1,4), Ny, Nz, trigzr,
     &                      k_z, kz_max, wave_z, work )
               call deriv1( u(1,1,3), du(1,1,5), Ny, Nz, trigyr,
     &                      k_y, ky_max, wave_y, work )
               call deriv2( u(1,1,3), du(1,1,6), Ny, Nz, trigzr,
     &                      k_z, kz_max, wave_z, work )

               fname = 'yz'//ext1//'_'//ext//'.vtk'
               open(unit=16,file=fname,form='binary',
     &              convert='big_endian',access='stream',action='write')

               write(16) "# vtk DataFile Version 2.0"//newline//
     &                   "written by planes2vtk t ="//time_string//
     &                   newline//"BINARY"//newline
               write(16) "DATASET RECTILINEAR_GRID"//newline//
     &                   "DIMENSIONS 1 "//Ny_str//Nz_str//newline
               write(16) "X_COORDINATES 1 float"//newline,
     &                   x(i_yz_plane(m)), newline
               write(16) "Y_COORDINATES "//Ny_str//" float"//newline,
     &                   y, newline
               write(16) "Z_COORDINATES "//Nz_str//" float"//newline,
     &                   z, newline
               write(16) "POINT_DATA"//Nyz_str//newline

               do n=1, Lu+1
                  do k=1, Nz
                     do j=1, Ny
                        u_out(j,k) = u(j,k,n)
                     end do
                  end do
                  call periodicity4( u_out, Ny, Nz )
                  write(16) "Scalars "//trim(label(n))//" float"//
     &                      newline//"LOOKUP_TABLE default"//newline,
     &                      u_out, newline
               end do

               do k=1, Nz
                  do j=1, Ny
                     u_out(j,k) = du(j,k,5) - du(j,k,4)
                  end do
               end do
               call periodicity4( u_out, Ny, Nz )
               write(16) "Scalars vort_x float"//
     &                   newline//"LOOKUP_TABLE default"//newline,
     &                   u_out, newline

c   wx =  (dwdy - dvdz)
c   wy = -(dwdx - dudz)
c   wz =  (dvdx - dudy)

               do k=1, Nz
                  do j=1, Ny
                     w_x =  (du(j,k,   5) - du(j,k,   4))
                     w_y = -( u(j,k,Lu+3) - du(j,k,   2))
                     w_z =  ( u(j,k,Lu+2) - du(j,k,   1))
                     u_out(j,k) = sqrt( w_x**2 + w_y**2 + w_z**2 )
                  end do
               end do
               call periodicity4( u_out, Ny, Nz )
               write(16) "Scalars vort_mag float"//
     &                   newline//"LOOKUP_TABLE default"//newline,
     &                   u_out, newline

            end if

            close(16)

         end do

         close(10)

      end do

40    continue

99    stop
      end

      subroutine periodicity4( u_out, Nx, Ny )

      real(4) u_out(Nx+1,Ny+1)

      Nxp1 = Nx+1
      Nyp1 = Ny + 1

      do j=1, Ny
         u_out(Nxp1,j) = u_out(1,j)
      end do

      do i=1, Nxp1
         u_out(i,Nyp1) = u_out(i,1)
      end do

      return
      end
