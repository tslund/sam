      program remesh

c  © 2002 – 2014 NorthWest Research Associates, Inc. All Rights Reserved
c    Author: Thomas S. Lund, lund@cora.nwra.com

      include 'sam.h'

      complex, allocatable, dimension (:,:)   :: u_old, u_new
      real          work(L_params)
      logical       write_params
      character( 6) ext
      integer(8) isize_found, isize_expect, n_rec


      l_root = .true.
      write_params = .false.

c   *** Read input parameters for the new mesh

      call input( 'input.dat', write_params, ierr )
      if( ierr .ne. 0 ) stop

      Nx_new = Nx
      Ny_new = Ny
      Nz_new = Nz

c   *** Get the old mesh dimensions from the old header file

      open( unit=2, file='header.in' )

      call grab_input( 2, 'Nx', dummy, Nx_old, ierr )
      call grab_input( 2, 'Ny', dummy, Ny_old, ierr )
      call grab_input( 2, 'Nz', dummy, Nz_old, ierr )

      close(2)

c   *** Get the minimal dimensions for the old mesh

      Nx = Nx_old
      Ny = Ny_old
      Nz = Nz_old

      call init( write_params, ierr )
      if( ierr .ne. 0 ) stop

      Nx_min_old = Nx_min
      Ny_min_old = Ny_min
      Nz_min_old = Nz_min
      kx_max_old = kx_max
      ky_max_old = ky_max
      kz_max_old = kz_max

c   *** Get the minimal dimensions for the new mesh

      Nx = Nx_new
      Ny = Ny_new
      Nz = Nz_new

      call init( write_params, ierr )
      if( ierr .ne. 0 ) stop

      Nx_min_new = Nx_min
      Ny_min_new = Ny_min
      Nz_min_new = Nz_min
      kx_max_new = kx_max
      ky_max_new = ky_max
      kz_max_new = kz_max

c   *** Allocate arrays

      allocate( u_old(Nz_min_old,Ny_min_old),
     &          u_new(Nz_min_new,Ny_min_new)  )

c   *** Read the entire old header file

      i_nx = index_param( 'nx', labels, n_inputs )
      values(i_nx) = Nx_old
      i_ny = index_param( 'ny', labels, n_inputs )
      values(i_ny) = Ny_old
      i_nz = index_param( 'nz', labels, n_inputs )
      values(i_nz) = Nz_old

      call read_header( 'header.in', work(1), ierr )
      if( ierr .ne. 0 ) stop
      if( nt .ne. nt_restart ) then
         print *, 'ERROR nt_restart from input.dat does not match nt ',
     &            'from header.in'
         stop
      end if

      values(i_nx) = Nx_new
      values(i_ny) = Ny_new
      values(i_nz) = Nz_new

c   *** compute grid index limits

      ie1 = min(Nx_min_old,Nx_min_new)
      ib2 = ie1+1

      je1 = min(ky_max_old+1,ky_max_new+1)
      jb2 = je1+1
      je2 = Ny_min_new - ky_max_old
      jb3 = max(je1,je2)+1
      js3 = Ny_min_old - Ny_min_new

      ke1 = min(kz_max_old+1,kz_max_new+1)
      kb2 = ke1+1
      ke2 = Nz_min_new - kz_max_old
      kb3 = max(ke1,ke2)+1
      ks3 = Nz_min_old - Nz_min_new

c   *** Open data files

      write(ext,'(i6.6)') nt_restart

      open(unit=1,file='vel.in'  ,form='unformatted',access='direct',
     &     recl=Nz_min_old*Ny_min_old*16,action='read')

      open(unit=2,file='vel.'//ext,form='unformatted',access='direct',
     &     recl=Nz_min_new*Ny_min_new*16,action='write')

c   *** Make sure the old data file is of the correct size

      isize_expect = Nz_min_old*Ny_min_old*Nx_min_old*Lu*16

      isize_found  = isize_expect
      call get_fsize( 'vel.in', Lu*8, isize_found, n_rec )

      if( isize_found .ne. isize_expect ) then
         print *, 'ERROR: Vel file is not of correct size'
         print *, 'isize_expect, isize_found = ',
     &             isize_expect, isize_found
         stop
      end if

c   *** Read the old data, resize it, and write it out

      do n=1, Lu
         do i=1, ie1
            i_rec_old = (n-1)*Nx_min_old + i
            i_rec_new = (n-1)*Nx_min_new + i
            read(1,rec=i_rec_old) u_old
            do j=1, je1
               do k=1, ke1
                  u_new(k,j) = u_old(k,j)
               end do
               do k=kb2, ke2
                  u_new(k,j) = cmplx(0.0,0.0)
               end do
               do k=kb3, Nz_min_new
                  u_new(k,j) = u_old(k+ks3,j)
               end do
            end do
            do j=jb2, je2
               do k=1, Nz_min_new
                  u_new(k,j) = cmplx(0.0,0.0)
               end do
            end do
            do j=jb3, Ny_min_new
               do k=1, ke1
                  u_new(k,j) = u_old(k,j+js3)
               end do
               do k=kb2, ke2
                  u_new(k,j) = cmplx(0.0,0.0)
               end do
               do k=kb3, Nz_min_new
                  u_new(k,j) = u_old(k+ks3,j+js3)
               end do
            end do
            write(2,rec=i_rec_new) u_new
         end do
         do i=ib2, Nx_min_new
            i_rec_new = (n-1)*Nx_min_new + i
            u_new = cmplx(0.0,0.0)
            write(2,rec=i_rec_new) u_new
         end do
      end do

      close(1)
      close(2)

      call write_header( 'header.'//ext, labels, values )

      stop
      end
