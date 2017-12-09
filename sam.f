      program sam

      include 'sam.h'
      include 'mpif.h'

      complex, allocatable, dimension (:,:,:,:) :: u, r, u_sav
      real,    allocatable, dimension (:,:,:,:) :: uu
      real,    allocatable, dimension (:,:)     :: mean
      real,    allocatable, dimension (:)       :: sample, ek, tk, dk,
     &                                             vt, ek_0, 
     &                                             x_spec, y_spec,
     &                                             z_spec, work
      real jump
      character( 1) ext1
      character( 6) ext
      character(12) fname
      logical       write_params
      integer       i_symm1(4), i_symm2(4), amode

      ierr = 0
      write_params = .true.

c------ Initialize MPI, get myid, numprocs, and test if on root process

      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myid,ierr)
      call mpi_comm_size(mpi_comm_world,numprocs,ierr)

      i_root = 0
      l_root = .false.
      if( myid .eq. i_root ) l_root = .true.

c------ Write a copyright message

      if(l_root) then
         write(6,6)
6        format(/,'© 2017 NorthWest Research Associates, '
     &            'Inc.  All Rights Reserved',/,
     &            'Author: Thomas S. Lund, lund@cora.nwra.com',/)
      end if

c------ Get input parameters

      call input_p( 'input.dat', write_params )

c------ Initialize constants and other parameter values.  Open files.

      call init( write_params, ierr )
      call stop_on_error( ierr, 0 )

      call fft_init( )

      i_symm1 = i_symm

      i_symm2(1) =  i_symm(1)
      i_symm2(2) = -i_symm(1)
      i_symm2(3) = -i_symm(2)
      i_symm2(4) = -i_symm(4)

c ----- Determine index ranges for each process.

      call set_range( Nze, Nx_min )

c ----- Allocate arrays.  The parameter ipad_r is padding used to
c ----- ensure that the r array is large enough to be used as workspace
c ----- in the rhs and write_planes routines.

      k_sq = kx_max**2 + ky_max**2 + kz_max**2
      Lm = int( sqrt(float(k_sq)) + 0.5 ) + 1
      Nu = Lu
      Nr = Lr
      n_words_r_main = Nze*Ny_min*nxp*Lr
      n_words_r_rhs  = Nx_min*Ny_min*nzp*Lr
      n_words_r_wp   = (Nx*Ny)/2 + max(Nx,Ny)*nzp*Lu
      n_words_r = max(n_words_r_main,n_words_r_rhs,n_words_r_wp)
      ipad_r = ceiling((float(n_words_r)/float(n_words_r_main)-1.0)*Lr)
      n_wp = ceiling( float(Nx*Ny)/float(Nze*Ny_min*nxp*2) ) + 1
c      print *, 'myid, ipad_r, n_wp = ', myid, ipad_r, n_wp

      Lw = max( L_params,                               ! read_header
     &          2*Nz_min*Ny_min,                        ! read_field
     &          2*Nze,                                  ! z_trans_b
     &          2*Nze,                                  ! get_mean
     &          (Nx+1)*Ny,                              ! xy_trans_b
     &          2*Nze,                                  ! write_planes
     &          3*Lu,                                   ! vel_max
     &          Nx*Ny,                                  ! xy_trans_f
     &          2*Nze,                                  ! z_trans_f
     &          Lm*12,                                  ! spectra
     &          2*Nz_min*Ny_min,                        ! write_field
     &          Lm*2     )                              ! force

      allocate( u(Nze,Ny_min,nxp,Lu), r(Nze,Ny_min,nxp,Lr+ipad_r),
     &          uu(Nx,Ny,nzp,Lr), work(Lw), mean(Nze,Lu),
     &          sample(Lm), ek(Lm), tk(Lm), dk(Lm), vt(Lm), ek_0(Lm) )
      if( i_strat .eq. 1 ) then
         allocate( u_sav(Nz_min,Ny_min,nxp,3:4) )
      end if
 
c------ Either read or generate a new velocity field

      if( i_restart .eq. 1 ) then
         write(ext,'(i6.6)') nt_restart
         call read_header_p( 'header.'//ext, work )
         call read_field( 'vel.'//ext, u, work )
      else
         call initial_field( u, uu, ierr )
      end if

      call open_files( )

      nt_end = nt_start + n_steps

c ----- Begin time stepping loop.  Similar to the time variable, the 
c ----- time step index, nt, is associated with the start of the time
c ----- step and becomes incremented once the time step is complete.
c ----- Thus a fresh run starts at nt=0.  Note that we include
c ----- one additional time step in order to compute and write 
c ----- statistics in a unified way at the end of the run.  The extra
c ----- time step is aborted after the rhs is called and ouput is done.

      do nt=nt_start, nt_end

         do nrk=1, nrk_max

c         *** Compute the energy in the max shell for the C.-L. model

            if( nrk .eq. 1 .and. i_les .eq. 1 ) then
c            if( i_les .eq. 1 ) then
               call cl_energy( u )
               fac_cl = sqrt( 2*E_cl/float(k_cl) )
            end if

c ----------- Start the process of transforming the velocity to
c ----------- physical space.  Use the second half of r() as workspace.
c ----------- If planes are to be written at this time step, we also 
c ----------- transform the z derivatives and the pressure.  These 
c ----------- quantities are loaded into the first half of r since the 
c ----------- entire array is availble as workspace when nrk=1.  
c ----------- Note that the pressure is at a slightly different time
c ----------- level here, but this is probably ok for the purpose of
c ----------- writing visualization data.
            
            Nu1 = Nu
            if( mod(nt,n_skip_p) .eq. 0 .and. nrk .eq. 1 ) then
               do i=1, nxp
                  do j=1, Ny_min
                     do k=1, Nz_min
c                       r(k,j,i,Lu+1) = the pressure is already here
                        r(k,j,i,Lu+2) = iunit*wave_z*k_z(k)*u(k,j,i,1)
                        r(k,j,i,Lu+3) = iunit*wave_z*k_z(k)*u(k,j,i,2)
                        if( i_strat .eq. 1 ) then
                          r(k,j,i,Lu+4) = iunit*wave_z*k_z(k)*u(k,j,i,4)
                        end if
                     end do
                  end do
               end do
               call z_trans_b( r(1,1,1,Lu+1), r(1,1,1,1), Nu, i_symm2,
     &                         work(1), work(1), work(Nze+1) )
               Nu1 = 2*Nu
               if(l_root) then
                  do n=1, Lu
                     call get_mean( u(1,1,1,n), i_symm(n), 0,
     &                              mean(1,n), jump, work )
                  end do
               end if
            end if

            call z_trans_b( u, r(1,1,1,Lu+1), Nu, i_symm1,
     &                      work(1), work(1), work(Nze+1) )

c ----------- Save a copy of the w velocity and the temperature for
c ----------- stratified cases.  Also compute the lapse rate.

            if( i_strat .eq. 1 ) then
               if( i_prob .eq. 4 ) then
                  do i=1, nxp
                     do j=1, Ny_min
                        do k=1, Nz_min
                           u_sav(k,j,i,3) = cos_theta*u(k,j,i,3) -
     &                                      sin_theta*u(k,j,i,1)
                           u_sav(k,j,i,4) =           u(k,j,i,4)
                        end do
                     end do
                  end do
               else
                  do i=1, nxp
                     do j=1, Ny_min
                        do k=1, Nz_min
                           u_sav(k,j,i,3) = u(k,j,i,3)
                           u_sav(k,j,i,4) = u(k,j,i,4)
                        end do
                     end do
                  end do
               end if
c               if( nrk .eq. 1 .and. i_symm(4) .eq. 0 ) then
c                  if(l_root) then
c                     call get_mean( u(1,1,1,4), i_symm(4), 1,
c     &                              work(1), jump, work(Nze+1) )
c                     lapse = lapse + jump*zL_inv
c      print *, nt, lapse
c                  end if
c                  call mpi_bcast( lapse, 1, mpi_double_precision,
c     &                            i_root, mpi_comm_world, ierr )
c               end if
            end if

c ----------- Update the velocity with the previous right hand side.

            dt_fac2 = zeta(nrk)*dt
            if( nrk .gt. 1 ) then
               do n=1, Lu
                  do i=1, nxp
                     do j=1, Ny_min
                        do k=1, Nz_min
                           u(k,j,i,n) = u(k,j,i,n) + dt_fac2*r(k,j,i,n)
                        end do
                     end do
                  end do
               end do
            end if

c ----------- Finish transforming the velocity to physical space. It
c ----------- will be stored in the uu array.

            call x2z_decomp( r(1,1,1,Lu+1), uu, Nu, uu )

            if( mod(nt,n_skip_p) .eq. 0 .and. nrk .eq. 1 ) then
               call x2z_decomp( r, uu(1,1,1,Lu+1), Nu, uu(1,1,1,Lu+1) )
            end if

            call xy_trans_b( uu, uu, Nu1, work )

c ----------- Write planes data, using r as workspace.

            if( mod(nt,n_skip_p) .eq. 0 .and. nrk .eq. 1 ) then
               call write_planes( uu, r, r(1,1,1,n_wp), r(1,1,1,n_wp),
     &                            mean, work )
            end if

c ----------- Compute the current right hand side and pressure.  The 
c ----------- latter is returned in r(:,:,:,Lu+1).

            call rhs( uu, uu, r, u_sav, work )

c ----------- Compute and write spectra and statistics.

            if( ( mod(nt,n_skip_h) .eq. 0 .or.
     &            nt .eq. nt_end ) .and. nrk .eq. 1 ) then
               call spectra( u, r, sample, ek_0, tk, dk, vt, energy_w,
     &                       div_rms, div_max, n_div, nt, work )
               if(l_root) then
                  write(2,10) nt, time, energy_w
                  write(6,10) nt, time, 
     &                         0.5*(u_var(1)+u_var(2)+u_var(3))
                  write(3,10) nt, time, div_max, u_max(1:Lu)
                  write(7,10) nt, time, dt, cfl_x, cfl_y, cfl_z,
     &                        cfl_vis
                  write(12,10) nt, time, u_var(1:3),
     &                         0.5*(u_var(1)+u_var(2)+u_var(3)),
     &                         u_var(4:Lu)
10                format(i6,1p11e12.4)
                  if( i_prob .eq. 2 ) then
                     write(9,20) time*tfact+t0_cbc, energy_w
20                   format(1p,2e12.4)
                  end if
               end if
            end if

c ----------- Write velocity the field.

            if( (nt .gt. 0 .and. mod(nt,n_skip_v) .eq. 0) .or.
     &           nt .eq. nt_end ) then
               write(ext,'(i6.6)') nt
               if(l_root) then
                  call write_header( 'header.'//ext )
               end if
               call write_field( 'vel.'//ext, u, work )
            end if

c ----------- We can now exit the nt_end time step since all data is written.

            if( nt .eq. nt_end ) goto 50

c ----------- Update the velocity with the current right hand side.
c ----------- Use integrating factors to advance the viscous terms.

            dt_fac1 = gamma(nrk)*dt
            dt_fac = dt_fac1 + dt_fac2
            vis_f = vis*dt_fac
            do n=1, Lu
               if( n .eq. 4 ) vis_f = vis_f*Pr_inv
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
                        vis_fac = exp(-vis_f*rk_sq)
                        u(k,j,i,n) = (u(k,j,i,n) + dt_fac1*r(k,j,i,n))*
     &                               vis_fac
                     end do
                  end do
               end do
            end do

         end do

         time = time + dt

         if( i_force .eq. 1 ) then
           call force( u, ek_0, work(1), work(Lm+1) )
         end if

50       continue

      end do

      call mpi_finalize( ierr )

      stop
      end
