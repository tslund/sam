      program sam

      include 'sam.h'
      include 'mpif.h'

      complex, allocatable, dimension (:,:,:,:) :: u, r
      real,    allocatable, dimension (:,:,:,:) :: uu
      real,    allocatable, dimension (:)       :: sample, ek, tk, dk,
     &                                             vt, ek_0,
     &                                             x_spec, y_spec,
     &                                             z_spec, work
      character( 6) ext
      character(12) labels(L_params)
      real          values(L_params)
      logical        fixed(L_params)
      integer       i_symm1(4)

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

      call input_p( 'input.dat', labels, values, fixed )

c------ Initialize constants and other parameter values.  Open files.

      call init( )

      call fft_init( )

      i_symm1 = i_symm

c ----- Determine index ranges for each process.

      call set_range( Nze, Nx_min )

c ----- Allocate arrays.  The parameter ipad_r is padding used to
c ----- ensure that the r array is large enough to be used as workspace
c ----- in the rhs routine.

      k_sq = kx_max**2 + ky_max**2 + kz_max**2
      Lm = int( sqrt(float(k_sq)) + 0.5 ) + 1
      Nu = Lu
      Nr = Lr
      n_words_r_main = Nze*Ny_min*nxp*Lr
      n_words_r_rhs  = Nx_min*Ny_min*nzp*Lr
      n_words_r = max(n_words_r_main,n_words_r_rhs)
      ipad_r = ceiling((float(n_words_r)/float(n_words_r_main)-1.0)*Lr)

      Lw = max( 2*Nze,                                  ! main
     &          2*L_params,                             ! read_header
     &          (Nx+1)*Ny,                              ! xy_trans_f
     &          Ny*Nx,                                  ! xy_trans_b
     &          2*Nze,                                  ! z_trans_f
     &          2*Nze,                                  ! z_trans_b
     &          Lm*12      )                            ! spectra

      allocate( u(Nze,Ny_min,nxp,Lu), r(Nze,Ny_min,nxp,Lr+ipad_r),
     &          uu(Nx,Ny,nzp,Lr), work(Lw), sample(Lm), ek(Lm), tk(Lm),
     &          dk(Lm), vt(Lm), ek_0(Lm) )
 
c------ Either read or generate a new velocity field

      if( i_restart .eq. 1 ) then
         write(ext,'(i6.6)') nt_restart
         call read_header( 'header.'//ext, labels, values, fixed,
     &                     work(1), work(n_params+1) )
         call read_field( 'vel.'//ext, u, work )
      else
         call initial_field( u, uu, ierr )
      end if

      if( l_root ) call open_files( )

      nt_end = nt_start + n_steps

c ----- Begin time stepping loop.  We include one additional time step
c ----- in order to compute statistics and spectra at the end of the run.

      do nt=nt_start+1, nt_end+1

         do nrk=1, nrk_max

            dt_fac1 = gamma(nrk)*dt
            dt_fac2 = zeta( nrk)*dt
            dt_fac = dt_fac1 + dt_fac2

c         *** Compute the energy in the max shell for the C.-L. model

            if( nrk .eq. 1 .and. i_les .eq. 1 ) then
c            if( i_les .eq. 1 ) then
               call cl_energy( u )
               fac_cl = sqrt( 2*E_cl/float(k_cl) )
            end if

c ----------- Start the process of transforming the velocity to
c ----------- physical space.  Use the second half of r() as workspace

            
            call z_trans_b( u, r(1,1,1,Lu+1), Nu, i_symm1,
     &                      work(1), work(1), work(Nze+1) )

c ----------- Update the velocity with the previous right hand side.

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

            call xy_trans_b( uu, uu, Nu, work )

c ----------- Compute the current right hand side and pressure.  The 
c ----------- latter is returned in r(:,:,:,Lu+1).

            call rhs( uu, uu, r, work )

c ----------- Compute and write spectra and statistics.

            if( ( mod(nt-1,n_skip_h) .eq. 0 .or.
     &            nt-1 .eq. nt_end ) .and. nrk .eq. 1 ) then
               call spectra( u, r, sample, ek_0, tk, dk, vt, energy_w,
     &                       div_rms, div_max, n_div, nt-1, work )
               if(l_root) then
                  write(3,10) nt-1, time, dt, div_max, energy_w,
     &                        u_max(1:3)
10                format(i3,1p,7e11.3)
                  write(6,10) nt-1, time, energy_w
                  if( i_prob .eq. 2 ) then
                     write(9,20) time*tfact+t0_cbc, energy_w
20                   format(1p,2e12.4)
                  end if
               end if
               if( nt-1 .eq. nt_end ) goto 50
            end if

c ----------- Update the velocity with the current right hand side.
c ----------- Use integrating factors to advance the viscous terms.

            do n=1, Lu
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
                        vis_fac = exp(-vis*rk_sq*dt_fac)
c                        vis_fac = 1.0
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

c -------- Write velocity field

         if( mod(nt,n_skip_v) .eq. 0 .or. nt .eq. nt_end ) then
            write(ext,'(i6.6)') nt
            if(l_root) then
               call write_header( 'header.'//ext, labels, values )
            end if
            call write_field( 'vel.'//ext, u, work )
         end if

      end do

50    call mpi_finalize( ierr )

      stop
      end
