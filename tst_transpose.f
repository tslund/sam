      program tst_transpose

      include 'sam.h'
      include 'mpif.h'

      complex, allocatable, dimension (:,:,:,:) :: u, r
      complex(4), allocatable, dimension (:,:,:) :: u_out
      real,    allocatable, dimension (:,:,:,:) :: uu
      real,    allocatable, dimension (:)       :: sample, ek, tk, dk,
     &                                             vt, ek_0,
     &                                             x_spec, y_spec,
     &                                             z_spec, work
      complex diff
      character( 2) nx_str, ny_str, nz_str, ext2
      character( 6) ext
      character(12) labels(L_params)
      real          values(L_params)
      logical        fixed(L_params), write_params
      integer       i_symm1(4)
      integer(kind=mpi_offset_kind) :: offset, i_recl_out_8
      integer amode, status(mpi_status_size)

      write_params = .false.

c------ Initialize MPI, get myid, numprocs, and test if on root process

      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myid,ierr)
      call mpi_comm_size(mpi_comm_world,numprocs,ierr)

      i_root = 0
      l_root = .false.
      if( myid .eq. i_root ) l_root = .true.

c------ Get input parameters

      call input_p( 'input.dat', labels, values, fixed, write_params )

c------ Initialize constants and other parameter values.  Open files.

      call init( )

      call fft_init( )

      i_symm1 = i_symm

      if( l_root ) call open_files( )

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
      print *, 'myid, frac_r, i_pad_r = ', myid,
     &          float(n_words_r)/float(n_words_r_main), ipad_r

      Lw = max( 2*Nze,                                  ! main
     &          2*L_params,                             ! read_header
     &          Nz_min*Ny_min*Lu,                       ! initial_field
     &          (Nx+1)*Ny,                              ! xy_trans_f
     &          Ny*Nx,                                  ! xy_trans_b
     &          2*Nze,                                  ! z_trans_f
     &          2*Nze,                                  ! z_trans_b
     &          Lm*12      )                            ! spectra

      allocate( u(Nze,Ny_min,nxp,Lu), r(Nze,Ny_min,nxp,Lr+ipad_r),
     &          uu(Nx,Ny,nzp,Lr), work(Lw), sample(Lm), ek(Lm), tk(Lm),
     &          dk(Lm), vt(Lm), ek_0(Lm) )
      allocate( u_out(Nz_min+1,Ny_min+1,Nx_min) )
 
c------ Either read or generate a new velocity field

      if( i_restart .eq. 1 ) then
         write(ext,'(i6.6)') nt_restart
         call read_header( 'header.'//ext, labels, values, fixed,
     &                     work(1), work(n_params+1) )
         call read_field( 'vel.'//ext, u, work )
      else
         call initial_field( u, work, ierr )
      end if

      r = cmplx(0.0,0.0)

      call spectra( u, r, sample, ek_0, tk, dk, vt, energy_w,
     &              div_rms, div_max, n_div, nt, work )

      write(3,10) nt, time, dt, div_max, energy_w, u_avg(1:Lu)
10    format( i3, 1p,7e11.3 )

c ----- Do a round trip fft and write the data out

      call z_trans_b( u, r(1,1,1,Lu+1), Nu, i_symm1,
     &                work(1), work(1), work(Nze+1) )

      call x2z_decomp( r(1,1,1,Lu+1), uu, Nu, uu )

      call xy_trans_b( uu, uu, Nu, work )

      open(unit=17,file='sam_vel.out',form='unformatted',
     &     access='direct',recl=Nx*Ny*8,action='write')
      do n=1, 3
         do k=1, nzp
            kk = izs + k-1
            ind = (n-1)*Nz + kk
            write(17,rec=ind) uu(1:Nx,1:Ny,k,n)
         end do
      end do
      close(17)

      call xy_trans_f( uu, r, Nu, work )

      call z2x_decomp(  r, r, Nu,   uu )

      call z_trans_f( r, Nr, i_symm1, work(1), work(Nze+1) )

      work(1:2)=0.0
      do n=1, 3
         do i=1, nxp
            do j=1, Ny_min
               do k=1, Nz_min
                  diff = r(k,j,i,n) - u(k,j,i,n)
                  if( cabs(diff) .gt. 1.0e-15 ) then
                     work(1) = work(1) + 1.0
                  else
                     work(2) = work(2) + 1.0
                  end if
               end do
            end do
         end do
      end do

      call mpi_allreduce( work(1), work(3), 2, mpi_double_precision,
     &                    mpi_sum, mpi_comm_world, ierr )

      n_bad  = work(3)
      n_good = work(4)

      if(l_root) print *, 'n_good, n_bad = ', n_good, n_bad

      write(ext,'(i6.6)') nt+1
      if(l_root) then
         call write_header( 'header.'//ext, labels, values )
      end if
      call write_field( 'vel.'//ext, r, work )

      call mpi_finalize(ierr)

      stop
      end
