      subroutine open_files( )

      include 'sam.h'
      include 'mpif.h'

      integer amode
      character( 1) ext1
      character( 6) ext
      character(10) file_access

c ----- unit   file
c -----  1     stat.nt_restart (closed after read)
c -----  1     stat.nt_final   (closed after write)
c -----  2     energy.out      (remains open)
c -----  3     umax.out        (remains open)
c -----  4     spectra.out     (remains open)
c -----  7     tstep.out       (remains open)
c -----  8     c.out           (remains open)
c -----  9     cbc.decay       (remains open)
c ----- 10     time.out        (remains open)
c ----- 12     variance.out    (remains open)
c ----- 20+m   xy[m].out       (remains open)

c ----- Open files for the xy planes data.

      do k=1, nzp
         kk = izs + k-1
         do m=1, n_xy_planes
            write(ext1,'(i1.1)') m
            if( kk .eq. k_xy_plane(m) ) then
              open(unit=20+m,file='xy'//ext1//'.out',form='unformatted',
     &             access='direct',recl=Nx*Ny*8,action='write')
            end if
         end do
      end do

c ----- Open files for the xz and yz planes data.

      amode = ior( mpi_mode_create, mpi_mode_wronly )

      do m=1, n_xz_planes
         write(ext1,'(i1.1)') m
         call mpi_file_open( mpi_comm_world, 'xz'//ext1//'.out', amode,
     &                       mpi_info_null, fh(30+m), ierr )
      end do

      do m=1, n_yz_planes
         write(ext1,'(i1.1)') m
         call mpi_file_open( mpi_comm_world, 'yz'//ext1//'.out', amode,
     &                       mpi_info_null, fh(40+m), ierr )
      end do

c ----- Open a file for the horizontal means

      if(l_root) then
         open(unit=11,file='mean.out',form='unformatted',
     &        access='direct',recl=Nze*Lu*8,action='write')
      end if

c ----- Only roots initializes statistics and opens ascii output files.

      if( .not. l_root ) return

c------ Set file acces in order to append to files from a restart run.

      if( i_restart .eq. 1 ) then
         file_access = 'append'
      else
         file_access = 'sequential'
      end if

c------ Initialize statistics arrays.

      if( i_stat .gt. 1 ) then
         write(ext,'(i6.6)') nt_restart
         open(unit=1,file='stat.'//ext,form='unformatted')
         do n=1, 3
            read(1)
         end do
         read(1) t_stat0
         read(1) stat
         close(1)
      else
         t_stat0 = 0.0
         do n=1, Ls
            stat(n) = 0.0
         end do
      endif

c ----- Open ascii output files.

      open(unit=2,file='energy.out',access=file_access)

      if( i_restart .eq. 0 ) write(2,2)
2     format('#   nt    time       energy')

      open(unit=3,file='umax.out',access=file_access)

      if( i_restart .eq. 0 ) write(3,3)
3     format('#   nt    time       div_max      u_max       v_max',
     &       '       w_max       T_max' )

      open(unit=4,file='spectra.out',access=file_access)

      open(unit=7,file='tstep.out',access=file_access)

      if( i_restart .eq. 0 ) write(7,7)
7     format('#   nt    time         dt         cfl_x       cfl_y',
     &       '       cfl_z      cfl_vis')

      if( i_les  .eq. 1 ) then
         open(unit= 8,file='c.dat',access=file_access)
         if( i_restart .eq. 0 ) write(8,8)
8        format('   nt      t           c')
      end if

      if( i_prob .eq. 2 ) then
         open(unit=9,file='cbc.decay',access=file_access)
         if( i_restart .eq. 0 ) write(9,9)
9        format(' (U0*t/M)+t0  E/E_total_0')
      end if

      open(unit=10,file='time.out',access=file_access)

      if( i_restart .eq. 0 ) write(10,10)
10    format('#   nt    time  ')

      open(unit=12,file='variance.out',access=file_access)

      if( i_restart .eq. 0 ) write(12,12)
12    format('#   nt    time        u_var       v_var       w_var',
     &       '      0.5u_sq')

c ----- Write the maximum wavenumbers and individual processor 
c ----- index ranges

      open(newunit=i_unit,file='dims.out')

      write(i_unit,20) kx_max, ky_max, kz_max,
     &                 Nx_min, Ny_min, Nz_min, Nze
20    format('   kx_max ky_max kz_max'    ,/,3i7,//,
     &       '   Nx_min Ny_min Nz_min  Nz',/,4i7,//,
     &       '    myid   ixs    ixe    izs    ize    nxp    nzp')

      do m=0, numprocs-1
         write(i_unit,25) m, ix_s(m), ix_e(m), iz_s(m), iz_e(m),
     &                    nx_p(m), nz_p(m)
25       format(7i7)
      end do

      close(i_unit)

      return
      end
