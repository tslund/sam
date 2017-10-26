      subroutine open_files( )

      include 'sam.h'

      character( 6) ext
      character(10) file_access

c ----- unit   file
c -----  1     stat.nt_restart (closed after read)
c -----  1     stat.nt_final   (closed after write)
c -----  2     energy.out     (remains open)
c -----  3     umax.out       (remains open)
c -----  4     spectra.out     (remains open)
c -----  7     tstep.out       (remains open)
c -----  8     c.out           (remains open)
c -----  9     cbc.decay       (remains open)

c------ Set file acces in order to append to files from a restart run

      if( i_restart .eq. 1 ) then
         file_access = 'append'
      else
         file_access = 'sequential'
      end if

c------ Initialize statistics arrays

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

      open(unit=12,file='variance.out',access=file_access)

      if( i_restart .eq. 0 ) write(12,12)
12    format('#   nt    time        u_var       v_var       w_var',
     &       '      0.5u_sq')

      end
