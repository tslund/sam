      subroutine open_files( )

      include 'sam.h'

      character(6) ext

c ----- unit   file
c -----  1     stat.nt_restart (closed after read)
c -----  1     stat.nt_final   (closed after write)
c -----  3     umax.out       (remains open)
c -----  4     spectra.out     (remains open)
c -----  8     c.out           (remains open)
c -----  9     cbc.decay       (remains open)

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

      open(unit=3,file='umax.out')
      write(3,3)
3     format(' nt     t          dt       max_div     u_sq',
     &       '       u_max      v_max      w_max' )
      write(6,6)
6     format(' nt     t        energy')

      open(unit=4,file='spectra.out')

      if( i_les  .eq. 1 ) then
         open(unit= 8,file='c.dat')
         write(8,8)
8        format('   nt      t           c')
      end if

      if( i_prob .eq. 2 ) then
         open(unit=9,file='cbc.decay')
         write(9,9)
9        format(' (U0*t/M)+t0  E/E_total_0')
      end if

      end
