      program tst_pointer

      common/inputs/ xL, yL, zL, Nx, Ny, Nz
      common/inputso/ xLo, yLo, zLo, Nxo, Nyo, Nzo

      pointer (ipr, rdat)
      pointer (ipi, idat)
      real    rdat(10)
      integer idat(10)
      character*12 labels(6)
      logical required(6), found(6)

      ipr = loc(xL)
      ipi = loc(Nx)
c      ipi = loc(Nx) - (loc(Nx)-loc(xL))/2
c      print *, 'ipr, ipi = ', ipr, ipi
      

      labels(1) = 'xl';  rdat(1) = 10.0;  required(1) = .true.
      labels(2) = 'yl';  rdat(2) = 20.0;  required(2) = .true.
      labels(3) = 'zl';  rdat(3) = 30.0;  required(3) = .true.
      labels(4) = 'nx';  idat(1) = 1000;  required(4) = .true.
      labels(5) = 'ny';  idat(2) = 2000;  required(5) = .true.
      labels(6) = 'nz';  idat(3) = 3000;  required(6) = .true.

      n_vars_r = 3
      n_vars_i = 3
      n_vars   = n_vars_r + n_vars_i

      call assign_vals( labels, rdat, idat, found, n_vars, n_vars_r )

      j = index_param( 'nz', labels, n_vars )
      print *, 'position of nx = ', j

      ipr = loc(xLo)
      ipi = loc(Nxo)

c      call assign_vals( labels, rdat, idat, found, n_vars, n_vars_r )

      print *, xL
      print *, yL
      print *, zL
      print *, Nx
      print *, Ny
      print *, Nz

      print *, xLo
      print *, yLo
      print *, zLo
      print *, Nxo
      print *, Nyo
      print *, Nzo

      stop
      end

      subroutine assign_vals( labels, values_r, values_i, found,
     &                        n_vars, n_vars_r )

      real    values_r(n_vars_r)
      integer values_i(n_vars-n_vars_r)
      character*12 labels(n_vars), label
      logical found(n_vars), done

      found = .false.

      i_unit=1
      open(unit=i_unit,file='input.test')

      do L=1, 8
         call read_line( i_unit, label, value, ierr, done )
         if( done ) goto 90
         do i=1, n_vars
            if( label .eq. labels(i) ) then
               found(i) = .true.
               i1 = i - n_vars_r
               if( i1 .gt. 0 ) then
                  values_i(i1) = nint(value)
               else
                  values_r(i ) = value
               end if
               goto 20
            end if
         end do
         print *, 'Warning: extraneous variable found ', label
20       continue
      end do

90    close(1)

      return
      end 
