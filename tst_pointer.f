      program tst_pointer

      parameter( Lp=10 )

      common/inputs/ xL, yL, zL, Nx, Ny, Nz
      common/inputso/ xLo, yLo, zLo, Nxo, Nyo, Nzo

      pointer (ipr, rdat)
      pointer (ipi, idat)
      real    rdat(10)
      integer idat(10)
      real         values(Lp)
      character*12 labels(Lp)
      logical    required(Lp), found(Lp)

      ipr = loc(xL)
      ipi = loc(Nx)
c      ipi = loc(Nx) - (loc(Nx)-loc(xL))/2
      print *, 'ipr, ipi = ', ipr, ipi
      print *, 'loc(xL), loc(Nx) = ', loc(xL), loc(Nx)

      labels(1) = 'xl';  rdat(1) = 10.0;  required(1) = .true.
      labels(2) = 'yl';  rdat(2) = 20.0;  required(2) = .true.
      labels(3) = 'zl';  rdat(3) = 30.0;  required(3) = .true.
      labels(4) = 'nx';  idat(1) = 1000;  required(4) = .true.
      labels(5) = 'ny';  idat(2) = 2000;  required(5) = .true.
      labels(6) = 'nz';  idat(3) = 3000;  required(6) = .true.

      n_inputs_r = 3
      n_inputs_i = 3
      n_inputs   = n_inputs_r + n_inputs_i

c      call assign_vals( labels, rdat, idat, found, n_inputs,
c    &                   n_inputs_r )

      call read_params( 'input.dat', labels, values, found, n_inputs,
     &                   ierr)
      if( ierr .ne. 0 ) stop

c   *** Assign input data to the /params/ common block

      print *, 'about to call assign_params'

      call assign_params( loc(xL), loc(Nx), values,
     &                    n_inputs_r, n_inputs_i )

      j = index_param( 'nz', labels, n_inputs )
      print *, 'position of nx = ', j

      print *, 'xL = ', xL
      print *, 'yL = ', yL
      print *, 'zL = ', zL
      print *, 'Nx = ', Nx
      print *, 'Ny = ', Ny
      print *, 'Nz = ', Nz

      stop
      end

      subroutine assign_vals( labels, values_r, values_i, found,
     &                        n_inputs, n_inputs_r )

      real    values_r(n_inputs_r)
      integer values_i(n_inputs-n_inputs_r)
      character*12 labels(n_inputs), label
      logical found(n_inputs), done

      found = .false.

      i_unit=1
      open(unit=i_unit,file='input.test')

      do L=1, 8
         call read_line( i_unit, label, value, done, ierr )
         if( done ) goto 90
         do i=1, n_inputs
            if( label .eq. labels(i) ) then
               found(i) = .true.
               i1 = i - n_inputs_r
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
