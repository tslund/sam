      program tst_grab

      open(unit=1,file='input.dat')
      call grab_input( 1, 'Nx', dummy, Nx, ierr )
      close(1)

      print *, 'Nx = ', Nx

      stop
      end
