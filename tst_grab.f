      program tst_grab

      open(unit=1,file='input.test'
      call grab_input( 1, 'nx', dummy, Nx, ierr )
      close(1)

      print *, 'Nx = ', Nx

      stop
      end
