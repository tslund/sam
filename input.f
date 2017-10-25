      subroutine input( fname, labels, values, fixed )

      include 'sam.h'

      character( *) fname
      character( 8) default_flag
      character(12) labels(*)
      real          values(*)
      logical        fixed(*)
      logical     required(L_params), found(L_params)
      pointer(ipr, rdat)
      pointer(ipi, idat)
      real    rdat(L_params)
      integer idat(L_params)

      call set_labels( labels, values, required, fixed )

      call read_params( fname, labels, values, found, n_inputs, ierr )
      if( ierr .ne. 0 ) stop

      ipr = loc(xL)
      ipi = loc(Nx)

      n_inputs_i = n_inputs - n_inputs_r

      rdat(1:n_inputs_r) = values(1:n_inputs_r)
      idat(1:n_inputs_i) = values(n_inputs_r+1:n_inputs)

      call check_inputs( labels, values, found, required, ierr )
      if( ierr .ne. 0 ) stop

      open(newunit=i_unit,file='params.out')

      do i=1, n_inputs
         default_flag = ' '
         if( .not. found(i) ) default_flag = 'DEFAULT'
         write(i_unit,*) labels(i), values(i), default_flag
      end do

      close(i_unit)

      return
      end
