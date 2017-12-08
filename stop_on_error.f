      subroutine stop_on_error( ierr, i_bcast )

      include 'sam.h'
      include 'mpif.h'

      if( i_bcast .eq. 1 ) then
         call mpi_bcast( ierr, 1, mpi_integer, i_root,
     &                mpi_comm_world, ierr_mpi )
      end if

      if( ierr .ne. 0 ) then
         call mpi_finalize( ierr_mpi )
         stop
      end if

      return
      end
