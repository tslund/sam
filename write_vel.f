      subroutine write_vel( u )

      include 'sam.h'
      include 'mpif.h'

      real u(Nx,Ny,nzp,Lu)
      integer(kind=mpi_offset_kind) :: offset, irecl_out
      integer amode, status(mpi_status_size)

      n_words_out = Nx*Ny
      irecl_out   = n_words_out*8

      amode = ior( mpi_mode_create, mpi_mode_wronly )

      call mpi_file_open( mpi_comm_world, 'vel.out', amode,
     &                    mpi_info_null, fh(4), ierr )

      do n=1, Lu
         do k=1, nzp
            kk = izs + k-1
            offset = int(((n-1)*Nz + (kk-1)),8)*irecl_out
            call mpi_file_write_at( fh(4), offset, u(1,1,k,n),
     &           n_words_out, mpi_double_precision, status, ierr )
         end do
      end do

      call mpi_file_close( fh(4), ierr )

      return
      end
