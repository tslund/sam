      subroutine write_vel( u )

      include 'sam.h'
      include 'mpif.h'

      real u(Nx,Ny,nzp,Lu)
      integer(kind=mpi_offset_kind) :: offset, i_recl_out_8
      integer amode, status(mpi_status_size)

      n_words_out  = Nx*Ny
      i_recl_out   = n_words_out*8
      i_recl_out_8 = i_recl_out

      amode = ior( mpi_mode_create, mpi_mode_wronly )

      call mpi_file_open( mpi_comm_world, 'vel.out', amode,
     &                    mpi_info_null, fh(4), ierr )

      do n=1, Lu
         do k=1, nzp
            kk = izs + k-1
            offset = ( (n-1)*Nz + (kk-1) )*i_recl_out_8
            call mpi_file_write_at( fh(4), offset, u(1,1,k,n),
     &           n_words_out, mpi_double_precision, status, ierr )
         end do
      end do

      call mpi_file_close( fh(4), ierr )

      return
      end
