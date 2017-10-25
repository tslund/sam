      subroutine write_field( fname, u, work )

      include 'sam.h'
      include 'mpif.h'

      complex u(Nze,Ny_min,nxp,Lu), work(Nz_min,Ny_min)
      character(*) fname
      character(6) ext
      integer(kind=mpi_offset_kind) :: offset, i_recl_out_8
      integer amode, status(mpi_status_size)

      n_words_out  = Nz_min*Ny_min
      i_recl_out   = n_words_out*16
      i_recl_out_8 = i_recl_out

      amode = ior( mpi_mode_create, mpi_mode_wronly )

      call mpi_file_open( mpi_comm_world, fname, amode,
     &                    mpi_info_null, fh(3), ierr )

      do n=1, Lu
         do i=1, nxp
            ii = ixs + i-1
            do j=1, Ny_min
               do k=1, Nz_min
                  work(k,j) = u(k,j,i,n)
               end do
            end do
            offset = ( (n-1)*Nx_min + (ii-1) )*i_recl_out_8
            call mpi_file_write_at( fh(3), offset, work,
     &           n_words_out, mpi_double_complex, status, ierr )
         end do
      end do

      call mpi_file_close( fh(3), ierr )

      if( l_root .and. i_stat .ne. 0 ) then
         write(ext,'(i6.6)') nt
         open(unit=1,file='stat.'//ext,form='unformatted')
         write(1) Nx, Ny, Nz
         write(1) time, dt, i_cfl, dt0, cfl0, nrk_max
         write(1) vis, i_les
         write(1) t_stat0+(t-t_start)
         write(1) stat
         close(1)
      end if

      return
      end
