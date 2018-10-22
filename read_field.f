      subroutine read_field( fname, u, work )

      include 'sam.h'
      include 'mpif.h'

      complex u(Nze,Ny_min,nxp,Lu), work(Nz_min,Ny_min)
      character(*) fname
      character(6) ext
      integer(kind=mpi_offset_kind) :: offset, irecl_in
      integer status(mpi_status_size)
      integer(8) isize_found, isize_expect, n_rec

      n_words_in = Nz_min*Ny_min
      irecl_in   = n_words_in*16

      isize_expect = irecl_in*int(Nx_min*Lu,8)

      if(l_root) then
         isize_found  = isize_expect
         call get_fsize( fname, Lu*8, isize_found, n_rec )
      end if

      call mpi_bcast( isize_found, 1, mpi_integer8,
     &                i_root, mpi_comm_world, ierr )

      if( isize_found .ne. isize_expect ) then
         if(l_root) print *, 'ERROR: Vel file is not of correct size'
         if(l_root) print *, 'isize_expect, isize_found = ',
     &                        isize_expect, isize_found
         call mpi_finalize( ierr )
         stop
      end if

      call mpi_file_open( mpi_comm_world, fname, mpi_mode_rdonly,
     &                    mpi_info_null, fh(2), ierr )

      do n=1, Lu
         do i=1, nxp
            ii = ixs + i-1
            offset = int(((n-1)*Nx_min + (ii-1)),8)*irecl_in
            call mpi_file_read_at( fh(2), offset, work,
     &           n_words_in, mpi_double_complex, status, ierr )
            do j=1, Ny_min
               do k=1, Nz_min
                  u(k,j,i,n) = work(k,j)
               end do
               do k=Nz_min+1, Nze
                  u(k,j,i,n) = cmplx(0.0,0.0)
               end do
            end do
         end do
      end do

      call mpi_file_close( fh(2), ierr )

      return
      end
