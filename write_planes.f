      subroutine write_planes( uu, du, u_xz, u_yz, mean, work )

c ----- The fields contained in uu array are ordered as u, v, w, [T],
c ----- p, dudz, dvdz, [dTdz], where the fields involving T are not
c ----- present for unstratified cases.

      include 'sam.h'
      include 'mpif.h'

      real   uu(Nx,Ny,nzp,Lr), mean(Nze,Lu), work(Lw)
      real   du(Nx,Ny), u_xz(Nx,nzp,2*Lu), u_yz(Ny,nzp,2*Lu)
      real   jump
      integer(kind=mpi_offset_kind) :: offset, offset_0, offset1,
     &                                 NxNze, xz_block_size,
     &                                 NyNze, yz_block_size
      integer status(mpi_status_size)

      n_frame_p = n_frame_p + 1

      n_out = 2*Lu

c ----- Write the time.

      if(l_root) then
         write(10,10) nt, time
10       format(i6,1p,e12.4)
         flush(10)
      end if

c ----- Write the means.

      if(l_root) then
         write(11,rec=n_frame_p) mean
      end if

c ----- Write the xy planes.

      do k=1, nzp
         kk = izs + k-1
         do m=1, n_xy_planes
            if( kk .eq. k_xy_plane(m) ) then
               do n=1, n_out
                  irec = (n_frame_p-1)*n_out + n
                  write(20+m,rec=irec) uu(1:Nx,1:Ny,k,n)
               end do
            end if
         end do
      end do

c ----- Gather and write the xz planes data.

      n_words_out = Nx*nzp
      NxNze = int(Nx,8)*int(Nze,8)
      xz_block_size = NxNze*int(n_out*8,8)
      offset_0 = int((n_frame_p-1),8)*xz_block_size
      offset_1 = offset_0 + int((izs-1),8)*int(Nx*8,8)

      do m=1, n_xz_planes
         j = j_xz_plane(m)

         do n=1, Lu+1
            do k=1, nzp
c         write(70+myid,'(i5,1p,4e12.4)') k, uu(40,j,k,n)
               do i=1, Nx
                  u_xz(i,k,n) = uu(i,j,k,n)
               end do
            end do
         end do

         do k=1, nzp
            L=1
            do 20 n=1, Lu

               if( n .eq. 2 ) goto 20

               L = L + 1

               call y_deriv( uu(1,1,k,n), du, work )

               do i=1, Nx
                  u_xz(i,k,Lu+L) = du(i,j)
               end do

20          continue
         end do

         do n=1, n_out
            offset = offset_1 + (n-1)*NxNze*8
c      print *, 'myid, offset = ', myid, offset
      i=10
c      do k=1, nzp
c         write(70+myid,'(i5,1p,4e12.4)') k, u_xz(i,k,n)
c      end do
            call mpi_file_write_at( fh(30+m), offset, u_xz(1,1,n),
     &              n_words_out, mpi_double_precision, status, ierr )
c            call mpi_finalize(ierr)
c            stop
         end do

      end do

c ----- Gather and write the yz planes data.

      n_words_out = Ny*nzp
      NyNze = int(Ny,8)*int(Nze,8)
      yz_block_size = NyNze*int(n_out*8,8)
      offset_0 = int((n_frame_p-1),8)*yz_block_size
      offset_1 = offset_0 + int((izs-1),8)*int(Ny*8,8)

      do m=1, n_yz_planes
         i = i_yz_plane(m)

         do n=1, Lu+1
            do k=1, nzp
               do j=1, Ny
                  u_yz(j,k,n) = uu(i,j,k,n)
               end do
            end do
         end do

         do k=1, nzp
            L=1
            do 30 n=1, Lu

               if( n .eq. 1 ) goto 30

               L = L + 1

               call x_deriv( uu(1,1,k,n), du, work )

               do j=1, Ny
                  u_yz(j,k,Lu+L) = du(i,j)
               end do

30          continue
         end do

         do n=1, n_out
            offset = offset_1 + (n-1)*NyNze*8
            call mpi_file_write_at( fh(40+m), offset, u_yz(1,1,n),
     &              n_words_out, mpi_double_precision, status, ierr )
         end do

      end do

      return
      end
