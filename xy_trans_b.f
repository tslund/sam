      subroutine xy_trans_b( u1, u2, n_var, work )

      include 'sam.h'

      complex u1(Ny,Nx/2,nzp,n_var)
      real    u2(Nx,Ny  ,nzp,n_var), work(Nx+1,Ny)

c -------- Transform in y and then x, tranposing the x and y array
c -------- locations in the process.

      do n=1, n_var
         do k=1, nzp

            do i=1, kx_max+1
               ir = 2*(i-1)
               ii = ir+1
c      if(k.eq.4.and.i.eq.5) then
c         do j=1, Ny
c            write(6,'(i5,1p,2e12.4)') j, u1(j,i,k,n)
c         end do
c      end if
               call cfftb( Ny, u1(1,i,k,n), trigy )
c      if(k.eq.4.and.i.eq.5) then
c         do j=1, Ny
c            write(6,'(i5,1p,2e12.4)') j, u1(j,i,k,n)
c         end do
c      end if
               if( i .eq. 1 ) then
                  do j=1, Ny
                     work(i,j) = real( u1(j,i,k,n))
                  end do
               else
                  do j=1, Ny
                     work(ir,j) = real( u1(j,i,k,n))
                     work(ii,j) = aimag(u1(j,i,k,n))
                  end do
               end if
            end do
            do j=1, Ny
               do i=2*kx_max+2, Nx
                  work(i,j) = 0.0
               end do
            end do

            do j=1, Ny
c      if(k.eq.4.and.j.eq.3) then
c         do i=1, Nx
c            write(6,'(i5,1p,2e12.4)') i, work(i,j)
c         end do
c      end if
               call rfftb( Nx, work(1,j), trigx )
c      if(k.eq.4.and.j.eq.3) then
c         do i=1, Nx
c            write(6,'(i5,1p,2e12.4)') i, work(i,j)
c         end do
c      end if
               do i=1, Nx
                  u2(i,j,k,n) = work(i,j)
               end do
            end do

         end do
      end do

      return
      end
