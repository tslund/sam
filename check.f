      program check

      include 'sam.h'

      complex, allocatable, dimension(:,:,:,:) :: u1, u2, u3
      complex, allocatable, dimension(:,:)     :: u_in1
      complex, allocatable, dimension(:)       :: u_in2
      complex  diff
      character(12) labels(L_params)
      real          values(L_params)
      logical        fixed(L_params)
      integer       i_symm1(4)

c------ Set mpi parameters

      print *, 'enter the number of processors used for the sam run'
      read(5,*) numprocs

      myid = 0
      r_root = .true.

c------ Get input parameters

      call input( 'input.dat', labels, values, fixed )

c------ Initialize constants and other parameter values.  Open files.

      call init( )

      call fft_init( )

      i_symm1 = i_symm

c ----- Determine index ranges for each process.

      call set_range( Nze, Nx_min )

c ----- Allocate arrays.

      allocate( u1(Nz,Ny_min,Nx_min,3), u2(Nz,Ny_min,Nx_min,3),
     &          u3(Ny,Nx/2,Nz,3), u_in1(Nz_min,Ny_min), u_in2(Ny_min) )

      js3 = Ny - Ny_min
      je1 = ky_max+1
      jb2 = je1+1
      je2 = je1 + js3
      jb3 = je2+1

      ks3 = Nz - Nz_min
      ke1 = kz_max+1
      kb2 = ke1+1
      ke2 = ke1 + ks3
      kb3 = ke2+1

      open(unit=1,file='vel.000000',form='unformatted',access='direct',
     &     recl=Nz_min*Ny_min*16,action='read')

      L = 0
      do n=1, 3
         do i=1, Nx_min
            L = L + 1
            read(1,rec=L) u_in1
            do j=1, Ny_min
               do k=1, ke1
                  u1(k,j,i,n) = u_in1(k,j)
               end do
               do k=kb2, ke2
                  u1(k,j,i,n) = cmplx(0.0,0.0)
               end do
               do k=kb3, Nz
                  u1(k,j,i,n) = u_in1(k-ks3,j)
               end do
               call cfftb( Nz, u1(1,j,i,n), trigz )
            end do
         end do
      end do

      close(1)

c      open(unit=2,file='spec_vel1.out',form='unformatted',
      open(unit=2,file='post_vel1.out',form='unformatted',
     &     access='direct',recl=Nz*Ny_min*Nx_min*3*16,action='read')
      read(2,rec=1) u2
      close(2)

      n_good = 0; n_bad = 0
      do n=1, 3
         do i=1, Nx_min
            do j=1, Ny_min
               do k=1, Nz
                  diff = u2(k,j,i,n) - u1(k,j,i,n)
                  if( cabs(diff) .lt. 1.0e-8 ) then
                     n_good = n_good + 1
c                     print *, k, j, i, n
                  else
                     n_bad = n_bad + 1 
c                     print *, k, j, i, n, cabs(diff)
c                     print *, u1(k,j,i,n)
c                     print *, u2(k,j,i,n)
c                     print *, diff
                  end if
c                  if( j .eq. 1 .and. k .eq. 1 )
c     & print *, u1(i,j,k,n), u2(i,j,k,n)
               end do
            end do
         end do
      end do

      print *, 'n_good, bad spec_vel1 = ', n_good, n_bad

c      i=1; j=1; n=1
c      do k=1, Nz
c         write(6,'(i5,1p,4e12.4)') k, u2(k,j,i,n), u1(k,j,i,n)
c      end do

      open(unit=2,file='vel1.out',form='unformatted',access='direct',
     &     recl=Ny_min*16,action='read')

      Ny_minNx_min = Ny_min*Nx_min
      Ny_minNx_minNze = Ny_min*Nx_min*Nze

      L = 0
      do mx=0, numprocs-1
         do mz=0, numprocs-1
            do n=1, Lu
               ind1 = Ny_minNx_minNze*(n-1)
               do k=1, nz_p(mz)
                  kk = iz_s(mz) + k-1
                  ind2 = ind1 + Ny_minNx_min*(kk-1)
                  do i=1, nx_p(mx)
                     ii = ix_s(mx) + i-1
                     ind3 = ind2 + Ny_min*(ii-1)
                     L = L + 1
                     read(2,rec=L) u_in2
                     do j=1, Ny_min
                        ind = ind3 + j
                        u2(kk,j,ii,n) = u_in2(j)
c      write(70+myid,70) j, ii, kk, n, mz, mx, L, ind,
c     &                  real(u_in2(j))
c70    format(8i5,f10.3)
                     end do
                  end do
               end do
            end do
         end do
      end do

      close(2)

      n_good = 0
      do n=1, 3
         do i=1, Nx_min
            do j=1, Ny_min
               do k=1, Nz
                  diff = u2(k,j,i,n) - u1(k,j,i,n)
                  if( cabs(diff) .lt. 1.0e-10 ) then
                     n_good = n_good + 1
                  else
c                     print *, k, j, i, n, cabs(diff)
c                     print *, u1(k,j,i,n)
c                     print *, u2(k,j,i,n)
c                     print *, diff
                  end if
c                  if( j .eq. 1 .and. k .eq. 1 )
c     & print *, u1(i,j,k,n), u2(i,j,k,n)
               end do
            end do
         end do
      end do

      print *, 'n_good1 = ', n_good

      open(unit=3,file='vel2.out',form='unformatted',access='direct',
     &     recl=Ny*Nx/2*Nz*3*16,action='read')
      read(3,rec=1) u3
      close(3)

      n_good = 0
      L = 0
      do n=1, 3
         do k=1, Nz
            do i=1, Nx_min
               do j=1, Ny_min
                  jj = j
                  if( j .ge. jb2 ) jj = j + js3
                  L = L + 1
                  diff = u2(k,j,i,n) - u3(jj,i,k,n)
                  if( cabs(diff) .lt. 1.0e-10 ) then
                     n_good = n_good + 1
c                     print *, 'good', k, j, i, n
c                     print *, u2(k,j,i,n)
                  else
                     write(6,'(5i6,2f10.3)') k, j, i, n, L,
     &                        real(real(u2(k ,j,i,n)),kind=4),
     &                        real(real(u3(jj,i,k,n)),kind=4)
                  end if
               end do
            end do
         end do
      end do

      print *, 'n_good2 = ', n_good

      stop
      end
