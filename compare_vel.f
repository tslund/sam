      program compare_vel

      include 'sam.h'

      complex, allocatable, dimension(:,:,:,:) :: u1, u2
      complex, allocatable, dimension(:,:)     :: u_in1
      complex, allocatable, dimension(:,:)     :: u_in2
      complex  diff
      character(12) labels(L_params)
      real          values(L_params)
      logical        fixed(L_params), write_params

      write_params = .false.

c------ Set mpi parameters

      myid = 0
      r_root = .true.

c------ Get input parameters

      call input( 'input.dat', labels, values, fixed, write_params )

c------ Initialize constants and other parameter values.  Open files.

      call init( )

      call fft_init( )

c ----- Allocate arrays.

      allocate( u1(Nz_min,Ny_min,Nx_min,3), u2(Nz_min,Ny_min,Nx_min,3) )

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
     &     recl=Nz_min*Ny_min*Nx_min*Lu*16,action='read')
      open(unit=2,file='vel.000001',form='unformatted',access='direct',
     &     recl=Nz_min*Ny_min*Nx_min*Lu*16,action='read')

      open(unit=3,file='good.out')
      open(unit=4,file='bad.out')

      read(1,rec=1) u1
      read(2,rec=1) u2

      close(1)
      close(2)

      print *, 'enter tolerance'
      read(5,*) tol

      n_bad=0;  n_good=0
      do n=1, Lu
         do i=1, Nx_min
            do j=1, Ny_min
               do k=1, Nz_min
                  diff = u1(k,j,i,n) - u2(k,j,i,n)
                  if( cabs(diff) .gt. tol ) then
                     n_bad = n_bad + 1
                     write(4,10), k, j, i, n, cabs(diff)
10                   format(4i5,1p,2e12.4)
                  else
                     n_good = n_good + 1
                     write(3,10), k, j, i, n, u1(k,j,i,n)
                  end if
               end do
            end do
         end do
      end do

      close(3)
      close(4)

      print *, 'n_good, n_bad = ', n_good, n_bad

      stop
      end
