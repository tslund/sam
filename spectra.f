      subroutine spectra( u, r, sample, ek, tk, dk, vt, energy_w,
     &                    div_rms, div_max, n_div, nt_write, work )

      include 'sam.h'
      include 'mpif.h'

      parameter( nw=6 )
      real k_mag
      complex u(Nze,Ny_min,nxp,Lu), r(Nze,Ny_min,nxp,Lr), div
      real    sample(Lm), ek(Lm), tk(Lm), dk(Lm), vt(Lm),
     &        work(Lm,2*nw)

      if( ixs .eq. 1 ) then
         do n=1, Lu
            u_avg(n) = u(1,1,1,n)
         end do
      end if

      work(1:Lm,1:nw) = 0.0

      if( i_les .eq. 1 ) fac_cl = sqrt( 2*E_cl/float(k_cl) )
      do i=1, nxp
         ii = ixs + i-1
         factor = 2.0
         if( ii .eq. 1 .or. ii .eq. Nx2+1 ) factor = 1.0
         kx2 = k_x(ii)**2
         rkx = wave_x*float(k_x(ii))
         rkx2 = rkx**2
         do j=1, Ny_min
            k_sq2 = kx2 + k_y(j)**2
            rky = wave_y*float(k_y(j))
            rk_sq2 = rkx2 + rky**2
            do k=1, Nz_min
               k_sq  = k_sq2 + k_z(k)**2
               k_mag = sqrt(float(k_sq))
               is = 1 + int( k_mag + 0.5 )
               rkz = wave_z*float(k_z(k))
               rk_sq = rk_sq2 + rkz**2
               if( k_sq .eq. 0 .or. i_les .eq. 0 ) then
                  vist = 0.0
               else
                  vist = fac_cl*( a_cl + b_cl*exp(-(c_cl*k_cl)/k_mag) )
               end if
               u_sq = u(k,j,i,1)*conjg(u(k,j,i,1)) +
     &                u(k,j,i,2)*conjg(u(k,j,i,2)) +
     &                u(k,j,i,3)*conjg(u(k,j,i,3))
               udu  = r(k,j,i,1)*conjg(u(k,j,i,1)) +
     &                r(k,j,i,2)*conjg(u(k,j,i,2)) +
     &                r(k,j,i,3)*conjg(u(k,j,i,3))
               work(is,1) = work(is,1) + factor*u_sq
               work(is,2) = work(is,2) + factor*udu
               work(is,3) = work(is,3) + factor*rk_sq*(vis+vist)*u_sq
               work(is,4) = work(is,4) + factor*vist
               work(is,5) = work(is,5) + factor
               div = rkx*u(k,j,i,1) +
     &               rky*u(k,j,i,2) +
     &               rkz*u(k,j,i,3)
               div_sq = div*conjg(div)
               work(1,6) = work(1,6) +      div_sq
               work(2,6) = max( work(2,6),  div_sq )
c               if( div_sq .gt. 1.e-24 ) then
c                  print *, 'divergence ',k,j,i,div
c                  work(3,6) = work(3,6) + 1
c               end if
            end do
         end do
      end do

      call mpi_allreduce( work(1,1), work(1,nw+1), Lm*nw,
     &     mpi_double_precision, mpi_sum, mpi_comm_world, ierr )

      sum_ek = 0.0
      sum_tk = 0.0
      sum_dk = 0.0

      do i=1, Lm
         ek(    i) = work(i,nw+1)
         tk(    i) = work(i,nw+2)
         dk(    i) = work(i,nw+3)
         vt(    i) = work(i,nw+4)
         sample(i) = work(i,nw+5)
         sum_ek = sum_ek + ek(i)
         sum_tk = sum_tk + tk(i)
         sum_dk = sum_dk + dk(i)
      end do
      energy_w = 0.5*sum_ek

      div_rms = sqrt(work(1,nw+6))
      div_max = sqrt(work(2,nw+6))
      n_div   = work(3,nw+6)

      if(l_root) then

         write(4,10) nt_write, time, 0.5*sum_ek, sum_dk
10       format('nt =', i3, '  t =', 1pe11.4, '  KE =', 1pe11.4,
     &          '  Diss =', 1pe11.4, //,
     &          '           k     E(k)')

         do i=1, Lm
            if( sample(i) .ne. 0.0 ) then
               write(4,20) i-1, 4.0*pi*((i-1)**2)/sample(i)*(0.5*ek(i)),
     &                     tk(i)/sample(i), dk(i)/sample(i), 
     &                     vt(i)/sample(i)
c        &                  tk(i)
20             format( i12, 1p, 4e12.4 )
            end if
         end do

         write(4,*)

      end if

      return
      end
