      subroutine get_mean( u, i_symm1, i_remove, u_mean, jump, work )

      include 'sam.h'

      complex u(Nze,Ny_min,nxp)
      real u_mean(Nze), work(Nze), jump, jump1

      select case( i_symm1 )
      case(0)
         u_mean(1) = real(u(1,1,1))
         do k=2, kz_max+1
            kr = 2*(k-1)
            ki = kr+1
            u_mean(kr) = real( u(k,1,1))
            u_mean(ki) = aimag(u(k,1,1))
         end do
         do k=2*kz_max+2, Nz
            u_mean(k) = 0.0
         end do
         call rfftb( Nz, u_mean, trigzr )
         u_n = 3.0*u_mean( 1) - 3.0*u_mean(   2) + u_mean(   3)
         u_1 = 3.0*u_mean(Nz) - 3.0*u_mean(Nz-1) + u_mean(Nz-2)
         jump = 0.5*( (u_mean(Nz)-u_n) + (u_1-u_mean(1)) )
         if( i_remove .eq. 1 ) then
            do k=1, Nz
               u_lin = jump*( float(k-1)*Nz_inv - 0.5 )
               work(k) = u_lin
            end do
            call rfftf( Nz, work, trigzr )
            do k=1, kz_max+1
               work(k) = work(k)*Nz_inv
            end do
            u(1,1,1) = u(1,1,1) - cmplx(work(1),0.0)
            do k=2, kz_max+1
               k2 = Nz_min+2 - k 
               kr = 2*(k-1)
               ki = kr+1
               u(k ,1,1) = u(k ,1,1) - cmplx(work(kr), work(ki))
               u(k2,1,1) = u(k2,1,1) - cmplx(work(kr),-work(ki))
            end do
            u_mean(1) = real(u(1,1,1))
            do k=2, kz_max+1
               kr = 2*(k-1)
               ki = kr+1
               u_mean(kr) = real( u(k,1,1))
               u_mean(ki) = aimag(u(k,1,1))
            end do
            do k=2*kz_max+2, Nz
               u_mean(k) = 0.0
            end do
            call rfftb( Nz, u_mean, trigzr )
            u_n = 3.0*u_mean( 1) - 3.0*u_mean(   2) + u_mean(   3)
            u_1 = 3.0*u_mean(Nz) - 3.0*u_mean(Nz-1) + u_mean(Nz-2)
            jump1 = 0.5*( (u_mean(Nz)-u_n) + (u_1-u_mean(1)) )
            jump = jump - jump1
         end if
      case(-1)
         do k=2, kz_max+1
            u_mean(k) = real(u(k,1,1))
         end do
         do k=kz_max+2, Nze
            u_mean(k) = 0.0
         end do  
         call sint( Nz-1, u_mean(2), trigzs )
         u_mean(   1) = 0.0
         u_mean(Nz+1) = 0.0
         jump = 0.0
      case( 1 )
         do k=1, kz_max+1
            u_mean(k) = real(u(k,1,1))
         end do
         do k=kz_max+2, Nze
            u_mean(k) = 0.0
         end do  
         call cost( Nze, u_mean, trigzc )
c         jump = u_mean(Nze) - u_mean(1)
         jump = 0.0
      end select

      return
      end
