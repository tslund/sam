      subroutine assign_params( loc_real_start, loc_int_start, values, 
     &                          n_real, n_int )

      pointer(ipr, rdat)
      pointer(ipi, idat)
      real    rdat(n_real), values(n_real+n_int)
      integer idat(n_int)

      ipr = loc_real_start
      ipi = loc_int_start

      rdat(1:n_real) = values(1:n_real)
      idat(1:n_int ) = values(n_real+1:n_real+n_int)

      return
      end
