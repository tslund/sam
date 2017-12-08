      subroutine assign_params( real_start, int_start, values, 
     &                          n_real, n_int )

      include 'sam.h'

      pointer(ipr, rdat)
      pointer(ipi, idat)
      real    rdat(n_real), values(n_real+n_int)
      integer idat(n_int)

      ipr = loc(real_start)
      ipi = loc(int_start)

      rdat(1:n_real) = values(1:n_real)
      idat(1:n_int ) = values(n_real+1:n_real+n_int)

      return
      end
