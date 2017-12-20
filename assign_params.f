      subroutine assign_params( ipr, ipi, values, n_real, n_int )

      pointer(ipr, rdat)
      pointer(ipi, idat)
      real    rdat(n_real), values(n_real+n_int)
      integer idat(n_int)

      rdat(1:n_real) = values(1:n_real)
      idat(1:n_int ) = values(n_real+1:n_real+n_int)

      return
      end
