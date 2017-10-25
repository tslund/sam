      subroutine write_header( fname, labels, values )

      include 'sam.h'

      character( *) fname
      character(12) labels(*)
      real          values(*)
      pointer(ipr, rdat)
      pointer(ipi, idat)
      real    rdat(L_params)
      integer idat(L_params)

      ipr = loc(xL)
      ipi = loc(Nx)

      i_beg=1;        i_dat=n_inputs_r;  i_end=i_beg+i_dat-1
         values(i_beg:i_end) = rdat(1:i_dat)
      i_beg=i_end+1;  i_dat=n_inputs_i;  i_end=i_beg+i_dat-1
         values(i_beg:i_end) = idat(1:i_dat)

      ipr = loc(time)
      ipi = loc(nt)

      i_beg=i_end+1;  i_dat=n_dynpar_r;  i_end=i_beg+i_dat-1
         values(i_beg:i_end) = rdat(1:i_dat)
      i_beg=i_end+1;  i_dat=n_dynpar_i;  i_end=i_beg+i_dat-1
         values(i_beg:i_end) = idat(1:i_dat)

      open(newunit=i_unit,file=fname)

      do i=1, n_params
         write(i_unit,*) labels(i), values(i)
      end do

      close(i_unit)

      return
      end
