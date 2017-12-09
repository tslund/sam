      subroutine assign_run_time( values_h )

      include 'sam.h'

      real values_h(L_params)

      dt_input = dt

      call assign_params( loc(time), loc(nt), values_h(n_inputs+1), 
     &                    n_dynpar_r, n_dynpar_i )

      t_start  = time
      nt_start = nt
      dt = dt_input
      n_frame_p = 0

      return
      end
