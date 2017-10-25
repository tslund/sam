      subroutine set_labels( lab, val, req, fix )

      include 'sam.h'

      character(12) lab(*)
      real          val(*)
      logical       req(*), fix(*), F

c   *** All parameter labels must be in lower case!
c   *** To add a new parameter, add it both to the list below and to the
c   *** /params/ block in sam.h.  The ordering of parameters here must
c   *** match the ordering in the /params/ common block!

      F= .false.
      req(1:L_params) = .true.
      fix(1:L_params) = .true.
      val(1:L_params) = 0.0

c   *** Static input file parameters.  

      i = 1
c-------------------------------List reals first------------------------
      lab(i) = 'xl'                                              ; i=i+1
      lab(i) = 'yl'                                              ; i=i+1
      lab(i) = 'zl'                                              ; i=i+1
      lab(i) = 'c_smag'       ; fix(i)=F; req(i)=F; val(i)=0.01  ; i=i+1
      lab(i) = 'dt0'          ; fix(i)=F; req(i)=F; val(i)=0.01  ; i=i+1
      lab(i) = 'cfl0'         ; fix(i)=F; req(i)=F; val(i)=1.0   ; i=i+1
      lab(i) = 'vis'          ; fix(i)=F;                        ; i=i+1
      lab(i) = 'k_truncate'   ; fix(i)=F; req(i)=F; val(i)=1.0e+6; i=i+1
      n_inputs_r = i-1
c-------------------------------Integers below this line----------------
      lab(i) = 'nx'                                              ; i=i+1
      lab(i) = 'ny'                                              ; i=i+1
      lab(i) = 'nz'                                              ; i=i+1
      lab(i) = 'i_restart'    ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'nt_restart'   ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'i_strat'      ;           req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'i_les'        ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'n_steps'      ; fix(i)=F;                        ; i=i+1
      lab(i) = 'n_skip_h'     ; fix(i)=F; req(i)=F; val(i)=1e9   ; i=i+1
      lab(i) = 'n_skip_p'     ; fix(i)=F; req(i)=F; val(i)=1e9   ; i=i+1
      lab(i) = 'n_skip_v'     ; fix(i)=F; req(i)=F; val(i)=1e9   ; i=i+1
      lab(i) = 'n_skip_s'     ; fix(i)=F; req(i)=F; val(i)=1e9   ; i=i+1
      lab(i) = 'i_stat'       ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'i_cfl'        ; fix(i)=F; req(i)=F; val(i)=1     ; i=i+1
      lab(i) = 'nrk_max'      ; fix(i)=F; req(i)=F; val(i)=3     ; i=i+1
      lab(i) = 'n_dealias'    ; fix(i)=F; req(i)=F; val(i)=1     ; i=i+1
      lab(i) = 'iubc_z'       ;           req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'i_prob'       ;           req(i)=F; val(i)=2     ; i=i+1
      lab(i) = 'i_force'      ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'k_force'      ; fix(i)=F; req(i)=F; val(i)=2     ; i=i+1
      n_inputs = i-1

c   *** Dynamic parameters for the header file

c-------------------------------List reals first------------------------
      lab(i) = 'time'                                            ; i=i+1
      lab(i) = 'dt'                                              ; i=i+1
      lab(i) = 't_stat'                                          ; i=i+1
      n_dynpar_r = i-n_inputs-1
c-------------------------------Integers below this line----------------
      lab(i) = 'nt'                                              ; i=i+1
      lab(i) = 'n_frame_p'                                       ; i=i+1
      lab(i) = 'n_hist'                                          ; i=i+1
      n_dynpar = i-n_inputs-1

      n_params = n_inputs + n_dynpar

      n_inputs_i = n_inputs - n_inputs_r
      n_dynpar_i = n_dynpar - n_dynpar_r

      return
      end
