      subroutine set_labels( lab, val, req, fix )

      include 'sam.h'

      character(12) lab(L_params)
      real          val(L_params)
      logical       req(L_params), fix(L_params), F

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
      lab(i) = 'k_truncate'   ; fix(i)=F; req(i)=F; val(i)=1e6   ; i=i+1
      lab(i) = 'z0'                     ; req(i)=F; val(i)=0.0   ; i=i+1
      lab(i) = 'to'                     ; req(i)=F; val(i)=1.0   ; i=i+1
      lab(i) = 'lapse0'                 ; req(i)=F; val(i)=0.0   ; i=i+1
      lab(i) = 'grav'         ; fix(i)=F; req(i)=F; val(i)=0.0   ; i=i+1
      lab(i) = 'vis'          ; fix(i)=F                         ; i=i+1
      lab(i) = 'pr'           ; fix(i)=F; req(i)=F; val(i)=1.0   ; i=i+1
      lab(i) = 'amplitude'              ; req(i)=F; val(i)=0.0   ; i=i+1
      lab(i) = 'lambda_x'               ; req(i)=F; val(i)=0.0   ; i=i+1
      lab(i) = 'gam'                    ; req(i)=F; val(i)=0.0   ; i=i+1
      lab(i) = 'omega'                  ; req(i)=F; val(i)=0.0   ; i=i+1
      lab(i) = 'uo'                     ; req(i)=F; val(i)=0.0   ; i=i+1
      lab(i) = 'shear'                  ; req(i)=F; val(i)=1.0   ; i=i+1
      lab(i) = 'shear_ratio'            ; req(i)=F; val(i)=2.5   ; i=i+1
      lab(i) = 'flct_u'       ; fix(i)=F; req(i)=F; val(i)=0.0   ; i=i+1
      lab(i) = 'flct_t'       ; fix(i)=F; req(i)=F; val(i)=0.0   ; i=i+1
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
      lab(i) = 'k_xy_plane1'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'k_xy_plane2'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'k_xy_plane3'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'k_xy_plane4'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'k_xy_plane5'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'k_xy_plane6'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'k_xy_plane7'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'k_xy_plane8'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'k_xy_plane9'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'j_xz_plane1'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'j_xz_plane2'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'j_xz_plane3'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'j_xz_plane4'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'j_xz_plane5'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'j_xz_plane6'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'j_xz_plane7'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'j_xz_plane8'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'j_xz_plane9'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'i_yz_plane1'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'i_yz_plane2'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'i_yz_plane3'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'i_yz_plane4'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'i_yz_plane5'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'i_yz_plane6'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'i_yz_plane7'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'i_yz_plane8'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
      lab(i) = 'i_yz_plane9'  ; fix(i)=F; req(i)=F; val(i)=0     ; i=i+1
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
