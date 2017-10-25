      parameter( mLx=2050, mLy=2050, mLz=2050 )
      parameter( Ls=54, L_params=100 )
      parameter( maxp=1024)

      real Nz_inv, NxNy_inv, NxNyNz_inv
      real l_cbc, m_cbc
      real k_truncate, k_truncate_sq
      complex iunit
      integer fh

      common/params/ xL, yL, zL, c_smag, dt0, cfl0, vis, k_truncate,
     &               Nx, Ny, Nz, i_restart, nt_restart, i_strat, i_les,
     &               n_steps, n_skip_h, n_skip_p, n_skip_v, n_skip_s,
     &               i_stat, i_cfl, nrk_max, n_dealias, iubc_z, 
     &               i_prob, i_force, k_force

      common/run_time/ time, dt, t_stat, nt, n_frame_p, n_hist

      common/consts/ pi, two_pi, iunit

      common/dimensions/ Lu, Lr, Lw, Lm, Nze, ipad_r,
     &                   Nx_min, Ny_min, Nz_min,
     &                   n_inputs_r, n_inputs_i, n_inputs,
     &                   n_dynpar_r, n_dynpar_i, n_dynpar,
     &                   n_params

      common/mesh/ dx, dy, dz, dx_inv, dy_inv, dz_inv

      common /mpi_stuff/ myid, numprocs, l_root, i_root,
     &                   ixs, ixe, izs, ize,
     &                   nxp, nzp, fh(100),
     &                   ix_s( 0:maxp), ix_e( 0:maxp), nx_p(0:maxp),
     &                   iz_s( 0:maxp), iz_e( 0:maxp), nz_p(0:maxp)

      common/tstep/ gamma(4), zeta(4), beta(4),
     &              vis_cfl, w_cfl, cfl_by_dt, t_start,
     &              vis_fac_x(mLx), vis_fac_y(mLy), vis_fac_z(mLz),
     &              nrk, nt_start

      common/les/ E_cl, a_cl, b_cl, c_cl, k_cl

      common/cbc/ l_cbc, m_cbc, t0_cbc, u_inf_cbc, u_0_cbc, tfact

      common/de_alias/ z_shift, k_truncate_sq

      common/fft/ trigx( 2*mLx+15), trigy( 4*mLy+15), trigz(4*mLz+15),
     &            trigzs(3*mLz+15), trigzc(3*mLz+15),
     &            wave_x, wave_y, wave_z,
     &            k_x(mLx), k_y(mLy), k_z(mLz),
     &            i_symm(4), Nz_inv, NxNy_inv, NxNyNz_inv, 
     &            kx_max, ky_max, kz_max,
     &            Nx2, NxNy, NxNyNz

      common/statistics/ t_stat0, 
     &                   u_avg(4), u_var(4), u_max(4), vis_max, energy,
     &                   cfl_x, cfl_y, cfl_z, cfl_vis,
     &                   stat(Ls), stat_tmp(Ls)