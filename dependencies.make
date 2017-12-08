sam: assign_params.o assign_run_time.o assign_static.o check_inputs.o cl_energy.o div_check.o fft_init.o fft.o fileop_subs.o force.o get_mean.o initial_field.o init.o input.o open_files.o read_field.o read_header.o read_line.o read_params.o rhs.o sam.o set_labels.o set_range.o spectra.o stop_on_error.o vel_max.o write_field.o write_header.o write_planes.o x2z_decomp.o x_deriv.o xy_trans_b.o xy_trans_f.o y_deriv.o z2x_decomp.o z_trans_b.o z_trans_f.o

tst_input: assign_params.o assign_run_time.o assign_static.o check_inputs.o input.o read_header.o read_line.o read_params.o set_labels.o set_range.o stop_on_error.o tst_input.o write_header.o

tst_transpose: assign_params.o assign_run_time.o assign_static.o check_inputs.o div_check.o fft_init.o fft.o fileop_subs.o initial_field.o init.o input.o open_files.o read_field.o read_header.o read_line.o read_params.o set_labels.o set_range.o spectra.o stop_on_error.o tst_transpose.o write_field.o write_header.o x2z_decomp.o xy_trans_b.o xy_trans_f.o z2x_decomp.o z_trans_b.o z_trans_f.o

check: assign_params.o check_inputs.o check.o fft_init.o fft.o init.o input.o read_line.o read_params.o set_labels.o set_range.o

compare_vel: assign_params.o check_inputs.o compare_vel.o fft_init.o fft.o init.o input.o read_line.o read_params.o set_labels.o

planes2vtk: assign_params.o check_inputs.o deriv1.o deriv2.o fft_init.o fft.o init.o input.o planes2vtk.o read_line.o read_params.o set_labels.o

remesh: assign_params.o check_inputs.o fileop_subs.o grab_input.o init.o input.o read_header.o read_line.o read_params.o remesh.o set_labels.o to_lower.o write_header.o

tst_pointer: read_line.o tst_pointer.o

assign_params.o assign_run_time.o assign_static.o check_inputs.o check.o cl_energy.o compare_vel.o div_check.o fft_init.o force.o get_mean.o initial_field.o init.o input.o open_files.o planes2vtk.o read_field.o read_header.o remesh.o rhs.o sam.o set_labels.o set_range.o spectra.o stop_on_error.o tst_input.o tst_transpose.o vel_max.o write_field.o write_header.o write_planes.o x2z_decomp.o x_deriv.o xy_trans_b.o xy_trans_f.o y_deriv.o z2x_decomp.o z_trans_b.o z_trans_f.o: sam.h

