sam: check_inputs.o cl_energy.o div_check.o fft_init.o fft.o fileop_subs.o force.o initial_field.o init.o input_p.o open_files.o read_field.o read_header.o read_line.o read_params.o rhs.o sam.o set_labels.o set_range.o spectra.o vel_max.o write_field.o write_header.o write_vel.o x2z_decomp.o xy_trans_b.o xy_trans_f.o z2x_decomp.o z_trans_b.o z_trans_f.o

tst_input: check_inputs.o input_p.o read_header.o read_line.o read_params.o set_labels.o set_range.o tst_input.o write_header.o

tst_transpose: check_inputs.o div_check.o fft_init.o fft.o fileop_subs.o initial_field.o init.o input_p.o open_files.o read_field.o read_header.o read_line.o read_params.o set_labels.o set_range.o spectra.o tst_transpose.o write_field.o write_header.o x2z_decomp.o xy_trans_b.o xy_trans_f.o z2x_decomp.o z_trans_b.o z_trans_f.o

tst_pointer: read_line.o tst_pointer.o

check: check_inputs.o check.o fft_init.o fft.o init.o input.o read_line.o read_params.o set_labels.o set_range.o

compare_vel: check_inputs.o compare_vel.o fft_init.o fft.o init.o input.o read_line.o read_params.o set_labels.o

check_inputs.o check.o cl_energy.o compare_vel.o div_check.o fft_init.o force.o initial_field.o init.o input.o input_p.o open_files.o read_field.o read_header.o rhs.o sam.o set_labels.o set_range.o spectra.o tst_input.o tst_transpose.o vel_max.o write_field.o write_header.o write_vel.o x2z_decomp.o xy_trans_b.o xy_trans_f.o z2x_decomp.o z_trans_b.o z_trans_f.o: sam.h

