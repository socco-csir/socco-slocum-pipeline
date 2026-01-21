% GLIDER
%
% Files
%   func_aanderaa_optode_calculation            - function [o] = oxy_cal_glider(p,t,s,bphase,tim,cal_coeff)
%   func_analyze_tso_data                       - function [] = func_analyze_tso_data(glider,ref,comp_type,op,reference_file)
%   func_append_struct                          - function [newstruct] = append_struct(struct1,struct2)
%   func_apply_optode_aanderaa_coeff            - function [o] = oxy_cal_glider(p,t,s,dphase,tim,cal_coeff,t_orig,o_orig)
%   func_apply_optode_calibration               - function [o] = oxy_cal_glider(p,t,s,bphase,tim,cal_coeff)
%   func_cell_flow_speed_from_glider_speed      - function [flow_speed] = func_cell_flow_speed_from_glider_speed(glider_speed,op);
%   func_check_result                           - function [bad_dens_sum,sali_dev_sum,devi3,devi3_first_half,s_field] = func_check_result(op);
%   func_correct_thermal_lag_tomeu              - function [ctd_temperature_in_cell,coefa,coefb] = correctThermalLag(datenum,ctd_temperature,flow_speed,correctionParams,doplot)
%   func_ctd_cell_flow_speed_estimate           - function [flow_speed] = func_ctd_cell_flow_speed_estimate(pressure,pitch,op)
%   func_derive_temp_for_cond_and_get_salinity  - function [] = func_derive_temp_for_cond_and_get_salinity(op)
%   func_determine_tso_ctd_offsets_on_dens      - function [] = func_determine_tso_ctd_offsets_on_dens(op,reference_file)
%   func_determine_tso_ctd_offsets_on_press     - function [] = func_determine_tso_ctd_offsets_on_press(op,reference_file)
%   func_determine_tso_glider_offsets_on_dens   - function [] = func_determine_tso_glider_offsets(op,reference_file)
%   func_determine_tso_glider_offsets_on_press  - function [] = func_determine_tso_glider_offsets(op,reference_file)
%   func_determine_tso_mooring_offsets          - function [] = func_determine_mooring_tso_offsets()
%   func_display_t4c_results                    - function [] = func_display_t4c_results()
%   func_dynamic_model_deviation                - function [res]=dyn_optim(par);
%   func_exp_filter                             - function [y] = func_exp_filter(x,alpha)
%   func_extract_comparison_data_set_on_dens    - function [] = func_extract_comparison_data_set_on_dens(op)
%   func_extract_comparison_data_set_on_press   - function [] = func_extract_comparison_data_set_on_press(op)
%   func_find_pressure_jumps                    - function [jump_ind,jump_correction,down_or_up] = func_find_pressure_jumps(nav_pressure_old,m_present_time,sci_pressure,ctd_time)
%   func_get_deployment_statistics              - function [] = func_get_deployment_statistics()
%   func_glider_acceleration                    - function [a_h,a_v,forces] = func_glider_acceleration(params,temp,rho,ballast,glide_angle,vel,p,aoa,fin)
%   func_grid_glider_profiles                   - function [gridded] = func_grid_glider_profiles(data,mode)
%   func_integrate_glider_acceleration          - function [x,z,aoa,vel,scaling,devi,z_of_p] = func_integrate_glider_acceleration(p,temp,pitch,rho,ballast,params,is_in_surface_mode,batt_pos,fin,air_pump,internal_pressure,shear_lift,no_plot,index,op);
%   func_load_variable_list                     - function [slocum_variable_list,allglider_variable_list] = func_load_variable_lists()
%   func_magdev                                 - function [dev]=magdev(flat,flon,elevkm,year);
%   func_minimization_options                   - BUILDMINIMIZATIONOPTIONS - Builds a set of options for minimization
%   func_optimization_dummy_script              - this is a script that is being used twice
%   func_optimize_dynamical_model               - function [] = func_optimize_dynamical_model()
%   func_optimize_optode_delay                  - function [res] = func_optimize_optode_delay(data,cal_coeff,delays_t,is_aanderaa_coeff)
%   func_optimize_temp_for_cond_params          - function [] = func_optimize_temp_for_cond_params(op)
%   func_optimizer_tomeu                        - function [correctionParams,devi] = adjustThermalLagParams(downcasts, upcasts, firstGuess, mode,op )
%   func_phase_from_optode_oxygen               - function [phase] = func_phase_from_optode_oxygen(p,t_optode,o,cal_coeff)
%   func_plot_delay_function                    - function [] = func_plot_delay_function(params)
%   func_prep4calc                              - function [upcasts,downcasts,meanflow,data] = func_prep4calc(data,op);
%   func_rotate_uv                              - function [ur,vr] = rotate_uv(u,v,ang)
%   func_rough_angle_of_attack                  - function [angle_of_attack,flight_angle] = func_rough_angle_of_attack(pitch)
%   func_sfigure                                - SFIGURE  Create figure window (minus annoying focus-theft).
%   func_simple_aanderaa_delay                  - function [oxygen_undelayed] = func_simple_aanderaa_undelay(cal_coeff,tim,o_orig,t_ctd)
%   func_slocum_best_position                   - function [data] = func_slocum_best_position(data,op);
%   func_slocum_degmin                          - function [decdegs] = func_slocum_degmin(webbdegs)
%   func_slocum_fill_vector                     - function [newdata] = func_slocum_fill_vector(olddata)
%   func_slocum_irregular_geomar_names          - function [] = func_slocum_irregular_geomar_names()
%   func_slocum_rename_all_bd_files             - function [] = func_slocum_rename_all_bd_files()
%   func_slocum_sparse                          - function [dataout] = func_slocum_sparse(datain)
%   func_ts_deviation_tomeu                     - function [sum_of_areas,sum_of_areas_first_half] = func_ts_deviation_tomeu(...
%   func_webbtime2mattime                       - function [datn] = func_webbtime2mattime(mpt,rev)
%   optimize                                    - Optimize general constrained problems using Nelder-Mead algorithm
%   step_00_prepare_processing                  - function [] = step_00_prepare_processing()
%   step_01_copy_bd_to_processing_folder        - function [] = step_01_copy_bd_to_processing_folder()
%   step_02_binary_bd_to_ascii                  - function [] = step_02_binary_bd_to_ascii()
%   step_03_load_ascii_save_mat                 - function [] = step_03_load_ascii_save_mat()
%   step_04_merge_mat_files                     - function [] = step_04_merge_mat_files()
%   step_05_best_position_vector                - function [] = step_05_best_position_vector()
%   step_06_cut_up_down_parts                   - function [] = step_06_cut_up_down_parts()
%   step_07_cleanup_science                     - function [] = step_07_cleanup_science()
%   step_08_correct_pressure                    - function [] = step_08_ccorrect_pressure()
%   step_09_interp_to_1sec                      - function [] = step_09_interp_to_1sec()
%   step_10_prepare_flow_speed                  - function [] = step_10_prepare_flow_speed()
%   step_11_optimize_temp_for_cond_params_dpdt  - function [] = step_11_optimize_temp_for_cond_params_dpdt(no_optim)
%   step_12_prepare_shearlift                   - function [] = step_12_prepare_shearlift()
%   step_13_optimize_dynamical_model            - function [] = step_13_optimize_dynamical_model()
%   step_14_glider_flight_model                 - function [] = step_14_glider_flight_model()
%   step_15_optimize_temp_for_cond_params_model - function [] = step_15_optimize_temp_for_cond_params_model(no_optim)
%   step_16_pick_model_or_dpdt                  - function [] = step_16_pick_model_or_dpdt()
%   step_17_correct_optode                      - function [] = step_17_correct_optode(no_optim)
%   step_18_grid_glider_profiles                - function [] = step_18_grid_glider_profiles()
%   step_19_create_comparison_data              - function [] = step_19_create_comparison_data()
%   step_20_determine_tso_offsets               - function [] = step_20_determine_tso_offsets()
%   step_21_finalize_data                       - function [] = step_21_finalize_data()
%   step_22_summarize_deployment                - function [] = step_22_summarize_deployment()
%   step_23_data_for_others                     - function [] = step_23_data_for_others()
%   step_24_add_suna                            - function [] = step_24_add_suna()
