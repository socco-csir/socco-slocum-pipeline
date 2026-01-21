function [bad_dens_sum,sali_dev_sum,devi3,devi3_first_half,s_field] = func_check_result(op);
% function [bad_dens_sum,sali_dev_sum,devi3,devi3_first_half,s_field] = func_check_result(op);
% 
% GEOMAR SVN $Id: func_check_result.m 557 2019-01-25 10:53:26Z gkrahmann@geomar.de $
%
% check the optimization results on the full data set
%
% input  : op                          - option structure set in 'processing_parameters.m'
%                                        change  op.res  to try different parameters by hand
%                                        Use  op.res = [1,1,0,0] to test the uncorrected data.
%
% output : bad_dens_sum                - sum of all instable densities found in the data set
%          sali_dev_sum                - sum of 
%          devi3                       - sum of all areas underneath the TS-diagram of down-up
%                                        profile pairs (this is the parameter that is optimized)
%          devi3_first_half            - same as before but only for the first half of profiles
%          s_field                     - gridded salinity field
%
% version 2.1.0  last change 25.01.2023

% G.Krahmann, GEOMAR, Sep 2012

% s_field as output                                          GK, 21.12.2012  1-->2
% change function call names                                 
% change from CSIRO seawater to TEOS-10 library              GK, 25.01.2023  2-->2.1.0


%
% load the 1 second interpolated data set, on which all the optimization takes place
%
data = load([op.deplname,'_1sec']);
data2 = load([op.deplname,'_1sec_derived']);
data.pressure_filtered = data2.pressure_filtered;
data.glider_speed_in_cell_direction_dpdt = ...
    data2.glider_speed_in_cell_direction_dpdt;
data.glider_speed_in_cell_direction_model = ...
    data2.glider_speed_in_cell_direction_model;

%
% set parameters to do a check on the full data set
%
op.subsample_profile = 1;
op.min_depth4cond_temp = -10;
op.s_yo_numbers = [];
[downcasts,upcasts,meanflow] = func_prep4calc(data,op);

%
% calculate the deviation result
%
%[devi3,devi3_first_half] = func_ts_deviation_garau(op.res,downcasts,upcasts,op.res*0+1,op.area_function);
[devi3,devi3_first_half] = func_ts_deviation_garau(op.res,downcasts,upcasts,1,op.area_function);

%
% apply the flow speed calculation and derive salinities
%
for m=1:length(data.main_datenum)
  if strcmp(op.use_flow_speed,'dpdt')
    flow_speed = func_cell_flow_speed_from_glider_speed(...
       data.glider_speed_in_cell_direction_dpdt{m},op);
  else
    flow_speed = func_cell_flow_speed_from_glider_speed(...
       data.glider_speed_in_cell_direction_model{m},op);
  end
  good = find(~isnan(data.ctd_temperature{m}+flow_speed));
  data.tempInCell{m} = nan*data.ctd_temperature{m};
  if ~isempty(good)
    data.tempInCell{m}(good) = func_correct_thermal_lag_garau(data.main_datenum{m}(good),...
      data.ctd_temperature{m}(good),flow_speed(good),op.res);
  end
  data.ctd_salinity{m} = gsw_SP_from_R(data.ctd_conductivity{m}/gsw_C3515,data.tempInCell{m},...
    data.pressure_filtered{m});
end

%
% grid the data
%
uncalib = func_grid_glider_profiles(data,2);

%
% prepare gridded data for the calculation of the final deviation measure
%
% since this is just looking for density inversions we use lat=lon=0
%
SA = gsw_SA_from_SP(uncalib.ctd_salinity, uncalib.pressure*ones(1,size(uncalib.ctd_salinity,2)),0,0);
CT = gsw_CT_from_t(SA,uncalib.ctd_temperature,uncalib.pressure*ones(1,size(uncalib.ctd_salinity,2)));
uncalib.d = gsw_rho(SA,CT,uncalib.pressure*ones(1,size(uncalib.ctd_salinity,2)));
%uncalib.d = sw_pden( uncalib.ctd_salinity, uncalib.ctd_temperature,...
%  uncalib.pressure*ones(1,size(uncalib.ctd_temperature,2)),0);

%
% blank profiles that are much longer down or up
% usually these are problem cases such as aborts
%
for n=1:2:size(uncalib.ctd_temperature,2)
  good1 = length(find(~isnan(uncalib.ctd_temperature(:,n))));
  good2 = length(find(~isnan(uncalib.ctd_temperature(:,n+1))));
  if good1>5*good2
    uncalib.ctd_salinity(:,n+1) = nan;
  end
  if good2>5*good1
    uncalib.ctd_salinity(:,n) = nan;
  end
end
  
fprintf(1,'\n',[]);
if op.no_plot~=1
  figure
  pcolor(uncalib.ctd_salinity);
  title(num2str(op.res(1:min([length(op.res),4]))))
end

%
% sum up the different deviations
%
% final measure is  sali_dev_sum + abs(bad_dens_sum) + ts_diagr_first_half
% sali_dev_sum         is a value for 'good looking' in a plot
% abs(bad_dens_sum)    penalizes instabilities
% ts_diagr_first_half  is the minimized area in the TS diagram (excluding potentially bad data)
%
d = diff(uncalib.d);
bad = find(d<0);
bad_dens_sum = sum(d(bad));
ds = uncalib.ctd_salinity(:,2:2:end-1)-uncalib.ctd_salinity(:,3:2:end);
sali_dev_sum = nsum(abs(ds(:)));
disp(['bad_dens_sum = ',num2str(bad_dens_sum)])
disp(['sali_dev_sum = ',num2str(sali_dev_sum)])
disp(['ts_diagr_sum = ',num2str(devi3)])
disp(['ts_diagr_sum_first_half = ',num2str(devi3_first_half)])

if op.no_plot~=1
  figure
  plot(uncalib.ctd_salinity,uncalib.ctd_temperature)

  figure
  plot(uncalib.ctd_salinity(:,1:2:end),uncalib.ctd_temperature(:,1:2:end))
end
s_field = uncalib.ctd_salinity;
