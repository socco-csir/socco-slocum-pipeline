function [] = step_24_add_suna()
% function [] = step_24_add_suna()
%
% GEOMAR SVN $Id: step_23_data_for_others.m 557 2019-01-25 10:53:26Z gkrahmann@geomar.de $
%
% add calibrated suna data to the final files
%
% version 2.1.0  last change 16.06.2021

% G.Krahmann, GEOMAR,  Sep 2020

% change mat file name for OPUS/SUNA code         GK, 08.12.2020  1-->2.0.0
% catch nox duplicates                            GK, 16.06.2021   2.0.0-->2.1.0


%
% look for and load processing parameters
%
if exist('./processing_parameters.m') | exist('.\processing_parameters.m')
  processing_parameters
else
  error(['could not find file ''processing_parameters.m'' in current folder ',pwd])
end


%
% load the calibrated SUNA data
%
load([op.deplname,'_nox_calibrated.mat'])


%
% remove duplicates in nox data
%
[light_datenum,ind] = unique(light_datenum);
light_nitrate_um_calibrated = light_nitrate_um_calibrated(ind);


%
% load the 1 second glider data
%
load([op.deplname,'_1sec'])


%
% interpolate calibrated SUNA data to 1 second 
%
for n=1:length(main_datenum)

  good = find(~isnan(light_nitrate_um_calibrated+light_datenum) &...
    light_datenum>=main_datenum{n}(1) & light_datenum<=main_datenum{n}(end));
  if length(good)>1
   nitrate{n} = interp1(light_datenum(good),light_nitrate_um_calibrated(good),main_datenum{n});
  else
    nitrate{n} = nan*main_datenum{n};
  end

end


%
% grid SUNA data
%
data.nitrate = nitrate;
data.pressure = best_pressure;
[gridded] = func_grid_glider_profiles(data,5);


%
% add nitrate data to other final data
%
save([op.deplname,'_final_1sec.mat'],'nitrate','-append')
nitrate = gridded.nitrate;
save([op.deplname,'_final_gridded.mat'],'nitrate','-append')
