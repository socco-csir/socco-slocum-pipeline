function [] = step_18_grid_glider_profiles()
% function [] = step_18_grid_glider_profiles()
% 
% GEOMAR SVN $Id: step_18_grid_glider_profiles.m 557 2019-01-25 10:53:26Z gkrahmann@geomar.de $
%
% grid 1 sec interpolated glider profile data to 1dbar steps
%
% version 5  last change 15.01.2019

% G.Krahmann, GEOMAR,  Aug 2012

% modified handling of derived variables    GK, 16.07.2014  1-->2
% fixed problem with picking of S           GK, 01.02.2017  2-->3
% add processing diary, remove v6 mat comp  GK, 30.01.2018  3-->4
% extended diary                            GK, 15.01.2019  4-->5


%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 18  at  ',datestr(now)])
disp('grid glider data')
diary off


%
% look for and load processing parameters
%
if exist('./processing_parameters.m') | exist('.\processing_parameters.m')
  processing_parameters
else
  error(['could not find file ''processing_parameters.m'' in current folder ',pwd])
end


%
% load data
%
data = load([op.deplname,'_1sec']);
data2 = load([op.deplname,'_1sec_derived']);
data2.best_pressure = data.best_pressure;


%
% do the gridding
%
gridded = func_grid_glider_profiles(data,1);
gridded2 = func_grid_glider_profiles(data2,1);
fn = fieldnames(gridded2);
for n=1:length(fn)
  gridded = setfield(gridded,fn{n},getfield(gridded2,fn{n}));
end


%
% cleanup data near the surface, which stems from the non-perfect separation
% of down and up yos
%
for n=1:size(gridded.main_datenum,2)
  if all(isnan(gridded.ctd_temperature(11:end,n)))
    gridded.ctd_temperature(:,n) = nan;
    gridded.ctd_conductivity(:,n) = nan;
    gridded.ctd_pressure(:,n) = nan;
    gridded.ctd_salinity_dpdt(:,n) = nan;
    gridded.ctd_salinity(:,n) = nan;
    gridded.dphase(:,n) = nan;
    gridded.bphase(:,n) = nan;
    gridded.chlorophyll(:,n) = nan;
    gridded.turbidity(:,n) = nan;
    gridded.cdom(:,n) = nan;
    gridded.oxygen(:,n) = nan;
    gridded.nitrate(:,n) = nan;
    gridded.water_w(:,n) = nan;
    if isfield(gridded,'ctd_salinity_model')
      gridded.ctd_salinity_model(:,n) = nan;
    end
    if isfield(gridded,'aanderaa_oxygen_calculated_undelayed')
      gridded.aanderaa_oxygen_calculated_undelayed(:,n) = nan;
    end
    if isfield(gridded,'aanderaa_oxygen_calculated_undelayed_filtered')
      gridded.aanderaa_oxygen_calculated_undelayed_filtered(:,n) = nan;
    end
    if isfield(gridded,'aanderaa_oxygen_calculated')
      gridded.aanderaa_oxygen_calculated(:,n) = nan;
    end
    if isfield(gridded,'aanderaa_oxygen_undelayed')
      gridded.aanderaa_oxygen_undelayed(:,n) = nan;
    end
    if isfield(gridded,'aanderaa_oxygen_undelayed_filtered')
      gridded.aanderaa_oxygen_undelayed_filtered(:,n) = nan;
    end
    if isfield(gridded,'geomar_oxygen_calculated_undelayed')
      gridded.geomar_oxygen_calculated_undelayed(:,n) = nan;
    end
    if isfield(gridded,'geomar_oxygen_calculated_undelayed_filtered')
      gridded.geomar_oxygen_calculated_undelayed_filtered(:,n) = nan;
    end
    if isfield(gridded,'geomar_oxygen_calculated')
      gridded.geomar_oxygen_calculated(:,n) = nan;
    end
    if isfield(gridded,'geomar_oxygen_undelayed')
      gridded.geomar_oxygen_undelayed(:,n) = nan;
    end
    if isfield(gridded,'geomar_oxygen_undelayed_filtered')
      gridded.geomar_oxygen_undelayed_filtered(:,n) = nan;
    end
    gridded.oxygen_temperature(:,n) = nan;
  end
end



%
% copy the salinity and oxygen data
%
o_is_found = 0;
data2.oxygen_cal = 'none';
if isfield(gridded,'geomar_oxygen_calculated_undelayed_filtered')
  if any(~isnan(gridded.geomar_oxygen_calculated_undelayed_filtered(:)))
    gridded.good_oxygen = gridded.geomar_oxygen_calculated_undelayed_filtered;
    data2.good_oxygen = data2.geomar_oxygen_calculated_undelayed_filtered;
    o_is_found = 1;
    data2.oxygen_cal = 'geomar';
  end
end
if isfield(gridded,'aanderaa_oxygen_calculated_undelayed_filtered') & o_is_found==0
  if any(~isnan(gridded.aanderaa_oxygen_calculated_undelayed_filtered(:)))
    gridded.good_oxygen = gridded.aanderaa_oxygen_calculated_undelayed_filtered;
    data2.good_oxygen = data2.aanderaa_oxygen_calculated_undelayed_filtered;
    o_is_found = 1;
    data2.oxygen_cal = 'aanderaa';
  end
end
if isfield(gridded,'aanderaa_oxygen_undelayed_filtered') & o_is_found==0
  if any(~isnan(gridded.aanderaa_oxygen_undelayed_filtered(:)))
    gridded.good_oxygen = gridded.aanderaa_oxygen_undelayed_filtered;
    data2.good_oxygen = data2.aanderaa_oxygen_undelayed_filtered;
    o_is_found = 1;
    data2.oxygen_cal = 'simple_aanderaa';
  end
end
if o_is_found==0
  gridded.good_oxygen = gridded.oxygen;
  data2.good_oxygen = data.oxygen;
end
load([op.deplname,'_best_optimizer'],'flow_type')
if strcmp(flow_type,'_dpdt')
  gridded.good_salinity = gridded.ctd_salinity_dpdt;
  data2.good_salinity = data2.ctd_salinity_dpdt;
else
  gridded.good_salinity = gridded.ctd_salinity_model;
  data2.good_salinity = data2.ctd_salinity_model;
end


%
% save the data
%
save([op.deplname,'_gridded'],'-struct','gridded')
save([op.deplname,'_1sec_derived'],'-struct','data2')


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 18  at  ',datestr(now)])
disp(' ')
diary off

