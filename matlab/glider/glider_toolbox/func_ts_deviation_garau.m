function [sum_of_areas,sum_of_areas_first_half] = func_ts_deviation_garau(...
  params,downcasts,upcasts,ts_data_reduction,area_function)
%function [sum_of_areas,sum_of_areas_first_half] = func_ts_deviation_garau(...
%  params,downcasts,upcasts,ts_data_reduction,area_function)
% 
% GEOMAR SVN $Id: func_ts_deviation_tomeu.m 942 2022-09-08 12:11:01Z gkrahmann@geomar.de $
%
% Function to determine the area in TS-space enclosed by two temperature and salinity profiles
% when the salinity profile is derived from lagged temperature data.
%
% If more than one profile pair is given as part of the input structure, then the 
% sum of the areas of all pairs is given.
%
% Input  : params                           - optimization parameters
%          downcasts                        - structure of downcast data, can be a vector of structures
%          upcasts                          - similar as downcasts
%          ts_data_reduction       [1]      - to speed up the very slow computation, one can reduce
%                                             the number of data points in the analyzed TS-diagrams.
%                                             Set to higher than 1 to skip values.
%          area_function           [1]      - 2: fakearea.m   1:buildPolygon.m
%
% Output : sum_of_areas                     - area or sum of areas underneath the TS-diagrams
%          sum_of_areas_first_half          - area or sum of areas underneath the TS-diagrams
%                                             but only for the first half of the profiles given
%                                             This one is helpfull for cases of very slow gliders because of
%                                             bio-fouling.
%
% version 4.3.0  last change 25.01.2023

% G.Krahmann, GEOMAR, Sep 2012
%
% coded according to ideas of Bartolome Garau and colleagues
% Garau, Bartolomé & Ruiz, Simón & Zhang, Weifeng & Pascual, Ananda & 
% Heslop, Emma & Kerfoot, John & Tintoré, J.. (2011). 
% Thermal Lag Correction on Slocum CTD Glider Data. 
% Journal of Atmospheric and Oceanic Technology - 
% J ATMOS OCEAN TECHNOL. 28. 1065-1071. 10.1175/JTECH-D-10-05030.1. 

% G.Krahmann:
% This is modified from freely available code downloaded from http://www.socib.es/;glider/doco/gliderToolbox/ctdTools/thermalLagTools
% in 2011. See 'Availability of the code' statement in DOI:10.1175/JTECH-D-10-05030.1

% using CTD pressure now                                GK, 12.12.2012  1-->2
% added fakearea.m for optimization                     GK, 06.03.2017  2-->3
% remove dependency on own rms.m                        GK, 10.03.2020  3-->4.0.0
% added reference                                       GK, 28.09.2021  4.0.0-->4.1.0
% fix function name                                     GK, 12.11.2021  4.1.0-->4.1.1
% added code source statement                           GK, 23.08.2022  4.1.1-->4.2.0
% changed function name
% changed from CSIRO seawater to TEOS-10 library        GK, 25.01.2023  4.2.0-->4.3.0

%
% check input arguments
%
if nargin<5
  area_function = 1;
end

%
% initialize result
%
sum_of_areas = 0;
sum_of_areas_first_half = 0;
nprofiles = length(downcasts);


%
% loop over all downcasts
%
for n=1:length(downcasts)

  % extract a single down/upcast pair
  if length(downcasts)==1
    downcast1 = downcasts;
    upcast1 = upcasts;
  else
    downcast1 = downcasts(n);
    upcast1 = upcasts(n);
  end

  % Correct both profiles with the same parameters
  [downcast1.ctd_temperature_in_cell,down_coefa,down_coefb] =...
    func_correct_thermal_lag_garau( downcast1.datenum,...
    downcast1.ctd_temperature, downcast1.flow_speed, params);
  [upcast1.ctd_temperature_in_cell,up_coefa,up_coefb] =...
    func_correct_thermal_lag_garau( upcast1.datenum,...
    upcast1.ctd_temperature, upcast1.flow_speed, params);

if 0
if length(params)==4
  figure(1)
  clf
  subplot(2,2,1)
  plot(down_coefa)
  subplot(2,2,2)
  plot(down_coefb)
  subplot(2,2,3)
  plot(up_coefa)
  subplot(2,2,4)
  plot(up_coefb)
  drawnow
  pause
end
end
    
  % Compute salinity and temperature of each profile
  temperatureDowncast = downcast1.ctd_temperature;
  cndrDowncast = downcast1.ctd_conductivity_ratio;
  salinityDowncast = gsw_SP_from_R(cndrDowncast, downcast1.ctd_temperature_in_cell, downcast1.ctd_pressure);
            
  temperatureUpcast = upcast1.ctd_temperature;
  cndrUpcast = upcast1.ctd_conductivity_ratio;
  salinityUpcast = gsw_SP_from_R(cndrUpcast, upcast1.ctd_temperature_in_cell, upcast1.ctd_pressure);

  if length(params)==5 | length(params)==3
    good = find(~isnan(temperatureDowncast));
    temperatureDowncast(good) = interp1(good,temperatureDowncast(good),good+params(end),'linear','extrap');
    good = find(~isnan(temperatureUpcast));
    temperatureUpcast(good) = interp1(good,temperatureUpcast(good),good+params(end),'linear','extrap');
  end
    
  % catch problem cases
  % this should not happen
  if length(find(isnan(salinityDowncast)))~=length(find(isnan(temperatureDowncast)))
    disp('1 size difference in deviation')
    keyboard
  end
  if length(find(isnan(salinityUpcast)))~=length(find(isnan(temperatureUpcast)))
    disp('2 size difference in deviation')
    keyboard
  end

  % std of the S gradients
  dsalinityDowncast = diff(salinityDowncast);
  dsalinityUpcast = diff(salinityUpcast);
  good_down = find(~isnan(dsalinityDowncast));
  good_up = find(~isnan(dsalinityUpcast));
  single_gradient1 = rms(dsalinityDowncast(good_down)) + rms(dsalinityUpcast(good_up));
  single_gradient1 = single_gradient1*20;

  % apply the data reduction to speed up the calculation
  salinityDowncast = salinityDowncast(1:ts_data_reduction:end);
  salinityUpcast = salinityUpcast(1:ts_data_reduction:end);
  temperatureDowncast = temperatureDowncast(1:ts_data_reduction:end);
  temperatureUpcast = temperatureUpcast(1:ts_data_reduction:end);

  % do the actual area calculations
  good_down = find(~isnan(salinityDowncast+temperatureDowncast));
  good_up = find(~isnan(salinityUpcast+temperatureUpcast));

  % area of the TS curves
  if area_function==1
    [s, t, single_area1] = buildPolygon( ...
            salinityUpcast(good_up), temperatureUpcast(good_up), ...
            salinityDowncast(good_down), temperatureDowncast(good_down));
  elseif area_function==2
    [s, t, single_area1] = fakearea( ...
            salinityUpcast(good_up), temperatureUpcast(good_up), ...
            salinityDowncast(good_down), temperatureDowncast(good_down));
  end
  %single_top = abs(nmean(salinityDowncast(good_down(1:10)))-nmean(salinityUpcast(good_up(end-10:end))));

  % area of the surrounding hull
  %x = [salinityDowncast(good_down),salinityUpcast(good_up)];
  %y = [temperatureDowncast(good_down),temperatureUpcast(good_up)];
  %[k,single_area2] = convhull(x,y);

  %single_area = single_area1+single_area2;
  %single_area = single_area1+single_area2+single_gradient1;
  single_area = single_area1+single_gradient1;
  %single_area = single_area1;

if n==0
  func_sfigure(1);
  clf
  plot(salinityDowncast,temperatureDowncast,salinityUpcast,temperatureUpcast)
  drawnow
  disp(params)
  pause(0.1)
end

  % sum up the areas
  sum_of_areas = sum_of_areas + single_area;
  %sum_of_areas = sum_of_areas + single_area + single_top;
  if n<=nprofiles/2
    sum_of_areas_first_half = sum_of_areas_first_half + single_area;
  end
end
