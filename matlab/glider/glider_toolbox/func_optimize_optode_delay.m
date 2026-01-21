function [res] = func_optimize_optode_delay(data,cal_coeff,delays_t,is_aanderaa_coeff)
% function [res] = func_optimize_optode_delay(data,cal_coeff,delays_t,is_aanderaa_coeff)
% 
% GEOMAR SVN $Id: func_optimize_optode_delay.m 303 2017-03-06 08:52:37Z gkrahmann@geomar.de $
%
%
%
% input  : 
%
% output :
%
% version 2  last change 02.03.2017

% G.Krahmann, GEOMAR  2012

% introduced independent calculation for simple Aanderaa delay  GK, 02.03.2017  1-->2

%
% first we handle the case of only delaying the optode
% and jump back out when finished
%
if is_aanderaa_coeff==100
  data.oxygen_undelayed = func_simple_aanderaa_delay(cal_coeff,data.main_datenum,data.oxygen,data.ctd_temperature);
  if length(delays_t)>1
    gridded = func_grid_glider_profiles(data,4);
    res{1} = gridded;
  else
    res{1}.oxygen_undelayed = data.oxygen_undelayed;
  end
  return
end

%
% this is the case of provided calibration coefficients
%
for m=1:length(delays_t)

  for n=1:length(data.ctd_temperature)

    % grab the original CTD temperatures
    delayed_ctd_temperature{n} = data.ctd_temperature{n};

    % figure out which data points are ok
    good = find(~isnan(delayed_ctd_temperature{n}));

    % fill all gaps between the first and last good point
    if length(good)>1
      delayed_ctd_temperature{n}(good(1):good(end)) = interp1(good,delayed_ctd_temperature{n}(good),...
        [good(1):good(end)]);

      % now filter this complete part of the CTD temperature series
      delayed_ctd_temperature{n}(good(1):good(end)) = RC_filter_2points_vardt_vartau(...
        delayed_ctd_temperature{n}(good(1):good(end))',[1:length(good)]',delays_t(m),0)';
    end
  end

  if ~isempty(cal_coeff.foil_coeff)
    [data.oxygen_calculated,data.oxygen_calculated_undelayed] =...
      func_apply_optode_aanderaa_coeff(data.ctd_pressure,delayed_ctd_temperature,...
      data.ctd_salinity_dpdt,data.diffphase,data.main_datenum,cal_coeff,...
      data.ctd_temperature,data.oxygen,data.oxygen_temperature);
  else
    [data.oxygen_calculated,data.oxygen_calculated_undelayed] =...
      func_apply_optode_calibration(data.ctd_pressure,delayed_ctd_temperature,...
      data.ctd_salinity_dpdt,data.diffphase,data.main_datenum,cal_coeff,...
      data.ctd_temperature,data.oxygen);
  end

  if length(delays_t)>1
    gridded = func_grid_glider_profiles(data,4);

    res{m} = gridded;

  else

    res.oxygen_calculated = data.oxygen_calculated;
    res.oxygen_calculated_undelayed = data.oxygen_calculated_undelayed;

  end

end
