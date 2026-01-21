function [ctd_temperature_in_cell,coefa,coefb] = correctThermalLag(datenum,ctd_temperature,flow_speed,correctionParams,doplot)
% function [ctd_temperature_in_cell,coefa,coefb] = correctThermalLag(datenum,ctd_temperature,flow_speed,correctionParams,doplot)
% 
% GEOMAR SVN $Id: func_correct_thermal_lag_tomeu.m 942 2022-09-08 12:11:01Z gkrahmann@geomar.de $
%
%CORRECTTHERMALLAG - CTDs Thermal lag correction.
% This function receives as input parameter a CTD profile, and applies a
% thermal lag correction on it.
% The correction applied uses the recursive scheme described in
% Morison, J., R. Andersen, N. Larson, E. D'Asaro, and T. Boyd, 1994:
% The Correction for Thermal-Lag Effects in Sea-Bird CTD Data.
% Journal of Atmospheric and Oceanic Technology, vol. 11, pages 1151ï¿½1164.
%
% Syntax:
%    correctedProfileData = correctThermalLag(basicProfileData)
%    correctedProfileData = correctThermalLag(basicProfileData, correctionParams)
%    [correctedProfileData, correctionParams] = correctThermalLag(basicProfileData)
%
% Inputs:
%    basicProfileData - A profile structure*
%    correctionParams - The set of parameters to be used in the correction*
%
% Outputs:
%    correctedProfileData - A profile structure*
%    correctionParams - The set of parameters used in the correction*
%
% * Profile structure: A struct that contains several fields,
%   all of them column vectors with the same length:
%   - ptime: Present time instant at which this row was collected
%   - depth: Depth (pressure in decibars) measured by the CTD
%   - temp: Temperature measured by the CTD
%   - cond: Conductivity measured by the CTD
%   - pitch: Pitch angle of the glider (optional)
%   The output profile has the same information of the input plus
%   two fields with the corrected profile properties:
%   - condOutCell: corrected conductivity, removing the effects of the
%     temperature difference between the outer and inner parts of the
%     conductivity cell.
%   - tempInCell: corrected temperature, that is, the temperature of
%     the water mass lying inside the conductivity cell.
%
%   From this information, the user can choose which one of the two
%   corrections to use in order to compute salinity:
%   - Combine 'condOutCell' with 'temp' (Expected values outside of the
%     conductivity cell).
%   - Combine 'cond' with 'tempInCell' (Expected values inside of the
%     conductivity cell).
%
% * Correction parameters: A vector of four elements, consisting in
%   alpha_offset, alpha_slope, tau_offset and tau_slope.
%   These parameters are used to compute alpha and tau,
%   the amplitude and time constant respectively, which are inversely
%   proportional to the flow speed.
%
% Example:
%   correctedProfileData = correctThermalLag(basicProfileData) corrects the
%   profile information contained in the input 'basicProfileData',
%   using Morison parameters
%
%   correctedProfileData = correctThermalLag(basicProfileData, correctionParams)
%   corrects the profile information contained in 'basicProfileData',
%   using correctionParams as the parameters to be used
%   during the correction.
%
%   [correctedProfileData, correctionParams] = correctThermalLag(basicProfileData)
%   corrects the profile information contained in 'basicProfileData', and
%   provides in the output the correction parameters used for the correction.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ADJUSTTHERMALLAGPARAMS
%
% Author: Bartolome Garau
% Work address: Parc Bit, Naorte, Bloc A 2Âºp. pta. 3; Palma de Mallorca SPAIN. E-07121
% Author e-mail: tgarau@socib.es
% Website: http://www.socib.es
% Creation: 17-Feb-2011

% G.Krahmann:
% This is modified from freely available code downloaded from http://www.socib.es/;glider/doco/gliderToolbox/ctdTools/thermalLagTools
% in 2011. See 'Availability of the code' statement in DOI:10.1175/JTECH-D-10-05030.1

% GEOMAR version 1.0.0  last change 23.08.2022

% added code source statement                              GK, 23.08.2022  -->1.0.0


if 0
if length(correctionParams)==6
  correctionParams = correctionParams(1:4);
end
if length(correctionParams)==7
  correctionParams = correctionParams(1:5);
end
if length(find(isnan(correctionParams)))==1
  correctionParams = correctionParams(find(~isnan(correctionParams)));
end
end


  % Extract the information contained in the profile

if length(correctionParams)==5
  alpha_offset = correctionParams(1);
  alpha_slope  = correctionParams(2);
  tau_offset = correctionParams(3);
  tau_slope  = correctionParams(4);
  t_offset  = correctionParams(5);
elseif length(correctionParams)==4
  alpha_offset = correctionParams(1);
  alpha_slope  = correctionParams(2);
  tau_offset = correctionParams(3);
  tau_slope  = correctionParams(4);
  t_offset = 0;
elseif length(correctionParams)==3
  alpha_offset = correctionParams(1);
  alpha_slope  = 0;
  tau_offset = correctionParams(2);
  tau_slope  = 0;
  t_offset = correctionParams(3);
elseif length(correctionParams)==2
  alpha_offset = correctionParams(1);
  alpha_slope  = 0;
  tau_offset = correctionParams(2);
  tau_slope  = 0;
  t_offset = 0;
end

if t_offset~=0
  ind = [1:length(ctd_temperature)];
  good = find(~isnan(ctd_temperature));
  if length(good)>=2
    ctd_temperature(good) = interp1(ind(good),ctd_temperature(good),ind(good)+t_offset,...
      'linear','extrap');
  end
end
  
% Some initial precomputations
deltaTemp    =     gradient(ctd_temperature) ;
    
tau   =   tau_offset +   tau_slope ./ sqrt(flow_speed);
alpha = alpha_offset + alpha_slope ./      flow_speed ;
bad = find(alpha==0);
if ~isempty(bad)
  alpha(bad) = eps;
end

if any(alpha)<0
  disp(['alpha : ',num2str(min(alpha))])
  disp(['flowspeed : ',num2str(max(flow_speed))])
end
if any(alpha)>1
  disp(['alpha : ',num2str(max(alpha))])
  disp(['flowspeed : ',num2str(max(flow_speed))])
end
if any(tau)<0
  disp(['tau : ',num2str(min(tau))])
  disp(['flowspeed : ',num2str(max(flow_speed))])
end


% Relation between a and b coefficients with respect to alpha and tau
% 0.5 for the Nyquist freq of our sampling freq if 1 Hz
coefa = 0.5 * 4 .* alpha .* tau ./ (1 + 4 .* tau * 0.5);
%    coefa = nans(coefa,1.5,1.5,'>');
coefb = 1 - 2 .* coefa ./ alpha;

tempCorrection = zeros(size(ctd_temperature));

dummy1 = coefa.*deltaTemp;
dummy1 = nans(dummy1,0,nan,'==');
for depthLevel = 1:length(ctd_temperature)-1,
  if datenum(depthLevel+1)-datenum(depthLevel)>10/86400 | isnan(ctd_temperature(depthLevel)) |...
      isnan(ctd_temperature(depthLevel+1))
    tempCorrection(depthLevel+1) = 0;
  else
    tempCorrection(depthLevel+1) = ...
      - coefb(depthLevel) .* tempCorrection(depthLevel)...
         + dummy1(depthLevel);
        % + coefa * (ctd_temperature(depthLevel+1)-ctd_temperature(depthLevel));
  end
end

ctd_temperature_in_cell  = ctd_temperature - tempCorrection;
