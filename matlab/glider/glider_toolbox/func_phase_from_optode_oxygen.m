function [phase] = func_phase_from_optode_oxygen(p,t_optode,o,cal_coeff)
% function [phase] = func_phase_from_optode_oxygen(p,t_optode,o,cal_coeff)
% 
% GEOMAR SVN $Id: func_phase_from_optode_oxygen.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $

%
% function to reverse calculate phase information from optode oxygen and its coefficients
%
% input : p             - pressure cell vector or single value
%         t_optode      - delayed temperature cell vector or single value
%         o             - oxygen by the optode
%         cal_coeff     - structure with calibration coefficients
%
% ouput : phase         - optode phase that delivers the original oxygen reading
%
% version 1       last change 15.07.2014

% G.Krahmann, GEOMAR, Jul 2014


%
% handle cell input
%
if iscell(p)
  for n=1:length(p)
    phase{n} = func_phase_from_optode_oxygen(p{n},t_optode{n},o{n},cal_coeff);
  end
  return
end
    

%
% get Aanderaa coefficients
%
if length(cal_coeff.foil_coeff(:))==20
  modeltype = 'aanderaa';
elseif length(cal_coeff.foil_coeff(:))==21
  modeltype = '5x5b';
end
P_atm = 1013.25;
pcfactor = 3.2;


for n=1:length(p)
  oxygen_calculated{n} = optcalcO2(t{n},dphase{n},cal_coeff.foil_coeff,...
    modeltype,s{n},P_atm,p{n},pcfactor);
if 0
  oxygen_calculated2{n} = optcalcO2(t{n},dphase{n},cal_coeff.foil_coeff,...
    modeltype,0*s{n},P_atm,0*p{n},pcfactor);
  oxygen_calculated3{n} = optcalcO2(t{n},dphase{n},cal_coeff.foil_coeff,...
    modeltype,0*s{n}+35,P_atm,0*p{n},pcfactor);
  oxygen_calculated4{n} = optcalcO2(t_optode{n},dphase{n},cal_coeff.foil_coeff,...
    modeltype,0*s{n}+35,P_atm,0*p{n},pcfactor);
  oxygen_calculated5{n} = optcalcO2(t_optode{n},dphase{n},cal_coeff.foil_coeff,...
    modeltype,0*s{n},P_atm,0*p{n},pcfactor);
end
end


%
% apply reverse filter to optode data
%

%%% parameters for backfolding
tau_O2 = cal_coeff.calib_param.tau_O2;
tau_fact = cal_coeff.calib_param.tau_fact;
tempref_tauO2 = cal_coeff.calib_param.tempref_tauO2;

for n=1:length(p)
  % reverse delay the calibrated o
  good = find(~isnan(oxygen_calculated{n}.*t{n}));
  oxygen_calculated_undelayed{n} = nan*oxygen_calculated{n};
  if length(good)>1
    oxygen_calculated_undelayed{n}(good) = RC_filterback_2points_vardt_vartau(oxygen_calculated{n}(good),...
      tim{n}(good)*86400,tau_O2,tau_fact,t_orig{n}(good),tempref_tauO2);
  end
  % reverse delay the uncalibrated o
  good = find(~isnan(o_orig{n}.*t{n}));
  oxygen_undelayed{n} = nan*o_orig{n};
  if length(good)>1
    oxygen_undelayed{n}(good) = RC_filterback_2points_vardt_vartau(o_orig{n}(good),...
      tim{n}(good)*86400,tau_O2,tau_fact,t_orig{n}(good),tempref_tauO2);
  end
end
