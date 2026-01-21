function [oxygen_aanderaa] = oxy_cal_glider(p,t,s,bphase,dphase,foil_coeff)
% function to apply an optode calibration to glider optode data
% 
% GEOMAR SVN $Id: func_aanderaa_optode_calculation.m 890 2021-12-15 15:38:48Z gkrahmann@geomar.de $
%
% function [o] = oxy_cal_glider(p,t,s,bphase,tim,cal_coeff)
%
% input : p             - pressure cell vector
%         t             - delayed temperature cell vector
%         s             - salinity cell vector
%         bphase        - Aanderaa b-phase cell vector
%         dphase        - Aanderaa d-phase cell vector
%         foil_coeff    - structure with foil coefficients
%
% ouput : o             - Aanderaa formula oxygen 
%
% version 1.1.0       last change 17.11.2021

% G.Krahmann, GEOMAR, Apr 2014
% adapted from
% J.Hahn, GEOMAR, Kiel, Jan. 2012

% changed header                                  GK, 17.11.2021  1-->1.1.0


%
% calculate oxygen
%
if length(foil_coeff)==20
  oxygen_aanderaa = optcalcO2(temp,phase,foilcoef,'aanderaa',sal,P_atm,P_dbar,pcfactor)
elseif length(foil_coeff)==21
  oxygen_aanderaa = optcalcO2(temp,phase,foilcoef,'5x5b',sal,P_atm,P_dbar,pcfactor)
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
  good = find(~isnan(oxygen_calibrated{n}.*t{n}));
  oxygen_calibrated_undelayed{n} = nan*oxygen_calibrated{n};
  if length(good)>1
    oxygen_calibrated_undelayed{n}(good) = RC_filterback_2points_vardt_vartau(oxygen_calibrated{n}(good),...
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
