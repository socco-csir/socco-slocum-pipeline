function [oxygen_undelayed] = func_simple_aanderaa_undelay(cal_coeff,tim,o_orig,t_ctd)
% function [oxygen_undelayed] = func_simple_aanderaa_undelay(cal_coeff,tim,o_orig,t_ctd)
% 
% GEOMAR SVN $Id: func_simple_aanderaa_delay.m 591 2019-03-20 11:42:59Z gkrahmann@geomar.de $
%
% function to apply a simple reverse time filter to the original oxygen data
%
% input : cal_coeff          - calibration coefficient structure
%         tim                - time vector
%         o_orig             - 'factory-calibrated' oxygen
%         t_ctd              - CTD temperature 
%
% ouput : oxygen_undelayed   - undelayed oxygen
%
% version 1       last change 02.03.2017

% G.Krahmann, GEOMAR, Mar 2017


%
% apply reverse filter to optode data
%

%%% parameters for backfolding
tau_O2 = cal_coeff.calib_param.tau_O2;
tau_fact = cal_coeff.calib_param.tau_fact;
tempref_tauO2 = cal_coeff.calib_param.tempref_tauO2;

for n=1:length(t_ctd)
  % reverse delay the uncalibrated o
  good = find(~isnan(o_orig{n}.*t_ctd{n}));
  oxygen_undelayed{n} = nan*o_orig{n};
  if length(good)>1
    oxygen_undelayed{n}(good) = RC_filterback_2points_vardt_vartau(o_orig{n}(good),...
      tim{n}(good)*86400,tau_O2,tau_fact,t_ctd{n}(good),tempref_tauO2);
  end
end
