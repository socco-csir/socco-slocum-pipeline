function [oxygen_calculated,oxygen_calculated_undelayed,oxygen_undelayed] = oxy_cal_glider(p,t,s,bphase,tim,cal_coeff,t_orig,o_orig)
% function to calibrate glider optode oxygen data
% 
% GEOMAR SVN $Id: func_apply_optode_calibration.m 890 2021-12-15 15:38:48Z gkrahmann@geomar.de $
%
% function [o] = oxy_cal_glider(p,t,s,bphase,tim,cal_coeff)
%
% input : p             - pressure cell vector
%         t             - delayed temperature cell vector
%         s             - salinity cell vector
%         bphase        - Aanderaa b-phase cell vector
%         cal_coeff     - structure with calibration coefficients
%         t_orig        - undelayed temperature cell vector
%         o_orig        - 'factory-calibrated' oxygen
%
% ouput : o             - calibrated oxygen 
%
% version 3       last change 02.03.2017

% G.Krahmann, GEOMAR, Dec 2013
% adapted from
% J.Hahn, GEOMAR, Kiel, Jan. 2012

% added o_orig to be able to apply 'only' the reverse filter       GK, 07.03.2014  1-->2
% removed oxygen_undelayed calc                                    GK, 02.03.2017  2-->3


%
% calibrate oxygen
%
if ~isempty(cal_coeff.calib_param.cal_method)
  cal_method = cal_coeff.calib_param.cal_method;
  foilcoef = cal_coeff.calib_param.foilcoef;
  phasefunstr = cal_coeff.calib_param.phasefunstr;
  beta = cal_coeff.calib_param.beta;
  pcfactor = cal_coeff.calib_param.pcfactor;

  for n=1:length(p)
    out = opt_evalcal(t{n}(:),bphase{n}(:),cal_method,foilcoef,p{n}(:),s{n}(:),phasefunstr,beta,pcfactor);
    oxygen_calculated{n} = out.O2fit';
  end
else
  for n=1:length(p)
    oxygen_calculated{n} = nan*o_orig{n};
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
  % reverse delay the re-calculated o
  good = find(~isnan(oxygen_calculated{n}.*t{n}));
  oxygen_calculated_undelayed{n} = nan*oxygen_calculated{n};
  if length(good)>1
    oxygen_calculated_undelayed{n}(good) = RC_filterback_2points_vardt_vartau(oxygen_calculated{n}(good),...
      tim{n}(good)*86400,tau_O2,tau_fact,t_orig{n}(good),tempref_tauO2);
  end
end
