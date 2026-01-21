function [oxygen_calculated,oxygen_calculated_undelayed] = ...
  func_apply_optode_aanderaa_coeff(p,t,s,dphase,tim,cal_coeff,t_orig,o_orig,t_optode)
% function to apply aanderaa optode coefficients to glider optode phase data
% 
% GEOMAR SVN $Id: func_apply_optode_aanderaa_coeff.m 890 2021-12-15 15:38:48Z gkrahmann@geomar.de $
%
% function [o] = oxy_cal_glider(p,t,s,dphase,tim,cal_coeff,t_orig,o_orig)
%
% input : p             - pressure cell vector
%         t             - delayed temperature cell vector
%         s             - salinity cell vector
%         dphase        - Aanderaa b-phase cell vector
%         cal_coeff     - structure with calibration coefficients
%         t_orig        - undelayed temperature cell vector
%         o_orig        - 'factory-calibrated' oxygen
%
% ouput : o             - calibrated oxygen 
%
% version 2.1.0       last change 17.11.2021

% G.Krahmann, GEOMAR, Apr 2014
% adapted from
% J.Hahn, GEOMAR, Kiel, Jan. 2012

% removed oxygen_undelayed calc                   GK, 02.03.2017  1-->2
% changed header                                  GK, 17.11.2021  2-->2.1.0



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
  % reverse delay the re-calculated o
  good = find(~isnan(oxygen_calculated{n}.*t{n}));
  oxygen_calculated_undelayed{n} = nan*oxygen_calculated{n};
  if length(good)>1
    oxygen_calculated_undelayed{n}(good) = RC_filterback_2points_vardt_vartau(oxygen_calculated{n}(good),...
      tim{n}(good)*86400,tau_O2,tau_fact,t_orig{n}(good),tempref_tauO2);
  end
end
