function O2sat_sal=O2conctoO2sat(O2conc_sal,T,S,P_atm)
% function pO2conc_sal=O2sattoO2conc(O2sat_sal,T,S,P_atm)
% 
% GEOMAR SVN $Id: O2conctoO2sat.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% calculate O2 saturation in % from O2 concentration / umol/l with salinity S and
% atmospheric/GTD pressure / mbar correction 
%
% part of optcalc-toolbox
% H. Bittig, IFM-GEOMAR
% 31.03.2011

pH2O = watervapor(T,S); % pH2O / atm
atm_press=P_atm/1013.25; % atm. pressure / atm
th0=1-(0.999025+0.00001426.*T-0.00000006436.*T.^2); % theta0
oxy_sol=O2solubility(T,S); % O2 solubility / umol/l
% pressure corrected O2 solubility / umol atm / l
oxy_sol_pc=oxy_sol.*atm_press.*(((1-pH2O./atm_press).*(1-th0.*atm_press))./((1-pH2O).*(1-th0)));
% oxygen concentration in umol/l (salinity corrected)
O2sat_sal=O2conc_sal./oxy_sol_pc.*100;
