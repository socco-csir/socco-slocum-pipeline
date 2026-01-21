function O2conc_sal=pO2toO2conc(pO2,T,S,P_atm)
% function O2conc_sal=pO2toO2conc(pO2,T,S,P_atm)
% 
% GEOMAR SVN $Id: pO2toO2conc.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% calculate O2 concentration in umol/l from pO2 / mbar with salinity and
% atm. pressure / mbar (or GTD pressure, respectively) correction
%
% part of optcalc-toolbox
% H. Bittig, IFM-GEOMAR
% 31.03.2011

% pH2O / atm
pH2O = watervapor(T,S);
% atm. pressure / atm
atm_press=P_atm/1013.25; 
% theta0
th0=1-(0.999025+0.00001426.*T-0.00000006436.*T.^2);
% O2 solubility / umol/l
oxy_sol=O2solubility(T,S);
% pressure corrected O2 solubility / umol atm / l
oxy_sol_pc=oxy_sol.*atm_press.*(((1-pH2O./atm_press).*(1-th0.*atm_press))./((1-pH2O).*(1-th0)));
% oxygen concentration in umol/l (salinity corrected)
O2conc_sal=(pO2.*oxy_sol_pc)./((atm_press-pH2O).*0.20946.*1013.25);
