function pw=watervapor(T,S)
% function pw=watervapor(T,S)
% 
% GEOMAR SVN $Id: watervapor.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% calculating pH2O / atm after Weiss and Price 1980
% T in °C
%
% part of optcalc-toolbox
% H. Bittig, IFM-GEOMAR
% 31.03.2011

pw=(exp(24.4543-(67.4509*(100./(T+273.15)))-(4.8489*log(((273.15+T)./100)))-0.000544.*S));
