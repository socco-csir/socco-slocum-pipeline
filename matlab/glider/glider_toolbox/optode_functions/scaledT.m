function out=scaledT(in)
% function out=scaledT(in)
% 
% GEOMAR SVN $Id: scaledT.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% calculate scaled temperature
%
% part of optcalc-toolbox
% H. Bittig, IFM-GEOMAR
% 31.03.2011

out=log((298.15-in)./(273.15+in));
