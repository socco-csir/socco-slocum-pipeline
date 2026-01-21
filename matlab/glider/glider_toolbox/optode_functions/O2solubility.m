function out=O2solubility(T,S)
% function out=O2solubility(T,S)
% 
% GEOMAR SVN $Id: O2solubility.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% calculate oxygen solubilty / umol/l
%
% part of optcalc-toolbox
% H. Bittig, IFM-GEOMAR
% 31.03.2011

sca_T = scaledT(T);
out=((exp(2.00856+3.224.*sca_T+3.99063.*sca_T.^2+4.80299.*sca_T.^3+0.978188.*sca_T.^4+...
    1.71069.*sca_T.^5+S.*(-0.00624097-0.00693498.*sca_T-0.00690358.*sca_T.^2-0.00429155.*sca_T.^3)...
    -0.00000031168.*S.^2))./0.022391903);
