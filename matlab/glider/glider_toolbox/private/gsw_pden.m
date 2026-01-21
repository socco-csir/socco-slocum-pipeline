function [rho] = gsw_pden(SP,t,p,pr,lon,lat)
% function [rho] = gsw_pden(SP,t,p,pr,lon,lat)
%
% GEOMAR SVN $Id: nsum.m 942 2022-09-08 12:11:01Z gkrahmann@geomar.de $
% 
% near-drop-in replacement for seawater library function sw_pden 
% (potential density) but based on TEOS-10
% near, because it has two additional arguments lon and lat
%
% input  : SP             - practical salinity
%          t              - in situ temperature
%          p              - pressure
%          pr             - reference pressure
%          lon            - longitude
%          lat            - latitude
%
% output : rho            - TEOS-10 potential density in kg/m^3
%
% version 1.0.0  last change 25.01.2023

% G.Krahmann, GEOMAR  Jan 2023

if nargin<5
  lon = 0;
  lat = 0;
end
if nargin<4
  pr = 0;
end

SA = gsw_SA_from_SP(SP,p,lon,lat);
pt = gsw_pt_from_t(SA,t,p,pr);
CT = gsw_CT_from_t(SA,t,p);
rho = gsw_rho(SA,CT,pr);
%rho2 = gsw_pot_rho_t_exact(SA,t,p,pr)
