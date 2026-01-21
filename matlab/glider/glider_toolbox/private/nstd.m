function [y] = nstd(x,dim)
% Standard deviation, ignoring NaN.
%
% GEOMAR SVN $Id: nstd.m 942 2022-09-08 12:11:01Z gkrahmann@geomar.de $
%
% function [y] = nstd(x,dim)
%
% Similar to STD, but NaN's are ignored.
%
% The arguments of NSTD are more restrictive than STD. See 'help std'.
% The standard deviation here is the one that divides by N-1 .
%
% input  :      x              - data array to be 'std'-ed
%               dim            - dimension over which to calculate
%
% output :      y              - 'std'-ed array
%
% version 2.0.0  last change 08.09.2022

% G.Krahmann, GEOMAR, Sep 2022

% rewritten to use Matlab internal function     GK, 08.09.2022 -->2.0.0

if nargin==1
  y = std(x,0,'omitnan');
elseif nargin==2
  y = std(x,0,dim,'omitnan');
else
  error('wrong number of arguments')
end
