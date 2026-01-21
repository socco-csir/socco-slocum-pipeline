function [y] = nmedian(x,dim)
% Median value, ignoring NaN.
%
% GEOMAR SVN $Id: nmedian.m 942 2022-09-08 12:11:01Z gkrahmann@geomar.de $
%
% function [y] = nmedian(x,dim)
%
% Same as MEDIAN, but NaN's are ignored.
%
% input  :      x              - data array to be 'median'-ed
%               dim            - dimension over which to calculate the median
%
% output :      y              - 'median'-ed array
%
% version 2.0.0  last change 08.09.2022

% G.Krahmann, GEOMAR, Sep 2022

% rewritten to use Matlab internal function     GK, 08.09.2022 -->2.0.0

if nargin==1
  y = median(x,'omitnan');
elseif nargin==2
  y = median(x,dim,'omitnan');
else
  error('wrong number of arguments')
end
