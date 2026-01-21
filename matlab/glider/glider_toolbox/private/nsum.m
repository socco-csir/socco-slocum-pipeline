function [y] = nsum(x,dim)
% Sum up values, ignoring NaN.
%
% GEOMAR SVN $Id: nsum.m 942 2022-09-08 12:11:01Z gkrahmann@geomar.de $
%
% function [y] = nsum(x,dim)
%
% Same as SUM, but NaN's are ignored.
%
% input  :      x              - data array to be summed
%               dim            - dimension over which to sum
%
% output :      y              - summed array
%
% version 2.0.0  last change 08.09.2022

% G.Krahmann, GEOMAR, Sep 2022

% rewritten to use Matlab internal function     GK, 08.09.2022 -->2.0.0

if nargin==1
  y = sum(x,'omitnan');
elseif nargin==2
  y = sum(x,dim,'omitnan');
else
  error('wrong number of arguments')
end
