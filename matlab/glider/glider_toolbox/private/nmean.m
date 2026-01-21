function [y] = nmean(x,dim)
% Average or mean value, ignoring NaN.
%
% GEOMAR SVN $Id: nmean.m 942 2022-09-08 12:11:01Z gkrahmann@geomar.de $
%
% function [y] = nmean(x,dim)
%
% Same as MEAN, but NaN's are ignored.
%
% input  :      x              - data array to be averaged
%               dim            - dimension over which to average
%
% output :      y              - averaged array
%
% version 2.0.0  last change 08.09.2022

% G.Krahmann, GEOMAR, Sep 2022

% rewritten to use Matlab internal function     GK, 08.09.2022 -->2.0.0

if nargin==1
  y = mean(x,'omitnan');
elseif nargin==2
  y = mean(x,dim,'omitnan');
else
  error('wrong number of arguments')
end
