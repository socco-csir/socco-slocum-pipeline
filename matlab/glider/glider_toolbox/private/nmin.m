function [y,ind] = nmin(x,dim)
% Minimum value, ignoring NaN.
%
% GEOMAR SVN $Id: nmin.m 943 2022-09-12 14:39:42Z gkrahmann@geomar.de $
%
% function [y,ind] = nmin(x,dim)
%
% Same as MIN, but NaN's are ignored.
%
% input  :      x              - data array to be minimized
%               dim            - dimension over which to minimize
%
% output :      y              - minimized array
%               ind            - index of minimum
%
% version 2.0.1  last change 11.09.2022

% G.Krahmann, GEOMAR, Sep 2022

% rewritten to use Matlab internal function     GK, 08.09.2022 -->2.0.0
% wrong usage of min                            GK, 11.09.2022 2.0.0-->2.0.1

if nargin==1
  [y,ind] = min(x,[],'omitnan');
elseif nargin==2
  [y,ind] = min(x,dim,'omitnan');
else
  error('wrong number of arguments')
end
