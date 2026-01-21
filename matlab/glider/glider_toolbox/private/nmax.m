function [y,ind] = nmax(x,dim)
% Maximum value, ignoring NaN.
%
% GEOMAR SVN $Id: nmax.m 943 2022-09-12 14:39:42Z gkrahmann@geomar.de $
%
% function [y,ind] = nmin(x,dim)
%
% Same as MAX, but NaN's are ignored.
%
% input  :      x              - data array to be maximized
%               dim            - dimension over which to maximize
%
% output :      y              - maximized array
%               ind            - index of minimum
%
% version 2.0.1  last change 11.09.2022

% G.Krahmann, GEOMAR, Sep 2022

% rewritten to use Matlab internal function     GK, 08.09.2022 -->2.0.0
% wrong usage of max                            GK, 11.09.2022 2.0.0-->2.0.1

if nargin==1
  [y,ind] = max(x,[],'omitnan');
elseif nargin==2
  [y,ind] = max(x,[],dim,'omitnan');
else
  error('wrong number of arguments')
end
