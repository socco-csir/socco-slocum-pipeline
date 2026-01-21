function [y,y2] = func_exp_filter(x,alpha)
% function [y] = func_exp_filter(x,alpha)
% 
% GEOMAR SVN $Id: func_exp_filter.m 311 2017-03-22 10:41:52Z gkrahmann@geomar.de $
%
% exponential 'time' scale filter
% 
% y_n+1 = (1-alpha)*x_n+1 + alpha*y_n   or something like that
%
% input  :  x             - vector that is to be filtered
%           alpha         - coefficient
%
% output :  y             - filtered vector
%
% half life of weights is about  (2/alpha - 1)/2.8854  steps
%
% version 1  last change 23.04.2013

% G.Krahmann, GEOMAR April 2013


y = nan*x;
y(1) = x(1);
for n=2:length(x)
  y(n) = alpha*x(n)+(1-alpha)*y(n-1);
end

fac = (1-alpha).^[0:400];
x = x([ones(400,1);[1:length(x)]']);
fac = fac/sum(fac);
y2 = conv(x,fac);
y2 = y2(401:end);
