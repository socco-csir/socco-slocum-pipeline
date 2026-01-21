function [f_out] = RC_filter_2points_vardt_vartau(f_in,time,tau0,tau_fact,vardep_tau,vardep_tau0)
%% function [f_out] =
%% RC_filter_2points_vardt_vartau(f_in,time,tau0,tau_fact,vardep_tau,vardep_tau0)
% 
% GEOMAR SVN $Id: RC_filter_2points_vardt_vartau.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% function to calculate low pass filtered data assuming an exponential
% behaviour with time constant tau0 including a variable dependence 
%
% NOTE: - 'f_in(1)' is associated to be in an equilibrium state
%       - 'tau' and 'time' need to refer to the same unit e.g.: [s] or [days]
%       - 'tau' can be variable dependent, e.g. temperature dependent
%                   --> when slope of -0.4 s/degC is assumed, i.e.
%                           tau = tau(20deg) - 0.4 * (temp-20)
%                   --> if so, then set 'tau_fact', 'vardep_tau',
%                                       'vardep_tau0'
%
% INPUT:    f_in       ... data (nx1)
%           time       ... time corresponding to data ; (nx1)  [tau] = [time]
%           tau0       ... time constant at certain value vardep_tau0 (1x1)
%           tau_fact   ... linear dependence of tau (e.g. on temperature) is
%                          assumed with this factor ; (1x1)
%                          unit = [time]/[vardep_tau] , e.g. [s/degC]
%           vardep_tau ... variable, which affects tau ; corresponding to
%                          data vector ; (nx1) unit = [vardep_tau]
%           vardep_tau0 ... certain value of vardep_tau, where tau0 is
%                           given ; (1x1) unit = [vardep_tau]
% 
%   ------>     - tau_fact, vardep_tau, vardep_tau0 are optional; if not given,
%               then 'tau_fact' = 0 and no variable dependence is given
%               - if 'tau_fact' is set to 0 --> no variable dependence is given
%
% OUTPUT:   f_out   ... low pass filtered data (nx1)
%
%
%
%
% J.Hahn, GEOMAR
% 26.01.2012
%
% revised 24.02.2012, header information in 'NOTE' changed, J.Hahn
% catch tau_dep==0                                          GK, 07.03.2014
%

%%

len = length(f_in);   % length of data set
dt = diff(time);      % time difference between two data points

if exist('tau_fact','var')
    if isempty(tau_fact)|(tau_fact==0)
        tau_fact = 0;
        vardep_tau = zeros(len,1);
        vardep_tau0 = 0;
    end
else
    tau_fact = 0;
    vardep_tau = zeros(len,1);
    vardep_tau0 = 0;
end

tau_dep = tau0 + tau_fact*(vardep_tau(1:end-1)-vardep_tau0);    % set variable dependence of time constant

ef = repmat(0,size(dt));
good = find(tau_dep~=0);
if ~isempty(good)
  ef(good) = exp(-dt(good)./tau_dep(good));
end

f_out(1) = f_in(1);

for i=2:len
    f_out(i) = f_in(i)*(1-ef(i-1)) + f_out(i-1)*ef(i-1);
end
