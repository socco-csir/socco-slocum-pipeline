function [datn] = webbtime2mattime(mpt,rev)
% function [datn] = webbtime2mattime(mpt,rev)
%
% convert Webb's time variable to Matlab datenum variable
%
% input  :  mpt        - Webb time data (like m_present_time)
%           rev   ''   - if set to 'reverse' then mattime is converted to webbtime
%
% output :  datn       - Matlab time data (datenum)
%
% version 0.2	last change 16.08.2011

% G.Krahmann, IFM-GEOMAR, May 2008

% allow reverse operation

if nargin<2
  datn = mpt/86400+datenum(1970,1,1,0,0,0);
elseif strcmp(rev,'reverse')
  datn = (mpt-datenum(1970,1,1,0,0,0))*86400;
end

