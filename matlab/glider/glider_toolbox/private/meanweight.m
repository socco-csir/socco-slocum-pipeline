function [m]=meanweight(data,weight,dim,mode);
% calculates mean of data weighted with weight
%
% GEOMAR SVN $Id: meanweight.m 668 2020-02-11 10:19:35Z gkrahmann@geomar.de $
%
% function [m]=meanweight(data,weight,dim,[mode])
%
% accepts NaNs in data and weight
%
% input  :  data                - data array
%           weight              - weights array
%           dim                 - which dimension
%                                 if 0, mean of all points
%
% output :  m                   - mean value(s)
%
% uses :	nsum.m
%
% version 1.2.0		last change 04.10.2007

% Gerd Krahmann, IfM Kiel, Mar 1993
% added modes, G.Krahmann IfM Kiel 27.7.1994
% added compatibility to MATLAB 5	G.Krahmann, LODYC Paris, 
% removed small bug			G.Krahmann, LDEO 	1.1.0-->1.1.1
% removed compatibility to Matlab 4	GK, Sep 2007		1.1.1-->1.2.0


if nargin<4
  mode = 0;
end

if dim==0
  mode=1;
end

a=find(isnan(data)+isnan(weight));
if ~isempty(a)
  weight(a)=zeros(length(a),1);
end

if mode==1
  sweight=sum(weight(:));
  if sweight~=0
    m=nsum(data(:).*weight(:))/sweight;
  else
    m=nan;
  end
else
  sweight=sum(weight,dim);
  m=nsum(data.*weight,dim)./sweight;
end
