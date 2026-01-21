function [newdata] = func_slocum_fill_vector(olddata)
% function [newdata] = func_slocum_fill_vector(olddata)
%
% GEOMAR SVN $Id: func_slocum_fill_vector.m 188 2016-07-04 15:51:56Z gkrahmann@geomar.de $
% 
% fill a (sparse) Slocum data vector so that all spots are occupied
%
% input  : olddata                - old data vector that contains values only when they change
%
% output : newdata                - new data vector in which all elements contain values
%
% version 1  last change 22.03.2013

% G.Krahmann, GEOMAR, Mar 2013

%
% check for a sparse vector
%
newdata = olddata;
if issparse(newdata)
  newdata = func_slocum_sparse(newdata);
end

%
% fill vector
%
value = nan;
for n=1:length(newdata)
  if value~=newdata(n) & ~isnan(newdata(n))
    value = newdata(n);
  end
  newdata(n) = value;
end

