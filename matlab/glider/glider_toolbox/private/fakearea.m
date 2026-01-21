function [b,c,area] = fakearea(s1, t1, s2, t2)
% function [b,c,area] = fakearea(s1, t1, s2, t2)
%
% Calculate an approximate area between two TS lines.
% A good area is calculated when density is monotonous. If it is not
% monotonous a too large area will be calculated. Indeed this is a 
% penalty that we want in this case.
%
% Input :     s1          - salinity of the first line (should be an upcast)
%             t1          - temperature of the first line
%             s2          - salinity of the second line (should be a downcast)
%             t2          - temperature of the second line
%
% Output:     area        - area between the two TS curves
%
% version 1.1.0  last change 26.09.2022

% G.Krahmann, GEOMAR, Mar 2017

% remove NaNs in input data                           GK, 26.09.2022  1-->1.1.0

% dummy output to make it similar to buildPolygon.m
b = [];
c = [];

% remove NaNs
ind = find(~isnan(s1+t1));
s1 = s1(ind);
t1 = t1(ind);
ind = find(~isnan(s2+t2));
s2 = s2(ind);
t2 = t2(ind);

% Make sure input variables are column vectors
s1 = s1(:);
t1 = t1(:);
s2 = s2(:);
t2 = t2(:);

% join both profiles
s = [s1; s2];
t = [t1; t2];

% convert T and S to density and spiciness
% this is not an accurate calculation as we use 
sp = gsw_spiciness0(s,t);
sd = gsw_rho(s,t,0);

% Oversample the curves so that each connecting line has at least one
% point in the density steps. This here could be wrong for 
% large TS reductions.
over_sample_step = 0.003;
spo = interp1([1:length(sp)],sp,[1:over_sample_step:length(sp)]');
sdo = interp1([1:length(sd)],sd,[1:over_sample_step:length(sd)]');

% stick into other variable names
x = sdo;
y = spo;

% integration step in density units
istep = 0.005;

% figure out to which 'index' in density space each spice-density pair belongs
% then for each index get the minimum and maximum value of spice and calculate
% the difference. Multiplied by the density step width this is an area.
% Sum up all areas and we have something like the area between the two
% TS curves. In cases of instabilities the calculated area is an
% overestimation as it includes external areas. But this is actually good
% as this is an implicit penalty for instable parts of the TS diagram.
% see http://blogs.mathworks.com/loren/2008/02/20/under-appreciated-accumarray/
% for an intro into this function
xi = round(x/istep);
xi = xi-min(xi)+1;
xi = [xi,ones(size(xi))];
dminmax = accumarray(xi,y,[],@range);
area = sum(dminmax)*istep;

if 0
  figure(1)
  clf
  plot(sdo,spo,'.')
  ax = axis;
  hold on
  for xx=[nmin(x):istep:nmax(x)]
    ind = find(abs(x-xx)<=istep/2);
    if ~isempty(ind)
      plot(xx*[1,1],[nmin(y(ind)),nmax(y(ind))],'r')
    end
  end
  pause
end
