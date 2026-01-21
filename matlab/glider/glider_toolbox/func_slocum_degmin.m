function [decdegs] = func_slocum_degmin(webbdegs)
% function [decdegs] = func_slocum_degmin(webbdegs)
% 
% GEOMAR SVN $Id: func_slocum_degmin.m 942 2022-09-08 12:11:01Z gkrahmann@geomar.de $
%
% convert Webb's position numbers into decimal degrees
%
% input  :	webbdegs	- Webb's position number
%
% output :	decdegs		- decimal degrees
%
% version 2.1.0  last change 08.09.2022

% G.Krahmann, IFM-GEOMAR, Oct 2006

% allow for cell structures                                              GK, 23.08.2012  0.1-->2
% replace a leftover call to the old function name 'degmin' with the new
% function name 'func_slocum_degmin'                                     GK, 08.09.2022  2-->2.1.0

if iscell(webbdegs)
  for n=1:length(webbdegs)
    decdegs{n} = func_slocum_degmin(webbdegs{n});
  end
else
  si = sign(webbdegs);
  webbdegs = abs(webbdegs);
  degs = floor(webbdegs/100);
  mins = webbdegs - degs*100;
  decdegs = degs + mins/60;
  decdegs = si .* decdegs;

  bad = find(abs(decdegs)>360);
  if ~isempty(bad)
    decdegs(bad) = 0;
  end
end
