function h = func_sfigure(h)
% SFIGURE  Create figure window (minus annoying focus-theft).
% 
% GEOMAR SVN $Id: func_sfigure.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% Usage is identical to figure.
%
% Daniel Eaton, 2005
%
% See also figure

% The original function was downloaded from
% the Matlab File Exchange in 2010
% https://de.mathworks.com/matlabcentral/fileexchange/8919-smart-silent-figure?s_tid=prof_contriblnk
% no license is provided with the file (checked 2023).
% The function has been modified following comments in the Matlab File Exchange.

if nargin>=1 
	if ishandle(h)
		set(0, 'CurrentFigure', h);
	else
		h = figure(h);
	end
else
	h = figure;
end
if nargout == 1
        varargout{1} = h;
end
