function options = buildMinimizationOptions
%BUILDMINIMIZATIONOPTIONS - Builds a set of options for minimization
% Optional file header info (to give more details about the function than in the H1 line)
% Optional file header info (to give more details about the function than in the H1 line)
% 
% GEOMAR SVN $Id: func_minimization_options.m 942 2022-09-08 12:11:01Z gkrahmann@geomar.de $
%
% Syntax: options = buildMinimizationOptions
%
% Inputs: none
%
% Outputs:
%    options - an optimset output
%
% Example:
%    options = buildMinimizationOptions
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OPTIMSET
%
% Author: Bartolome Garau
% Work address: Parc Bit, Naorte, Bloc A 2Âºp. pta. 3; Palma de Mallorca SPAIN. E-07121
% Author e-mail: tgarau@socib.es
% Website: http://www.socib.es
% Creation: 17-Feb-2011
%
% Adapted for GEOMAR processing.  G.Krahmann GEOMAR 2012

% GEOMAR version 1.0.0  last change 23.08.2022

% G.Krahmann:
% This is modified from freely available code downloaded from http://www.socib.es/;glider/doco/gliderToolbox/ctdTools/thermalLagTools
% in 2011. See 'Availability of the code' statement in DOI:10.1175/JTECH-D-10-05030.1

% added code source statement                                      GK, 23.08.2022 -->1.0.0


tv = ver('optim');
if isempty(tv)
    defaultOptions = optimset('fminsearch');
else
    defaultOptions = optimset('fmincon');
end

    
    options = optimset(defaultOptions, ...
        'Display',    'iter',          ...
        'LargeScale', 'off',           ...
        'Algorithm',  'active-set',    ...
        'TolFun',     1e-6,            ...
        'TolX',       1e-7);%,            ...
%        'TolCon',     1e-10,            ...
%        'TolFun',     1e-4,            ...
%        'TolCon',     1e-5,            ...
%        'TolX',       1e-5);%,            ...
%        'Plotfcns',   []); % {@optimplotfval, @optimplotfirstorderopt, @optimplotx});

end
