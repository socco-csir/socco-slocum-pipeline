function [correctionParams,devi] = func_optimizer_garau(downcasts, upcasts, firstGuess, mode,op )
% function [correctionParams,devi] = func_optimizer_garau(downcasts, upcasts, firstGuess, mode,op )
% 
% GEOMAR SVN $Id: func_optimizer_tomeu.m 942 2022-09-08 12:11:01Z gkrahmann@geomar.de $
%
%ADJUSTTHERMALLAGPARAMS - CTDs Thermal lag parameters adjustment.
% This function receives as input parameters 'downcast' and 'upcast',
% two profiles which are supposed to be measuring the same water column
% in different directions.
% Based on the assumption that both profiles should be as similar as
% possible, it finds the thermal lag parameters related to alpha and tau, such that
% the area in a TS diagram between the corrected profiles is minimum.
%
% Syntax: correctionParams = adjustThermalLagParams(downcast, upcast)
%
% Inputs:
%    downcast - Profile structure*
%    upcast - Profile structure*
%
% * Profile structure: A struct that contains several fields,
%   all of them column vectors with the same length:
%   - ptime: Present time instant at which this row was collected
%   - depth: Depth (pressure in decibars) measured by the CTD
%   - temp: Temperature measured by the CTD
%   - cond: Conductivity measured by the CTD
%   - pitch: Pitch angle of the glider (optional)
%
% Outputs:
%    correctionParams - vector containing the correction parameters
%   used to correct the profile thermal lag.
%
% Example:
%    correctionParams = adjustThermalLagParams(downcast, upcast)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: CORRECTTHERMALLAG, FMINCON
%
% Author: Bartolome Garau
% Work address: Parc Bit, Naorte, Bloc A 2Âºp. pta. 3; Palma de Mallorca SPAIN. E-07121
% Author e-mail: tgarau@socib.es
% Website: http://www.socib.es
% Creation: 17-Feb-2011
%
% Adapted to GEOMAR processing   G.Krahmann GEOMAR 2012
%
% GEOMAR version 1.1.0  last change 25.01.2023

% G.Krahmann:
% This is modified from freely available code downloaded from http://www.socib.es/;glider/doco/gliderToolbox/ctdTools/thermalLagTools
% in 2011. See 'Availability of the code' statement in DOI:10.1175/JTECH-D-10-05030.1

% added code source statement                                         GK, 23.08.2022  -->1.0.0
% change function name                                                GK, 25.01.2023  1.0.0-->1.1.0

if nargin<4
  mode = 1;
end

     % Minimize the area between the two profiles in the TS Diagram
    options = func_minimization_options;
    
    if length(firstGuess)==5
      lowerBound = [0.01,eps,eps,eps,-2];
      upperBound = [0.2, op.minimum_flow_speed, 100, 100,2];
    elseif length(firstGuess)==4
      lowerBound = [0.01,eps,eps,eps];
      upperBound = [0.2, op.minimum_flow_speed, 100, 100];
    elseif length(firstGuess)==3
      lowerBound = [0.01,eps,-2];
      upperBound = [1, 100, 2];
    elseif length(firstGuess)==2
      lowerBound = [0.01,1];
      upperBound = [1, 100];
    end
    
    for nn=1:length(firstGuess)
      firstGuess(nn) = nmax([lowerBound(nn),firstGuess(nn)]);
      firstGuess(nn) = nmin([upperBound(nn),firstGuess(nn)]);
    end
    devi(1) = func_ts_deviation_garau(firstGuess,downcasts,upcasts,op.ts_data_reduction,op.area_function);
    if mode==1
      correctionParams = fmincon(@(x)func_ts_deviation_garau(x,downcasts,upcasts,...
        op.ts_data_reduction,op.area_function), firstGuess,...
	[], [], [], [], lowerBound, upperBound, [], options);
    else
      correctionParams = optimize(@(x)func_ts_deviation_garau(x,downcasts,upcasts,...
        op.ts_data_reduction,op.area_function), firstGuess,...
	lowerBound, upperBound,[],[],[],[],[],[],options);
    end

    devi(2) = func_ts_deviation_garau(correctionParams,downcasts,upcasts,op.ts_data_reduction,op.area_function);
