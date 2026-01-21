function [flow_speed] = func_cell_flow_speed_from_glider_speed(glider_speed,op);
% function [flow_speed] = func_cell_flow_speed_from_glider_speed(glider_speed,op);
% 
% GEOMAR SVN $Id: func_cell_flow_speed_from_glider_speed.m 942 2022-09-08 12:11:01Z gkrahmann@geomar.de $
%
% calculate speed of flow through conductivity cell from glider speed
%
% input  : glider_speed               - cell array or vector of glider speeds
%          op                         - processing option structure
%
% output : flow_speed                 - cell array or vector of cell flow speeds
%
% version 3.1.0  last change 23.08.2022

% G.Krahmann, GEOMAR, May 2014

% G.Krahmann:
% This is modified from freely available code downloaded from http://www.socib.es/;glider/doco/gliderToolbox/ctdTools/thermalLagTools
% in 2011. See 'Availability of the code' statement in DOI:10.1175/JTECH-D-10-05030.1

% add more text                                        GK, 22.02.2017  1-->2
% added default for flow_speed_filter                  GK, 02.09.2020  2-->3.0.0
% added code source statement                          GK, 23.08.2022  3.0.0-->3.1.0

%
% check whether we have a cell array. If yes, recurse
%
if iscell(glider_speed)
  for n=1:length(glider_speed)
    flow_speed{n} = func_cell_flow_speed_from_glider_speed(glider_speed{n},op);
  end
  return
end

%
% This is taken from Garau. But in their code there is an
% error in the next line as they use a row instead of a column
% of this matrix.
%
speedFactorPols = [0.00, 0.00, 0.40;  % 0th order degree
                   0.00, 0.03, 0.45;  % 1st order degree
                   1.58, 1.15, 0.70]; % 2nd order degree

speedFactor = polyval(speedFactorPols(:,op.flow_speed_degree+1), ...
      glider_speed);


if ~isfield(op,'flow_speed_filter')
  op.flow_speed_filter = 30;
end


flow_speed = speedFactor .* glider_speed;

if length(flow_speed)>op.flow_speed_filter
  flow_speed = meanfilt(flow_speed,op.flow_speed_filter,'triang-onesided');
end


flow_speed = nans(flow_speed,op.minimum_flow_speed,nan,'==');
flow_speed = nans(flow_speed,op.minimum_flow_speed,op.minimum_flow_speed,'<');

flow_speed = flow_speed(:)';


