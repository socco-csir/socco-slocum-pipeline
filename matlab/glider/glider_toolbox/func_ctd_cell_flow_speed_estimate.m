function [glider_speed_in_glider_direction,glider_speed_in_movement_direction,glider_direction,...
  movement_direction,glider_speed_in_cell_direction,angle_of_attack,...
  measured_w] =...
   func_ctd_cell_flow_speed_estimate(pressure,pitch,op,latitude)
% function [flow_speed] = func_ctd_cell_flow_speed_estimate(pressure,pitch,op)
% 
% GEOMAR SVN $Id: func_ctd_cell_flow_speed_estimate.m 942 2022-09-08 12:11:01Z gkrahmann@geomar.de $
%
% function to estimate the speed of the water flow through an unpumped CTD conductivity cell
%
% input  :  pressure              - pressure time series (for Slocum gliders a filtered 
%                                   navigation pressure series is recommended) at 1 second time
%                                   resolution
%           pitch                 - pitch angle time series of the glider in degrees
%           op                    - processing option structure
%
% output :  flow_speed            - flow speed time series in m/s
%           glider_speed_in_glider_direction
%           glider_speed_in_movement_direction
%           glider_direction
%           movement_direction
%           glider_speed_in_cell_direction
%           angle_of_attack
%
% version 6.2.0  last change 25.01.2023

% G.Krahmann, GEOMAR, Aug 2012

% G.Krahmann:
% This is modified from freely available code downloaded from http://www.socib.es/;glider/doco/gliderToolbox/ctdTools/thermalLagTools
% in 2011. See 'Availability of the code' statement in DOI:10.1175/JTECH-D-10-05030.1

% flow speed could be 0                                   GK, 11.12.2012  1-->2
% properly handle variable pitch, include flwospeed calculation here instead
% of in 1sec interpolation                                GK, 19.02.2013  2-->3
% add modifiable paramater op.flow_speed_filter           GK, 20.05.2014  3-->4
% projection of movement speed in cell direction          GK, 28.05.2014  4-->5
% use z instead of p                                      GK, 15.02.2016  5-->6
% added code source statement                             GK, 23.08.2022  6-->6.1.0
% changed from CSIRO seawater to TEOS-10 routines         GK, 25.01.2023  6.1.0-->6.2.0

if iscell(pressure)
  for n=1:length(pressure)
    fprintf(1,'.',[]);
    if iscell(pitch)
      [glider_speed_in_glider_direction{n},glider_speed_in_movement_direction{n},...
       glider_direction{n}, movement_direction{n},glider_speed_in_cell_direction{n},...
       angle_of_attack{n},measured_w{n}] =...
         func_ctd_cell_flow_speed_estimate(pressure{n},pitch{n},op,latitude{n});
    else
      [glider_speed_in_glider_direction{n},glider_speed_in_movement_direction{n},...
       glider_direction{n}, movement_direction{n},glider_speed_in_cell_direction{n},...
       angle_of_attack{n},measured_w{n}] =...
         func_ctd_cell_flow_speed_estimate(pressure{n},pitch,op,latitude{n});
    end
    aoa_estimate{n} = angle_of_attack{n};
  end
  fprintf(1,'\n',[]);
%  save([op.deplname,'_glider_speed'],'glider_speed_in_glider_direction',...
%    'glider_speed_in_movement_direction','glider_direction','movement_direction',...
%    'aoa_estimate','glider_speed_in_cell_direction','op','-v6')
else

  %
  % get a pitch angle that can be used for the flow speed calculation
  %
  if length(pitch)>20
    glider_direction = meanfilt(pitch,10);
  else
    glider_direction = pitch;
  end

  %speedFactorPols = [0.00, 0.00, 0.40;  % 0th order degree
  %                   0.00, 0.03, 0.45;  % 1st order degree
  %                   1.58, 1.15, 0.70]; % 2nd order degree

  mlat = nmedian(latitude(:));
  %z = sw_dpth(pressure,mlat);
  z = -gsw_z_from_p(pressure,mlat);
  good = find(~isnan(z));

  glider_speed_in_glider_direction = nan*z;
  glider_speed_in_movement_direction = nan*z;
  glider_speed_in_cell_direction = nan*z;

  [angle_of_attack,movement_direction] = func_rough_angle_of_attack(glider_direction);

  sinus_movement_direction = sin(movement_direction(good)/180*pi);
  sinus_glider_direction = sin(glider_direction(good)/180*pi);
  cosinus_cell_direction = cos((movement_direction(good)-glider_direction(good))/180*pi);

  deltaDepth = abs(gradient(z(good))); % does not matter if downcast or upcast
  depthRate = deltaDepth;  % for a 1 sec timestep
  measured_w = gradient(z);

  glider_speed_in_movement_direction(good) = abs(depthRate ./ sinus_movement_direction);
  glider_speed_in_glider_direction(good) = abs(depthRate ./ sinus_glider_direction);
  glider_speed_in_cell_direction(good) = glider_speed_in_movement_direction(good).*cosinus_cell_direction;

if 0
  figure(1)
  clf
  plot(glider_speed_in_movement_direction(good))
  hold on
  plot(glider_speed_in_glider_direction(good),'r')
  plot(glider_speed_in_cell_direction(good),'g')
  pause
end

%  if length(glider_speed_in_movement_direction)>op.flow_speed_filter*2
%    surgeSpeed = meanfilt(glider_speed_in_movement_direction,op.flow_speed_filter);
%  else
%    surgeSpeed = glider_speed_in_movement_direction;
%  end

  %selectedDegree = op.flow_speed_degree; % First order approximation, second row of the matrix
  %speedFactor = polyval(speedFactorPols(selectedDegree+1, :), surgeSpeed);
  %speedFactor = polyval(speedFactorPols(:,selectedDegree+1), surgeSpeed);
  %speedFactor = 1;
 
%  flow_speed = speedFactor .* surgeSpeed; % Avoid division by zero
%
%  flow_speed = nans(flow_speed,op.minimum_flow_speed,nan,'==');
%  flow_speed = nans(flow_speed,op.minimum_flow_speed,op.minimum_flow_speed,'<');
    
end

