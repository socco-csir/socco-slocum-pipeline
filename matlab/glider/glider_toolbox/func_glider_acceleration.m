function [a_h,a_v,forces] = func_glider_acceleration(params,temp,rho,ballast,glide_angle,vel,p,aoa,fin,shear_lift,pitch,air_bladder,u,w)
% function [a_h,a_v,forces] = func_glider_acceleration(params,temp,rho,ballast,glide_angle,vel,p,aoa,fin)
% 
% GEOMAR SVN $Id: func_glider_acceleration.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% calculate horizontal and vertical acceleration and forces from the state of the glider
% and the water
%
% input  : params                     - glider parameters such as volume, area etc
%          temp                       - water temperature
%          rho                        - water density
%          ballast                    - state of the glider's ballast pump
%          glide_angle                - angle at which the glider is moving through the water
%          vel                        - glider velocity
%          p                          - pressure
%          aoa                        - angle of attack of the glider
%          fin                        - steering fin angle
%          shear_lift                 - extra lift force from vertical current shear
%          pitch                      - pitch angle of the glider
%          air_bladder                - estimated volume of the airbladder
%          u                          - speed forward ?
%          w                          - speed upward ?
%
% output : a_h                        - horizontal acceleration
%          a_v                        - vertical acceleration
%          forces                     - vector of forces
%
% version 1.1.0   last change 31.01.2023

% G.Krahmann, GEOMAR
% after equations from
% Merckelbach, Lucas; Smeed, David; Griffiths, Gwyn. 2010 Vertical water velocities from 
% underwater gliders. Journal of Atmospheric and Oceanic Technology, 27 (3). 547-563. 
% 10.1175/2009JTECHO710.1 

% change from CSIRO seawater to TEOS library               GK, 31.01.2023  1-->1.1.0

if any(isnan(air_bladder))
  disp('nan airbladder')
  air_bladder = 0;
end

a_w = 3.7;
S = 0.10;
a_h = 2.4;  % this is 1.2 *2 
a_h = a_w * (params.A_h/S) + 1.2;
C_D_1w = 0.78;
C_D_1h = 2.1;
T_0 = 16;

added_mass_factor = 1;

p = p*1e4;

f_b = gsw_grav(params.latitude) * rho * ( params.V_g *...
  (1 - params.epsilon*p + params.alpha_T*(temp-T_0) ) + ballast/1e6 + air_bladder/1e6);
forces(1) = f_b;

f_g = params.m_g * gsw_grav(params.latitude);
forces(2) = f_g;

C_D_0 = params.C_D_0;

C_D = C_D_0 + (C_D_1w + C_D_1h)*aoa^2;

f_d = 1/2 * C_D * rho * S * vel^2;
forces(3) = f_d;

C_L = (a_h + a_w) * aoa;

f_l = 1/2 * C_L * rho * S * vel^2;
forces(4) = f_l;

% glider going forward but with nose down
if glide_angle>0 & glide_angle< pi/2 & pitch<0 & p/1e4<10
  force_h = -sin(glide_angle)*f_l - cos(glide_angle)*f_d;
  force_v = f_b - f_g - cos(glide_angle)*f_l - sin(glide_angle)*f_d;
  jj = 2;
% glider going backward with nose down
elseif glide_angle> pi/2 & pitch<0 & p/1e4<10
  glide_angle = pi - glide_angle;
  force_h = -( sin(glide_angle)*f_l - cos(glide_angle)*f_d );
  force_v = f_b - f_g - cos(glide_angle)*f_l - sin(glide_angle)*f_d + shear_lift;
  jj = 2;
% glider going normal
else
  force_h = sin(glide_angle)*f_l - cos(glide_angle)*f_d;
  force_v = f_b - f_g - cos(glide_angle)*f_l - sin(glide_angle)*f_d + shear_lift;
  jj = 3;
end
forces(5) = shear_lift;
forces(6) = force_h;
forces(7) = force_v;
forces(8) = aoa;
forces(9) = glide_angle/pi*180;
forces(10) = pitch/pi*180;

a_h = force_h / (params.m_g*added_mass_factor);
a_v = force_v / (params.m_g*added_mass_factor);

if air_bladder>0  & pitch<0
%  disp([a_h,a_v,aoa/pi*180,glide_angle/pi*180,f_l,f_d,ballast,air_bladder])
%  pause
end

%disp([f_b,f_g,f_d,f_l,rho/1000,temp,ballast,air_bladder,glide_angle,vel,jj])
%pause(0.01)
