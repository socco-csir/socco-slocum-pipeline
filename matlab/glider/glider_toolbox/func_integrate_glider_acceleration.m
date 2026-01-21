function [x,z,aoa,vel,scaling,devi,z_of_p,glide_angle,u,w] = func_integrate_glider_acceleration(p,temp,pitch,rho,ballast,params,is_in_surface_mode,batt_pos,fin,air_pump,internal_pressure,shear_lift,no_plot,index,op,bubbles,...
   start_vel,start_glide_angle,start_aoa,start_u,start_w);
% function [x,z,aoa,vel,scaling,devi,z_of_p] = func_integrate_glider_acceleration(p,temp,pitch,rho,ballast,params,is_in_surface_mode,batt_pos,fin,air_pump,internal_pressure,shear_lift,no_plot,index,op);
% 
% GEOMAR SVN $Id: func_integrate_glider_acceleration.m 893 2022-01-10 15:07:34Z gkrahmann@geomar.de $
%
% Integrate the glider acceleration to get the flight path
%
% input  :
%
% output :
%
% version 2.1.0  last change 31.01.2023

% G.Krahmann, GEOMAR 2012

% various changes to implement surface bubble release and inertia from 
% previous yo                                                              GK, 10.01.2019
% bug when no good internal_pressure data available                        GK, 10.01.2022
% change from CSIRO seawater to TEOS library                               GK, 31.01.2022 2.0.1-->2.1.0


substep = 0.1;

ind = [1:substep:length(p)];
ind_full = find(ind==round(ind));
p = interp1([1:length(p)],p,ind);
pitch = interp1([1:length(pitch)],pitch,ind);
rho = interp1([1:length(rho)],rho,ind);
ballast = interp1([1:length(ballast)],ballast,ind);
temp = interp1([1:length(temp)],temp,ind);
is_in_surface_mode = interp1([1:length(is_in_surface_mode)],is_in_surface_mode,ind);
fin = interp1([1:length(fin)],fin,ind);
shear_lift = interp1([1:length(shear_lift)],shear_lift,ind);
air_pump = interp1([1:length(air_pump)],air_pump,ind);
bubbles = interp1([1:length(bubbles)],bubbles,ind);
good = find(~isnan(internal_pressure));
if length(good)<2
  internal_pressure = 0*ind;
else
  internal_pressure = interp1(good,internal_pressure(good),ind,'linear',nan);
end

u = repmat(nan,1,length(p));
w = repmat(nan,1,length(p));
z = repmat(nan,1,length(p));
x = repmat(nan,1,length(p));
ah = repmat(nan,1,length(p));
av = repmat(nan,1,length(p));
vel = repmat(nan,1,length(p));
aoa = repmat(nan,1,length(p));
glide_angle = repmat(nan,1,length(p));
if 0
u(1) = 0;
w(1) = 0;
z(1) = 0;
x(1) = 0;
vel(1) = 0;
aoa(1) = 90/180*pi;
glide_angle(1) = -90/180*pi;
end
pitch = pitch/180*pi;

%
% determine whether it is a down or an up yo
% and, if it is an up yo, whether and when the airbladder is inflated
% In case of active airbladder pump the glider will be in a not well determined flight
% mode. We use the pressure at the time when the air pump is turned on, to determine
% a scaling factor. If the air pump is not active or if the glider is on a descent, we
% simply use the last value to determine the scaling factor as we assume a  regular
% flight characteristic.
%
ind_surface = length(z);
if p(end)>p(1)
  isup = 0;
  ends_in_surface_mode = 0;
else
  isup = 1;
  if any(is_in_surface_mode==1)
    ends_in_surface_mode = 1;
    ind_surface = find(is_in_surface_mode==1);
    ind_surface = ind_surface(1);
  else
    ends_in_surface_mode = 0;
  end
end

%
% try to estimate how much air is in the airbladder
% this will be used as extra buoyancy
%
ind = find(diff(air_pump)>0);
air_bladder = repmat(0,size(air_pump));
dp = [];
if ~isempty(ind) 
  if ind(1)>1
    internal_pressure(1:ind(1)) = internal_pressure(ind(1));
    dp = -(internal_pressure-internal_pressure(ind(1)));
    if op.glider_volume>60/1000 & op.microrider_installed==0 & op.suna_installed==0
      air_bladder = dp/1.013e5*25;   % 25 is the air volume inside short glider in liters
    else
      air_bladder = dp/1.013e5*30;   % 30 is the air volume inside long glider in liters
    end
  end
  air_bladder = nans(air_bladder,0,0,'<');
  air_bladder = nans(air_bladder,0,nan,'==');
  dd = diff(air_bladder);
end
air_bladder = air_bladder*1000;


%
% debugging plot for airbladder volume
%
if 0
figure(6)
subplot(2,1,1)
plot(air_bladder,'x')
ylabel('airbladder volume [ml]')
subplot(2,1,2)
plot(-p)
xlabel('index')
ylabel('pressure')
end

if no_plot~=1
  func_sfigure(1);
  clf
  func_sfigure(2);
  clf
end
devi = [];
count = 0;
u(1) = 0;
u(1) = start_u;
w(1) = 0;
w(1) = start_w;
z(1) = 0;
x(1) = 0;
vel(1) = start_vel;
aoa(1) = 90/180*pi;
aoa(1) = start_aoa;
glide_angle(1) = -90/180*pi;
glide_angle(1) = start_glide_angle;
all_forces = repmat(nan,length(u),10);
for n=2:length(u)
  [a_h,a_v,forces] = func_glider_acceleration(params,temp(n),rho(n),ballast(n),glide_angle(n-1),vel(n-1),p(n),aoa(n-1),...
                                       fin(n-1),shear_lift(n),pitch(n),air_bladder(n),u(n-1),w(n-1));
  ah(n) = a_h;
  av(n) = a_v;
  u(n) = u(n-1) + substep * a_h;
  w(n) = w(n-1) + substep * a_v;
  x(n) = x(n-1) + substep * u(n);
  z(n) = z(n-1) + substep * w(n);
  glide_angle(n) = atan2(w(n),u(n));
  all_forces(n,:) = forces;
  if isup==0
    if z(n)>=0
      z(n) = 0;
      x(n) = x(n-1);
      u(n) = 0;
      w(n) = 0;
      glide_angle(n) = -90/180*pi;
    end
  else
    if z(n)<=0
      z(n) = 0;
      x(n) = x(n-1);
      u(n) = 0;
      w(n) = 0;
      glide_angle(n) = -90/180*pi;
    end
  end
  vel(n) = sqrt(u(n)^2+w(n)^2);
  aoa(n) = glide_angle(n) - pitch(n);
  if glide_angle(n)>pi/2 & pitch(n)<0          % this is the case when the glider moves backward
    aoa(n) = (pi-glide_angle(n)) + pitch(n);
  end
end

%
% adjust flight model to observed start and end pressures
%
z_of_p = gsw_z_from_p(p,params.latitude);
z = z - (z(1)-z_of_p(1));
if z(ind_surface)-z(1)~=0
  scaling = (z_of_p(ind_surface)-z_of_p(1))/(z(ind_surface)-z(1));
else
  scaling = 1;
end
z = (z-z(1)) * scaling + z(1);
vel = vel * scaling;
x = x * scaling;
z = nans(z,0,0,'>');

%
% plot model integration result
%
if no_plot~=1

  %
  % plot modelled and observed depth
  %
  func_sfigure(1);
  plot(z)
  hold on
  plot(z_of_p,'r')
  plot(ind_surface,0,'x')
  xlabel('model index')
  ylabel('depth')
  title(['yo:',int2str(index),'   b: model   r: observed'])
  drawnow

  %
  % plot angles of attached
  func_sfigure(2);
  plot(all_forces(:,9),'b')
  hold on
  plot(all_forces(:,10),'r')
  xlabel('model index')
  ylabel('angle [degrees]')
  title(['yo:',int2str(index),'   b: pitch angle   r: movement angle'])
  drawnow

  func_sfigure(3);
  plot(batt_pos)
  xlabel('model index')
  ylabel('battery position [inches]')
  title(['yo:',int2str(index)])
  drawnow

  func_sfigure(5);
  clf
  if z(1)>z(end)
    ind = find(z<-10);
    if ~isempty(ind)
      ind = [1:ind(1)];
    end
  else
    ind = find(z<-10);
    if ~isempty(ind)
      ind = [ind(end):length(z)];
    end
  end
  if ~isempty(ind)
    plot(z(ind))
    hold on
    plot(z_of_p(ind),'r')
    ax = axis;
    plot(pitch(ind)*180/pi,'g')
    plot(-bubbles(ind)/10,'k')
    axis(ax);
  end
  xlabel('model index')
  ylabel('depth')
  title(['yo:',int2str(index),'   b: model   r: observed  g: pitch  Zoom'])
  drawnow
%  pause
%  sfigure(2);
%  plot(x)
%  xlabel('model index')
%  ylabel('distance')
%  title('b: model   r: observed')
%  drawnow
end
if 0
figure(3)
clf
subplot(5,1,1)
plot(z)
subplot(5,1,2)
plot(ah)
grid on
subplot(5,1,3)
plot(u)
subplot(5,1,4)
plot(pitch/pi*180)
hold on
plot(glide_angle/pi*180,'r')
plot(aoa/pi*180,'g')
grid on
subplot(5,1,5)
plot(av)
grid on
pause
end

z = z(ind_full);
x = x(ind_full);
aoa = aoa(ind_full);
vel = vel(ind_full);
z_of_p = z_of_p(ind_full);
