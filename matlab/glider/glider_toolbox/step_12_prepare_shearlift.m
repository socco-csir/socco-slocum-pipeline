function [] = step_12_prepare_shearlift()
% function [] = step_12_prepare_shearlift()
% 
% GEOMAR SVN $Id: step_12_prepare_shearlift.m 557 2019-01-25 10:53:26Z gkrahmann@geomar.de $
%
% prepare shearlift data in case a current profile was given
%
% version 4  last change 30.01.2018

% G.Krahmann, GEOMAR  Aug 2012

% properly handle variable pitch                        GK, 19.02.2013  2-->3
% add processing diary,  remove v6 mat comp             GK, 30.01.2018  3-->4

%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 12  at  ',datestr(now)])
disp('calculate additional lift force from glider inertia in sheared flow')
diary off




%
% look for and load processing parameters
%
if exist('./processing_parameters.m') | exist('.\processing_parameters.m')
  processing_parameters
else
  error(['could not find file ''processing_parameters.m'' in current folder ',pwd])
end


%
% load 1 sec data
%
load([op.deplname,'_1sec'],'pressure','pitch','main_datenum','heading');
load([op.deplname,'_1sec_derived'])
  

%
% create current profile, if an m-file is present
%
if isfield(op,'shear_lift_factor')
  shear_lift_factor = op.shear_lift_factor;
else
  shear_lift_factor = 0;
end
if exist('make_current_profile.m')


  make_current_profile

  load([op.deplname,'_currentprofile'])
  if ~exist('time_prof')
    time_prof = 0;
  end
  if size(u_prof,1)==1
    u_prof = u_prof';
    v_prof = v_prof';
    z_prof = z_prof';
  end
else
  u_prof = [0;0;0];
  v_prof = [0;0;0];
  z_prof = [-100;0;2000];
  time_prof = 0;
end
if size(u_prof,2)==1
  u_prof = [u_prof,u_prof];
  v_prof = [v_prof,v_prof];
  z_prof = [z_prof,z_prof];
  ex = 1;
else
  ex = 0;
end
[dummy,g_u_prof] = gradient(u_prof);
[dummy,g_v_prof] = gradient(v_prof);
[dummy,g_z_prof] = gradient(z_prof);
u_shear = g_u_prof./g_z_prof;
v_shear = g_v_prof./g_z_prof;

if ex==1
  u_shear = u_shear(:,1);
  v_shear = v_shear(:,1);
  z_prof = z_prof(:,1);
end


%
% derive the shear in the direction of the glider movement at any time of the deployment
%
for n=1:length(pressure)
  u_shear1 = nan*pressure{n};
  v_shear1 = nan*pressure{n};
  [dummy,time_ind] = nmin(abs(time_prof-nmean(main_datenum{n})));
  good = find(~isnan(pressure{n}));
  u_shear1(good) = interp1(z_prof(:,time_ind),u_shear(:,time_ind),pressure{n}(good));
  v_shear1(good) = interp1(z_prof(:,time_ind),v_shear(:,time_ind),pressure{n}(good));
  local_shear{n} = rotuv((270+heading{n})/180*pi,u_shear1,v_shear1);
  bad = find(isnan(local_shear{n}));
  if ~isempty(bad)
    local_shear{n}(bad) = 0;
  end
  shear_lift{n} = shear_lift_factor * local_shear{n} ...
    .*glider_speed_in_movement_direction_dpdt{n}/nmax(glider_speed_in_movement_direction_dpdt{n});
  shear_lift{n} = meanfilt(shear_lift{n},10);

if 0
  figure(1)
  clf
  subplot(3,2,1)
  disp(nmean(heading{n}))
  plot(shear_lift{n},-pressure{n})
  subplot(3,2,2)
  plot(pressure{n})
  subplot(3,2,3)
  plot(u_prof(:,time_ind),-z_prof(:,time_ind))
  subplot(3,2,4)
  plot(v_prof(:,time_ind),-z_prof(:,time_ind))
  subplot(3,2,5)
  plot(u_shear(:,time_ind),-z_prof(:,time_ind))
  subplot(3,2,6)
  plot(v_shear(:,time_ind),-z_prof(:,time_ind))
  hold on
  plot(local_shear{n},-pressure{n},'r')
  pause
end


end


%
% save the data
%
data = load([op.deplname,'_1sec_derived']);
data.shear_lift = shear_lift;
save([op.deplname,'_1sec_derived'],'-struct','data')


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 12  at  ',datestr(now)])
disp(' ')
diary off



function [ur,vr]=rotuv(alpha,u,v)
% function [ur,vr]=rotuv(alpha,u,v)
%
% rotate velocities by angle alpha (in radian)
%
% input  : alpha       - rotation angle in radian
%          u           - velocity component in x/east direction
%          v           - velocity component in y/north direction
%
% output : ur          - rotated u velocity
%          vr          - rotated v velocity
%
% version 1, last change 24.05.2013
%
%
% a positive alpha will rotate the COORDINATE-SYSTEM in 
% a mathematical negative sense (clockwise)
%
% i.e.  u=1,v=1 and alpha=45deg  gives  u=0,v=1.41
%
% If you have u,v velocities and you want to know the velocity component in
% direction of the angle ang in degrees (counted like a ship's heading in mathematically
% negative sense starting with 0=north), then enter
% velcomp = rotuv((270+ang)/180*pi,u,v)

%IFM Kiel Martin Visbeck
% added explanations, G.Krahmann, 24.05.2013

ur=cos(alpha).*u-sin(alpha).*v;
vr=sin(alpha).*u+cos(alpha).*v;


