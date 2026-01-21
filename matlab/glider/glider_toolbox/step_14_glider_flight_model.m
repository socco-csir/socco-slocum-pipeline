function [] = step_14_glider_flight_model()
% function [] = step_14_glider_flight_model()
% 
% GEOMAR SVN $Id: step_14_glider_flight_model.m 557 2019-01-25 10:53:26Z gkrahmann@geomar.de $
%
% apply the optimized parameters that model the gliders dynamical flight behaviour
%
% requires :   processing_parameters.m in the current working folder
%
% version 7.1.0  last change 31.01.2023

% G.Krahmann, GEOMAR, Apr 2013

% cleanup                                                     GK, 01.03.2017  1-->2
% use func_exp_filter instead of exp_filter                   GK, 22.03.2017  2-->3
% add processing diary, remove v6 mat comp                    GK, 30.01.2018  3-->4
% add first_dive_bubbles                                      GK, 03.01.2019  4-->5
% extended diary                                              GK, 15.01.2019  5-->6
% recoded air in oil calculation                              GK, 17.01.2019  6-->7
% change from CSIRO seawater to TEOS library                  GK, 31.01.2023  7-->7.1.0

% references: Merckelbach, Lucas, David Smeed, Gwyn Griffiths, 2010: Vertical Water 
% Velocities from Underwater Gliders. J. Atmos. Oceanic Technol., 27, 547–563. 
%
% In step 13 several parameters have been optimized to model the steady state
% glider flight behaviour.
% In this step these parameters are used to integrate the full flight of the glider
% using the forces that act on the glider. Any mismatches are likely caused by vertical
% water movements. In these cases of the modeled depths will not match the observed
% ones. As the flight duration must be ok we scale the modeled depth of the
% glider to match the observed one at the beginning and the end of the yo.


%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 14  at  ',datestr(now)])
disp('integrate flight model')
diary off



%
% look for and load processing parameters which contain the location of
% the flash card copy of the glider
%
if exist('./processing_parameters.m')
  processing_parameters
else
  error(['could not find file ''processing_parameters.m'' in current folder ',pwd])
end

%
% load data
%
load([op.deplname,'_1sec.mat']);
load([op.deplname,'_1sec_derived.mat']);
load([op.deplname,'_dynamics.mat']);
load([op.deplname,'_yos.mat'],'dos_id');


%
% check whether there has no flight model been optimized
% In that case, we fill the variables for the modeled glider flight with the
% ones estimated from the simple dp/dt calculation.
%
if isempty(fixed_params)
  data = load([op.deplname,'_1sec_derived.mat']);
%  data.angle_of_attack = aoa;
  data.glider_speed_in_movement_direction_model = glider_speed_in_movement_direction_dpdt;

%  data.modeled_position_x = x;
%  data.modeled_position_z = z;
%  data.modeled_w = modeled_w;
%  data.local_shear = local_shear;
  for n=1:length(data.glider_speed_in_movement_direction_model)
   data.glider_speed_in_cell_direction_model{n} = data.glider_speed_in_cell_direction_dpdt{n};
  end
%  data.aoa_modeled = aoa;
  save([op.deplname,'_1sec_derived.mat'],'-struct','data')
  return
end
if 0
if exist([op.deplname,'_currentprofile.mat'])
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
for n=1:length(best_pressure)
  u_shear1 = nan*best_pressure{n};
  v_shear1 = nan*best_pressure{n};
  [dummy,time_ind] = nmin(abs(time_prof-nmean(main_datenum{n})));
  good = find(~isnan(best_pressure{n}));
  u_shear1(good) = interp1(z_prof(:,time_ind),u_shear(:,time_ind),best_pressure{n}(good));
  v_shear1(good) = interp1(z_prof(:,time_ind),v_shear(:,time_ind),best_pressure{n}(good));
  local_shear{n} = rotuv((270+heading{n})/180*pi,u_shear1,v_shear1);
  bad = find(isnan(local_shear{n}));
  if ~isempty(bad)
    local_shear{n}(bad) = 0;
  end
end
clear good
end

for n=1:length(best_pressure)
  w{n} = -gradient(best_pressure{n})./gradient(main_datenum{n}*86400);
  rho{n} = gsw_dens(ctd_salinity_dpdt{n},ctd_temperature{n},ctd_pressure{n});
  good{n} = find(~isnan(ctd_temperature{n}+pitch{n}));
  is_in_surface_mode{n} = m_air_pump{n}.*(best_pressure{n}<5);
end

max_V_bp = nmax([slocum_ballast_volume{:}]);
min_V_bp = nmin([slocum_ballast_volume{:}]);

for n=1:length(best_pressure)
  fprintf(1,'.');
  if length(fixed_params.C_D_0)>1
    if n>length(fixed_params.C_D_0)
      params.C_D_0 = fixed_params.C_D_0(end);
      disp('using last C_D_0')
    else
      params.C_D_0 = fixed_params.C_D_0(n);
    end
  else
    params.C_D_0 = fixed_params.C_D_0;
  end
  if length(fixed_params.m_g)>1
    params.m_g = fixed_params.m_g(n);
  else
    params.m_g = fixed_params.m_g;
  end
  params.epsilon = fixed_params.epsilon;
  params.alpha_T = fixed_params.alpha_T;
  params.temperature_filter = fixed_params.temperature_filter;
  params.latitude = nmean(latitude{n});
  params.V_g = op.glider_volume;
  params.A_h = op.frontal_area;
  if length(good{n})>1
    filtered_ctd_temperature = func_exp_filter(ctd_temperature{n}(good{n}),params.temperature_filter);
    if op.air_V>0
v1=slocum_ballast_volume{n};
         uncompressed_oil_and_air_outside = slocum_ballast_volume{n} - min_V_bp;
         air_of_total_factor = op.air_V / (max_V_bp - min_V_bp);
         oil_of_total_factor = (max_V_bp - min_V_bp - op.air_V) / (max_V_bp - min_V_bp);
         uncompressed_oil_outside = uncompressed_oil_and_air_outside * oil_of_total_factor;
         uncompressed_air_outside = uncompressed_oil_and_air_outside * air_of_total_factor;
         compression_factor = 10./(best_pressure{n}+10); 
         compressed_air_outside = uncompressed_air_outside .* compression_factor;
         compressed_oil_outside = uncompressed_oil_outside;
         compressed_oil_and_air_outside = compressed_oil_outside + compressed_air_outside;
         slocum_ballast_volume{n} = compressed_oil_and_air_outside + min_V_bp;
if 0
figure(7)
clf
plot(v1)
hold on
plot(slocum_ballast_volume{n},'r')
end
    end
    bubbles = repmat(0,size(best_pressure{n}));
    is_first_dive_after_call = 0;
    if n==1
      is_first_dive_after_call = 1;
    elseif dos_id{n}(1) ~= dos_id{n-1}(1)
      is_first_dive_after_call = 1;
    end

    if is_first_dive_after_call==1
      ind_start_dive = find(best_pressure{n}>op.start_dive_parameter(1));
      surface_p = nmin(best_pressure{n}(1:ind_start_dive));
      ind_start_dive = find( best_pressure{n}(1:ind_start_dive) - surface_p < op.start_dive_parameter(2) );
      ind_start_dive = ind_start_dive(end) - op.start_dive_parameter(3);
      if ind_start_dive<=0
        ind_start_dive = 1;
      end
      bubbles(1:ind_start_dive) = 200;
      bubbles(ind_start_dive+[0:10]) = linspace(200,0,11);
      slocum_ballast_volume{n} = slocum_ballast_volume{n} + ...
        bubbles*10./(best_pressure{n}+10);
    end

if 0
    if op.first_dive_bubbles>0
      if is_first_dive_after_call==1
         bubbles = repmat(op.first_dive_bubbles,size(best_pressure{n}));
         minus_pitch = meanfilt(pitch{n}(good{n}),10);
         minus_pitch = nans(minus_pitch,0,0,'>');
         minus_pitch = nans(minus_pitch,0,nan,'==');
         minus_pitch = cumsum(minus_pitch);
%         ind = find(meanfilt(pitch{n},10)<op.bubble_leak_angle & m_air_pump{n}==0);
%         ind = find(meanfilt(pitch{n},20)<op.bubble_leak_angle );
         ind = find(minus_pitch<op.sum_of_negative_angles_to_release_bubble );
         ind = find(minus_pitch<op.sum_of_negative_angles_to_release_bubble & m_air_pump{n}(good{n})==0);
         bubbles(good{n}(ind(1)):end) = 0;
         slocum_ballast_volume{n} = slocum_ballast_volume{n} + ...
           bubbles*10./(best_pressure{n}+10);
      end
    end
end

    if is_first_dive_after_call==1
      start_vel = 0;
      start_glide_angle = -90/180*pi;
      start_aoa = 90/180*pi;
      start_u = 0;
      start_w = 0;
    else
      start_vel = vel1(end);
      start_glide_angle = glide_angle(end);
      start_aoa = aoa1(end);
      start_u = u(end);
      start_w = w(end);
    end
    [x1,z1,aoa1,vel1,scaling1,devi1,z_of_p1,glide_angle,u,w] =...
      func_integrate_glider_acceleration(best_pressure{n}(good{n}),...
      filtered_ctd_temperature,...
      pitch{n}(good{n}),rho{n}(good{n}),slocum_ballast_volume{n}(good{n}),params,...
      is_in_surface_mode{n}(good{n}),slocum_battery_position{n}(good{n}),fin{n}(good{n}),...
      m_air_pump{n}(good{n}),internal_pressure{n}(good{n}),...
      shear_lift{n}(good{n}),op.no_plot,n,op,bubbles(good{n}),start_vel,start_glide_angle,...
      start_aoa,start_u,start_w);

if 0
clf
subplot(2,1,1)
plot(diff(z1))
hold on
plot(diff(z2),'r')
plot(-diff(best_pressure{n}(good{n})),'g')
subplot(2,1,2)
plot(ctd_temperature{n}(good{n}))
hold on
plot(filtered_ctd_temperature,'r')
pause
end

    x{n} = nan*rho{n};
    x{n}(good{n}) = x1;
    z{n} = nan*rho{n};
    z{n}(good{n}) = z1;
    z_of_p{n} = nan*rho{n};
    z_of_p{n}(good{n}) = z_of_p1;
    aoa{n} = nan*rho{n};
    aoa{n}(good{n}) = aoa1;
    glider_speed_in_movement_direction_model{n} = nan*rho{n};
    glider_speed_in_movement_direction_model{n}(good{n}) = vel1;
    scaling(n) = scaling1; 
    devi{n} = devi1;

    modeled_w{n} = gradient(z{n});
    measured_w{n} = gradient(z_of_p{n});

%    speedFactor = polyval(speedFactorPols(:,selectedDegree+1), vel{n});
%    flow_speed = speedFactor .* vel{n}; % Avoid division by zero

%    flow_speed = nans(flow_speed,op.minimum_flow_speed,nan,'==');
%    flow_speed = nans(flow_speed,op.minimum_flow_speed,op.minimum_flow_speed,'<');

%    flow_speed2{n} = flow_speed;
  else
    x{n} = nan*rho{n};
    z{n} = nan*rho{n};
    modeled_w{n} = nan*rho{n};
    measured_w{n} = nan*rho{n};
    aoa{n} = nan*rho{n};
    glider_speed_in_movement_direction_model{n} = nan*rho{n};
%    flow_speed2{n} = nan*rho{n};
    scaling(n) = nan;
    devi{n} = nan;
  end

end
fprintf(1,'\n');

good = find(~isnan(scaling));
good_down = find(~isnan(scaling(1:2:end)));
good_up = find(~isnan(scaling(2:2:end)));
if op.no_plot~=1
  figure
  hist(scaling(good))
  title('Scaling coefficients  observed depths / modeled depths')
end
disp(['mean scaling all   : ',num2str(nmean(scaling))])
disp(['mean scaling down  : ',num2str(nmean(scaling(1:2:end)))])
disp(['mean scaling up    : ',num2str(nmean(scaling(2:2:end)))])

data = load([op.deplname,'_1sec_derived.mat']);
data.angle_of_attack_model = aoa;
data.glider_speed_in_movement_direction_model = glider_speed_in_movement_direction_model;

data.modeled_position_x = x;
data.modeled_position_z = z;
data.modeled_w = modeled_w;
data.measured_w = measured_w;
data.shear_lift = shear_lift;
data.glider_speed_in_movement_direction_model = glider_speed_in_movement_direction_model;
for n=1:length(aoa)
 glider_speed_in_cell_direction_model{n} = glider_speed_in_movement_direction_model{n}.*cos(aoa{n});
 aoa{n} = aoa{n}*180/pi;
end
data.glider_speed_in_cell_direction_model = glider_speed_in_cell_direction_model;
save([op.deplname,'_1sec_derived.mat'],'-struct','data')

  
%
% determine final deviation sum
%
dev_sum = 0;
for n=1:length(modeled_w)
  dev_sum = dev_sum + nsum( abs(modeled_w{n} - measured_w{n}) );
end
disp(['Sum of absolute vertical speed differences between model and observation:'])
disp(dev_sum)


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 14  at  ',datestr(now)])
disp(' ')
diary off

