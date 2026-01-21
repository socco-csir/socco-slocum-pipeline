function [] = step_13_optimize_dynamical_model()
% function [] = step_13_optimize_dynamical_model()
% 
% GEOMAR SVN $Id: step_13_optimize_dynamical_model.m 840 2021-10-12 15:15:23Z gkrahmann@geomar.de $
%
% optimize parameters that model the gliders dynamical flight behaviour
%
% requires :   processing_parameters.m in the current working folder
%
% version 6.2.0  last change 31.01.2023

% G.Krahmann, GEOMAR, Aug 2012

% code cleanup                                           GK, 01.03.2017  1-->2
% removed exp_filter and created func_exp_filter.m       GK, 22.03.2017  2-->3
% add processing diary, remove v6 mat comp               GK, 30.01.2018  3-->4
% extended diary                                         GK, 15.01.2019  4-->5
% recoded air in oil calculation                         GK, 17.01.2019  5-->6
% replaced linreg by polyfit                             GK, 16.09.2021  6-->6.1.0
% change from CSIRO seawater to TEOS library             GK, 31.01.2023  6.1.0-->6.2.0

% references: Merckelbach, Lucas, David Smeed, Gwyn Griffiths, 2010: Vertical Water 
% Velocities from Underwater Gliders. J. Atmos. Oceanic Technol., 27, 547–563. 
%
% The calculations in these module follow largely the ideas from the paper. In difference
% to the data they use we have to include the thermal expansion of the glider. The
% expansion coefficient is fitted as a fourth independent parameter and assumed not
% to change of the over the deployment. In difference to their approach, the temperature
% used to calculate the thermal expansion is not the ambient temperature measured by
% the CTD, but a EMA (exponential-moving-average) filtered version of the ambient
% temperature. This takes into account that the glider's hull has a significant
% thermal mass.
%
% Unfortunately the parameters are not orthogonal. They each have distinguishable influences
% on the flight of the glider, but if you change one of them the optimum for another one might
% change. So if there are any inconsistencies in the data for whatever reason and this has an
% influence on one of the parameters then others are likely different too. We try to mitigate this
% here by first optimizing all parameters, then fixing one of them, optimizing the rest, fixing
% the next, and so.
% The hope is that we are able to determine those parameters that are most likely to be really constant.
% I.e. compressibility and thermal expansion, and then again determine
% drag coefficient and glider mass which sometimes vary in time. For the variation in time of
% the drag coefficient and the mass we use a third order polynomial
% increase. In case we find that drag coefficient or glider mass do not vary significantly, we keep
% them constant like the other parameters.
% Order of optimizations
% 1 - drag & compr & mass & thermal
% 2 - drag & compr & ---- & thermal
% 3 - drag & ----- & ---- & thermal
% 4 - drag & ----- & ---- & -------
% determine constant drag coeff or third order increase
% 5 - ---- & compr & mass & thermal
% 6 - ---- & ----- & mass & thermal
% 7 - ---- & ----- & ---- & thermal
% keep thermal and tempfilter and drag
% 8  - ---- & compr & mass & -------
% 9  - ---- & ----- & mass & -------
% determine constant mass or linear increase
% 10 - drag & ----- & ---- & -------
% 11 - ---- & compr & ---- & -------

clear global
global val plotnew handle_model_graph w w_glider_model nnn facs fixed_params param_list no_plot
close all


%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 13  at  ',datestr(now)])
disp('optimize flight model parameters')
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
load([op.deplname,'_1sec.mat'],'best_pressure','main_datenum','ctd_temperature',...
  'ctd_pressure','latitude','pitch','slocum_ballast_volume',...
  'slocum_battery_position','fin','pressure','thruster_power');
load([op.deplname,'_1sec_derived.mat'],'ctd_salinity_dpdt','shear_lift');


%
% There is a one-sided exponential filter on the temperature used with the thermal
% expansion coefficient.  Optimizations have found to be mostly between 0.002 and 0.02
% with those for microrider gliders a bit lower
% So we do not use the optimization of this parameter anymore.
%
if op.microrider_installed~=0 | op.suna_installed~=0
  fixed_params.temperature_filter = 0.006;
  temperature_filter = 0.006;
else
  fixed_params.temperature_filter = 0.008;
  temperature_filter = 0.008;
end


%
% calculate some other data
for n=1:length(pressure)
  w{n} = -gradient(best_pressure{n})./gradient(main_datenum{n}*86400);
  w{n} = meanfilt(w{n},2);
  rho{n} = gsw_dens(ctd_salinity_dpdt{n},ctd_temperature{n},ctd_pressure{n});
  dummy = ctd_temperature{n};
  ind = find(~isnan(dummy));
  if ~isempty(ind)
    dummy(ind) = func_exp_filter(dummy(ind),fixed_params.temperature_filter);
  end
  filtered_ctd_temperature{n} = dummy;
end


data.g = gsw_grav(nmean([latitude{:}]));
data.S = op.effective_wing_area;
data.A_h = op.frontal_area;
data.V_g = op.glider_volume;
data.pitch = pitch;
data.w = w;
data.rho = rho;
data.pressure = best_pressure;
data.filtered_ctd_temperature = filtered_ctd_temperature;
data.ctd_salinity = ctd_salinity_dpdt;
data.delta_V_bp = slocum_ballast_volume;
data.battery_position = slocum_battery_position;
data.time = main_datenum;
data.latitude = latitude;
data.fin = fin;
data.thruster_power = thruster_power;
data.shear_lift = shear_lift;

fixed_params.C_D_0= [];
fixed_params.epsilon = [];
fixed_params.m_g = [];
fixed_params.alpha_T = [];
para_names = {'Temp filter','Drag coeff','Compressibility','Mass','Thermal exp coeff'};

%
% Apply an air in the oil correction, if op.air_V differs from 0.
% The extra buoyancy volume is compressed with ocean pressure so that
% it is really only effective during ascent in the last few meters.
%
if op.air_V>0
  min_V_bp = nmin([slocum_ballast_volume{:}]);
  max_V_bp = nmax([slocum_ballast_volume{:}]);
  air_V = op.air_V;
  for n=1:length(data.delta_V_bp)
    uncompressed_oil_and_air_outside = data.delta_V_bp{n} - min_V_bp;
    air_of_total_factor = op.air_V / (max_V_bp - min_V_bp);
    oil_of_total_factor = (max_V_bp - min_V_bp - op.air_V) / (max_V_bp - min_V_bp);
    uncompressed_oil_outside = uncompressed_oil_and_air_outside * oil_of_total_factor;
    uncompressed_air_outside = uncompressed_oil_and_air_outside * air_of_total_factor;
    compression_factor = 10./(data.pressure{n}+10);
    compressed_air_outside = uncompressed_air_outside .* compression_factor;
    compressed_oil_outside = uncompressed_oil_outside;
    compressed_oil_and_air_outside = compressed_oil_outside + compressed_air_outside;
    data.delta_V_bp{n} = compressed_oil_and_air_outside + min_V_bp;
  end
end

%
% Optimization #1
%
% fit all four parameters
% drag, thermal exp, compressibility, mass
%
res = func_optimize_dynamical_model(data,op);
if isempty(res)
  fixed_params = [];
  save([op.deplname,'_dynamics.mat'],'fixed_params')
  return
end
if op.no_plot~=1
  figure(2)
  clf
end
fn = fieldnames(fixed_params);
count = 1;
for n=1:length(fn)
  if isempty(getfield(fixed_params,fn{n}))
    if op.no_plot~=1
      subplot(2,2,count)
      hist(res(:,count),20)
      hold on
      ax = axis;
      plot(nmedian(res(:,count))*[1,1],ax(3:4),'r')
      title(para_names{n})
    end
    disp([fn{n},' : ',num2str(nmedian(res(:,count)))])
    eval([fn{n},'=nmedian(res(:,count));'])
    count = count+1;
  end
end


%
% Optimization #2
%
% take median mass and use that as a fixed parameter
% then fit the remaining 3 parameters
% drag, thermal exp, compressibility
%
fixed_params.m_g = m_g;
res = func_optimize_dynamical_model(data,op);
if op.no_plot~=1
  figure(3)
  clf
end
fn = fieldnames(fixed_params);
count = 1;
for n=1:length(fn)
  if isempty(getfield(fixed_params,fn{n}))
    if op.no_plot~=1
      subplot(2,2,count)
      hist(res(:,count),20)
      hold on
      ax = axis;
      plot(nmedian(res(:,count))*[1,1],ax(3:4),'r')
      title(para_names{n})
    end
    disp([fn{n},' : ',num2str(nmedian(res(:,count)))])
    eval([fn{n},'=nmedian(res(:,count));'])
    count = count+1;
  end
end


%
% Optimization #3
%
% take median compressibility and use that as a fixed parameter
% then fit the remaining 2 parameters
% drag, thermal exp
%
fixed_params.epsilon = epsilon;
[res,yo_number] = func_optimize_dynamical_model(data,op);
if op.no_plot~=1
  figure(4)
  clf
end
fn = fieldnames(fixed_params);
count = 1;
for n=1:length(fn)
  if isempty(getfield(fixed_params,fn{n}))
    if op.no_plot~=1
      subplot(2,2,count)
      hist(res(:,count),20)
      hold on
      ax = axis;
      plot(nmedian(res(:,count))*[1,1],ax(3:4),'r')
      title(para_names{n})
    end
    disp([fn{n},' : ',num2str(nmedian(res(:,count)))])
    eval([fn{n},'=nmedian(res(:,count));'])
    count = count+1;
  end
end


%
% Optimization #4
%
% take median thermal expansion coefficient and use that as a fixed parameter
% then fit the remaining 1 parameter
% drag
%
fixed_params.alpha_T = alpha_T;
res = func_optimize_dynamical_model(data,op);
if op.no_plot~=1
  figure(5)
  clf
end
fn = fieldnames(fixed_params);
count = 1;
for n=1:length(fn)
  if isempty(getfield(fixed_params,fn{n}))
    if op.no_plot~=1
      subplot(2,2,count)
      hist(res(:,count),20)
      hold on
      ax = axis;
      plot(nmedian(res(:,count))*[1,1],ax(3:4),'r')
      title(para_names{n})
    end
    disp([fn{n},' : ',num2str(nmedian(res(:,count)))])
    eval([fn{n},'=nmedian(res(:,count));'])
    count = count+1;
  end
end


%
% The last parameter that was fitted, was the pitch angle independent part of the drag coefficient.
% In cases of strong biofouling this parameter can increase exponentially.
% Here we first determine whether the drag coefficient increased significantly over the deployment.
% We assume that the drag coefficient can not decrease significantly. In case we find a
% significant increase we try a third order fit to predict the drag coefficient.
%
% We determine a significant increase not through statistics but simply by looking whether
% a linear fit over the whole deployment predicts an increase by more than 0.005. The
% typical value for the drag coefficient is 0.11.
%
fixed_params.C_D_0 = C_D_0;
if length(res)>3                    
 % [a,b]=linreg(res,[1:length(res)]');
  a = polyfit([1:length(res)]',res,1);
  if a(2)*length(res)>0.005
    if isfield(op,'limit_drag_fit')
      if op.limit_drag_fit==1
        coeff = polyfit(yo_number,log(res)',1);
      else
        coeff = polyfit(yo_number,log(res)',3);
      end
    else
      coeff = polyfit(yo_number,log(res)',3);
    end
    if op.no_plot~=1
      figure(1)
      clf
      plot(yo_number,res,'x')
      hold on
      plot([1:yo_number(end)],exp(polyval(coeff,[1:yo_number(end)])),'r')
      title('Drag coefficient over time')
    end
    fixed_params.C_D_0 = exp(polyval(coeff,[1:length(data.pressure)]));
    disp('Paused for 30 seconds. Drag coefficient fit does not match well.')
    pause(30)
  end
end


%
% Optimization #5
%
% redo the optimization of 3 parameters with the drag coefficient fixed
% thermal expansion, compressibility, mass
%
fixed_params.epsilon = [];
fixed_params.m_g = [];
fixed_params.alpha_T = [];
res = func_optimize_dynamical_model(data,op);
if op.no_plot~=1
  figure(6)
  clf
end
fn = fieldnames(fixed_params);
count = 1;
for n=1:length(fn)
  if isempty(getfield(fixed_params,fn{n}))
    if op.no_plot~=1
      subplot(2,2,count)
      hist(res(:,count),20)
      hold on
      ax = axis;
      plot(nmedian(res(:,count))*[1,1],ax(3:4),'r')
      title(para_names{n})
    end
    disp([fn{n},' : ',num2str(nmedian(res(:,count)))])
    eval([fn{n},'=nmedian(res(:,count));'])
    count = count+1;
  end
end


% 
% Optimization #6
%
% take median compressibility and use that as second fixed parameter
% then fit the remaining 2 parameters
% thermal expansion, mass
%
fixed_params.epsilon = epsilon;
res = func_optimize_dynamical_model(data,op);
if op.no_plot~=1
  figure(7)
  clf
end
fn = fieldnames(fixed_params);
count = 1;
for n=1:length(fn)
  if isempty(getfield(fixed_params,fn{n}))
    if op.no_plot~=1
      subplot(2,2,count)
      hist(res(:,count),20)
      hold on
      ax = axis;
      plot(nmedian(res(:,count))*[1,1],ax(3:4),'r')
      title(para_names{n})
    end
    disp([fn{n},' : ',num2str(nmedian(res(:,count)))])
    eval([fn{n},'=nmedian(res(:,count));'])
    count = count+1;
  end
end


%
% Optimization #7
%
% take median mass and use that as third fixed parameter
% then fit the remaining 1 parameter
% thermal expansion
%
% This is the final fit the thermal expansion coefficient.
%
fixed_params.m_g = m_g;
res = func_optimize_dynamical_model(data,op);
if op.no_plot~=1
  figure(8)
  clf
end
fn = fieldnames(fixed_params);
count = 1;
for n=1:length(fn)
  if isempty(getfield(fixed_params,fn{n}))
    if op.no_plot~=1
      subplot(2,2,count)
      hist(res(:,count),20)
      hold on
      ax = axis;
      plot(nmedian(res(:,count))*[1,1],ax(3:4),'r')
      title(para_names{n})
    end
    disp([fn{n},' : ',num2str(nmedian(res(:,count)))])
    eval([fn{n},'=nmedian(res(:,count));'])
    eval([fn{n},'_std=nstd(res(:,count));'])
    count = count+1;
  end
end
fixed_params.alpha_T = alpha_T;
fixed_params.alpha_T_std = alpha_T_std;


%
% Optimization #8
%
% Keep thermal expansion and drag coefficient and
% then fit the remaining 2 parameters
% compressibility, mass
%
fixed_params.epsilon = [];
fixed_params.m_g = [];
res = func_optimize_dynamical_model(data,op);
if op.no_plot~=1
  figure(9)
  clf
end
fn = fieldnames(fixed_params);
count = 1;
for n=1:length(fn)
  if isempty(getfield(fixed_params,fn{n}))
    if op.no_plot~=1
      subplot(2,2,count)
      hist(res(:,count),20)
      hold on
      ax = axis;
      plot(nmedian(res(:,count))*[1,1],ax(3:4),'r')
      title(para_names{n})
    end
    disp([fn{n},' : ',num2str(nmedian(res(:,count)))])
    eval([fn{n},'=nmedian(res(:,count));'])
    count = count+1;
  end
end


%
% Optimization #9
%
% take median compressibility and use that as fixed parameter
% then fit the remaining 1 parameter
% mass
%
% This is the final fit for mass.
%
fixed_params.epsilon = epsilon;
res = func_optimize_dynamical_model(data,op);
if op.no_plot~=1
  figure(10)
  clf
end
fn = fieldnames(fixed_params);
count = 1;
for n=1:length(fn)
  if isempty(getfield(fixed_params,fn{n}))
    if op.no_plot~=1
      subplot(2,2,count)
      hist(res(:,count),20)
      hold on
      ax = axis;
      plot(nmedian(res(:,count))*[1,1],ax(3:4),'r')
      title(para_names{n})
    end
    disp([fn{n},' : ',num2str(nmedian(res(:,count)))])
    eval([fn{n},'=nmedian(res(:,count));'])
    eval([fn{n},'_std=nstd(res(:,count));'])
    count = count+1;
  end
end
fixed_params.m_g = m_g;
fixed_params.std_m_g = m_g_std;
%
% check whether there is a significant increase (more than 20g over the deployment) in glider mass
%
if length(res)>3
  res = res(:);
%  [a,b]=linreg(res,[1:length(res)]');
  a = polyfit([1:length(res)]',res,1);
  if a(2)*length(res)>0.020
    coeff = polyfit(yo_number,res',1);
    if op.no_plot~=1
      figure(1)
      clf
      plot(yo_number,res,'x')
      hold on
      plot([1:yo_number(end)],polyval(coeff,[1:yo_number(end)]),'r')
      title('Mass over time')
    end
    fixed_params.m_g = polyval(coeff,[1:length(data.pressure)]);
    fixed_params.m_g_std2 = nstd(res - polyval(coeff,yo_number)');
    disp('Paused for 30 seconds. Weight fit does not match well.')
    pause(30)
  end
end


%
% Optimization #10
% 
% fit drag coefficient only
%
% This is the final fit for the drag coefficient.
%
fixed_params.C_D_0 = [];
res = func_optimize_dynamical_model(data,op);
if op.no_plot~=1
  figure(11)
  clf
end
fn = fieldnames(fixed_params);
count = 1;
for n=1:length(fn)
  if isempty(getfield(fixed_params,fn{n}))
    if op.no_plot~=1
      subplot(2,2,count)
      hist(res(:,count),20)
      hold on
      ax = axis;
      plot(nmedian(res(:,count))*[1,1],ax(3:4),'r')
      title(para_names{n})
    end
    disp([fn{n},' : ',num2str(nmedian(res(:,count)))])
    eval([fn{n},'=nmedian(res(:,count));'])
    eval([fn{n},'=nmedian(res(:,count));'])
    eval([fn{n},'_std=nstd(res(:,count));'])
    count = count+1;
  end
end
fixed_params.C_D_0 = C_D_0;
fixed_params.std_C_D_0 = C_D_0_std;
if length(res)>3
  res = res(:);
%  [a,b]=linreg(res,[1:length(res)]');
  a = polyfit([1:length(res)]',res,1);
  if a(2)*length(res)>0.005
    if isfield(op,'limit_drag_fit')
      if op.limit_drag_fit==1
        coeff = polyfit(yo_number,log(res)',1);
      else
        coeff = polyfit(yo_number,log(res)',3);
      end
    else
      coeff = polyfit(yo_number,log(res)',3);
    end
    if op.no_plot~=1
      figure(1)
      clf
      plot(yo_number,res,'x')
      hold on
      plot([1:yo_number(end)],exp(polyval(coeff,[1:yo_number(end)])),'r')
      title('Drag coefficient over time')
    end
    fixed_params.C_D_0 = exp(polyval(coeff,[1:length(pressure)]));
    fixed_params.C_D_0_std2 = nstd(res - exp(polyval(coeff,yo_number))');
    disp('Paused for 30 seconds. Drag coefficient fit does not match well.')
    pause(30)
  end
end
title('Drag coefficient')


%
% Optimization #11
%
% fit compressibility only
%
% This is the final fit for the compressibility.
%
fixed_params.epsilon = [];
res = func_optimize_dynamical_model(data,op);
if op.no_plot~=1
  figure(12)
  clf
end
fn = fieldnames(fixed_params);
count = 1;
for n=1:length(fn)
  if isempty(getfield(fixed_params,fn{n}))
    if op.no_plot~=1
      subplot(2,2,count)
      hist(res(:,count),20)
      hold on
      ax = axis;
      plot(nmedian(res(:,count))*[1,1],ax(3:4),'r')
      title(para_names{n})
    end
    disp([fn{n},' : ',num2str(nmedian(res(:,count)))])
    eval([fn{n},'=nmedian(res(:,count));'])
    eval([fn{n},'_std=nstd(res(:,count));'])
    count = count+1;
  end
end
title('Compressibility')
fixed_params.epsilon = epsilon;
fixed_params.epsilon_std = epsilon_std;

disp(' ')
disp('Optimal parameter set:')
disp(['Drag coeff : ',num2str(fixed_params.C_D_0(1)),' - ',num2str(fixed_params.C_D_0(end))])
disp(['epsilon    : ',num2str(fixed_params.epsilon(1))])
disp(['alpha_T    : ',num2str(fixed_params.alpha_T(1))])
disp(['mass       : ',num2str(fixed_params.m_g(1)),' - ',num2str(fixed_params.m_g(end)) ])

save([op.deplname,'_dynamics.mat'],'fixed_params')


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 13  at  ',datestr(now)])
disp(' ')
diary off

