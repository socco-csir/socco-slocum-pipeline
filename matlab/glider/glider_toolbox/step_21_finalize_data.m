function [] = step_21_finalize_data()
% function [] = step_21_finalize_data()
% 
% GEOMAR SVN $Id: step_21_finalize_data.m 760 2020-12-10 15:17:57Z gkrahmann@geomar.de $
%
% apply the offsets determined in the previous step and final mat files
%
% version 6.1.0  last change 22.07.2020

% G.Krahmann, GEOMAR,  Mar 2014

% added magnetic deviation correction and processing diary 
% removed v6 mat file comp                                              GK, 30.01.2018  1-->2
% store best pressure as pressure data                                  GK, 06.03.2018  2-->3
% save nan*oxygen when no oxygen available                              GK, 28.03.2018  3-->4
% extended diary                                                        GK, 15.01.2019  4-->5
% added oa and og factors                                               GK, 13.03.2020  5-->6.0.0
% convert s_offset into a c_factor and apply that as correction         GK, 22.07.2020  6.0.0-->6.1.0




%
% look for and load processing parameters
%
if exist('./processing_parameters.m') | exist('.\processing_parameters.m')
  processing_parameters
else
  error(['could not find file ''processing_parameters.m'' in current folder ',pwd])
end


%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 21  at  ',datestr(now)])
disp('save final mat files')
disp(['t-off  = ',num2str(op.t_offset)])
disp(['s-off  = ',num2str(op.s_offset)])
disp(['use_c_factor  = ',num2str(op.use_c_factor)])
disp(['oa-off  = ',num2str(op.oa_offset)])
disp(['oa-fac  = ',num2str(op.oa_factor)])
disp(['og-off  = ',num2str(op.og_offset)])
disp(['og-fac  = ',num2str(op.og_factor)])
disp(['og-pfac = ',num2str(op.og_pfactor)])
diary off


%
% load gridded data
%
load([op.deplname,'_gridded'])
load([op.deplname,'_1sec_derived'],'oxygen_cal')


%
% store data in one structure
% and apply offsets to T,S,O
%
final_gridded.pressure = pressure;
final_gridded.latitude = latitude;
final_gridded.longitude = longitude;
final_gridded.u = m_water_vx;
final_gridded.v = m_water_vy;
final_gridded.time_datenum = main_datenum;
final_gridded.turbidity = turbidity;
final_gridded.chlorophyll = chlorophyll;
if exist('rhodamin')
  final_gridded.rhodamin = rhodamin;
end
if exist('cdom')
  final_gridded.cdom = cdom;
end
if isfield(op,'t_offset')
  final_gridded.temperature = ctd_temperature - op.t_offset;
else
  final_gridded.temperature = ctd_temperature;
end
if isfield(op,'s_offset')
  if op.use_c_factor==0
    final_gridded.salinity = good_salinity - op.s_offset;
  else
    load cond_factor
    dummy_cond = sw_cndr(good_salinity,ctd_temperature,ctd_pressure);
    final_gridded.salinity = sw_salt(dummy_cond*cond_factor,ctd_temperature,ctd_pressure);
  end
else
  final_gridded.salinity = good_salinity;
end
if strcmp(oxygen_cal,'geomar')
  if isfield(op,'og_offset')
    final_gridded.oxygen = good_oxygen * op.og_factor + ctd_pressure * op.og_pfactor - op.og_offset;
  else
    final_gridded.oxygen = good_oxygen;
  end
elseif strcmp(oxygen_cal,'aanderaa') | strcmp(oxygen_cal,'none')
  if isfield(op,'oa_offset')
    final_gridded.oxygen = good_oxygen * op.oa_factor - op.oa_offset;
  else
    final_gridded.oxygen = good_oxygen;
  end
else
  final_gridded.oxygen = nan*final_gridded.salinity;
end


% ISABELLE EDIT

% apply magnetic deviation

[yy,mm,dd] = datevec( nmean( main_datenum ) );
yy = nmean( yy + (mm-0.5)/12 );
ind = find(~isnan(m_water_vx));
magnetic_deviation = nan*m_water_vx;
lat = nmean(latitude);
lon = nmean(longitude);
magnetic_deviation(ind) = func_magdev(lat(ind),lon(ind),0,yy);
[final_gridded.u,final_gridded.v] = func_rotate_uv(final_gridded.u,final_gridded.v,magnetic_deviation);


%
% save this data
%
save([op.deplname,'_final_gridded'],'-struct','final_gridded')

%
% load 1sec data
%
load([op.deplname,'_1sec'])
load([op.deplname,'_1sec_derived'])


%
% store data in one structure
% and apply offsets to T,S,O
%
final_1sec.pressure = best_pressure;
final_1sec.latitude = latitude;
final_1sec.longitude = longitude;
final_1sec.u = m_water_vx;
final_1sec.v = m_water_vy;
final_1sec.time_datenum = main_datenum;
final_1sec.turbidity = turbidity;
final_1sec.chlorophyll = chlorophyll;
if exist('cdom')
  final_1sec.cdom = cdom;
end
if exist('rhodamin')
  final_1sec.rhodamin = rhodamin;
end
if isfield(op,'t_offset')
  for n=1:length(ctd_temperature)
    final_1sec.temperature{n} = ctd_temperature{n} - op.t_offset;
  end
else
  final_1sec.temperature = ctd_temperature;
end
if isfield(op,'s_offset')
  if op.use_c_factor==0
    for n=1:length(good_salinity)
      final_1sec.salinity{n} = good_salinity{n} - op.s_offset;
      good_salinity{n} = good_salinity{n} - op.s_offset;
    end
  else
    load cond_factor
    for n=1:length(good_salinity)
      dummy_cond1 = sw_cndr(good_salinity{n},ctd_temperature{n},ctd_pressure{n});
      final_1sec.salinity{n} = sw_salt(dummy_cond1*cond_factor,ctd_temperature{n},ctd_pressure{n});
    end
    good_salinity = final_1sec.salinity;
  end
else
  final_1sec.salinity = good_salinity;
end
if strcmp(oxygen_cal,'geomar')
  if isfield(op,'og_offset')
    for n=1:length(good_oxygen)
      final_1sec.oxygen{n} = good_oxygen{n} - op.og_offset;
    end
  else
    final_1sec.oxygen = good_oxygen;
  end
elseif strcmp(oxygen_cal,'aanderaa') | strcmp(oxygen_cal,'aanderaa')
  if isfield(op,'oa_offset')
    for n=1:length(good_oxygen)
      final_1sec.oxygen{n} = good_oxygen{n} - op.oa_offset;
    end
  else
    final_1sec.oxygen = good_oxygen;
  end
else
  final_1sec.oxygen = final_1sec.salinity;
  for n=1:length(final_1sec.oxygen)
    final_1sec.oxygen{n} = nan*final_1sec.oxygen{n};
  end
end


% ISABELLE EDIT
% apply magnetic deviation

[yy,mm,dd] = datevec( nmean( [main_datenum{:}] ) );
yy = nmean( yy + (mm-0.5)/12 );
for n=1:length(latitude)
  lat(n) = nmean(latitude{n});
  lon(n) = nmean(longitude{n});
  magnetic_deviation(n) = func_magdev(lat(n),lon(n),0,yy);
  [final_1sec.u(n),final_1sec.v(n)] = func_rotate_uv(final_1sec.u(n),final_1sec.v(n),magnetic_deviation(n));
end



%
% save this data
%
save([op.deplname,'_final_1sec'],'-struct','final_1sec')






%
% apply additional oxygen drift correction
%
if exist('get_o_drift.m')
  disp('get_o_drift')
  get_o_drift
end



%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 21  at  ',datestr(now)])
disp(' ')
diary off

