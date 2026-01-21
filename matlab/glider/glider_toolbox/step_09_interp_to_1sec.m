function [] = step_09_interp_to_1sec()
% function [] = step_09_interp_to_1sec()
% 
% GEOMAR SVN $Id: step_09_interp_to_1sec.m 760 2020-12-10 15:17:57Z gkrahmann@geomar.de $
%
% interpolate all data to 1 second steps
% use the best available time information for the respective data (i.e. time
% stamps for science data)
%
% version 11  last change 23.01.2019

% G.Krahmann, GEOMAR  Aug 2012

% properly handle variable pitch                        GK, 19.02.2013  2-->3
% add thruster_power                                    GK, 07.07.2014  3-->4
% found a problem with bogus oxy4 data                  GK, 02.03.2017  4-->5
% previous problem was not fixed                        GK, 15.03.2017  5-->6
% add processing diary, remove v6 mat comp              GK, 30.01.2018  6-->7
% add another CHL sensor                                GK, 20.02.2018  7-->8
% extended diary                                        GK, 15.01.2019  8-->9
% removed option shift_ctd_temp_by_seconds as this is now always optimized
% in the salinity optimization                          GK, 17.01.2019  9-->10
% handle main clock before clock correction as additional variable
%                                                       GK, 23.01.2019  10-->11


%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 09  at  ',datestr(now)])
disp('interpolate data to 1 sec timesteps')
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
% list of variables that will be stored
% variables have to be added here and all steps 3 and higher should be rerun
%
[slocum_variable_list,allglider_variable_list] = func_load_variable_list;
  

%
% interpolate data to 1 sec
%

% load result from previous step
load([op.deplname,'_yos'])


% loop over yos
for n=1:length(m_present_time)

  fprintf(1,'.',[]);

  % handle best pressure data
  best_pressure_datenum{n} = func_webbtime2mattime(best_pressure_timestamp{n});
  best_pressure{n} = best_pressure{n}*10; 

  % handle navigational data
  main_datenum{n} = func_webbtime2mattime(m_present_time{n});
  main_datenum_before_clock_correction{n} = func_webbtime2mattime(n_present_time_before_clock_correction{n});
  science_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
  latitude{n} = n_lat{n};
  longitude{n} = n_lon{n};
  pressure{n} = n_pressure{n}*10;
  good = find(~isnan(m_heading{n}));
  heading{n} = m_heading{n};
  heading{n}(good) = unwrap(heading{n}(good))/pi*180;
  fin{n} = m_fin{n}/pi*180;
  pitch{n} = m_pitch{n}/pi*180;
  roll{n} = m_roll{n}/pi*180;
  m_air_pump{n} = m_air_pump{n};
  internal_pressure{n} = (29.5-m_vacuum{n})*3386;
  thruster_power{n} = m_thruster_power{n};
  if sum(isnan(m_ballast_pumped{n}))<sum(isnan(m_de_oil_vol{n}))
    slocum_ballast_volume{n} = m_ballast_pumped{n};
  else
    slocum_ballast_volume{n} = m_de_oil_vol{n};
  end
  slocum_battery_position{n} = m_battpos{n};

  % handle CTD data
  ctd_conductivity{n} = sci_water_cond{n}*10;
  ctd_temperature{n} = sci_water_temp{n};
  ctd_pressure{n} = n_water_pressure{n}*10;
  if all(isnan(sci_ctd41cp_timestamp{n}))
    ctd_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
  else
    ctd_datenum{n} = func_webbtime2mattime(sci_ctd41cp_timestamp{n});
  end

  % handle Optode data
  if any(~isnan(sci_oxy3835_oxygen{n}))
    oxygen{n} = sci_oxy3835_oxygen{n};
  elseif any(~isnan(sci_oxy3835_wphase_oxygen{n}))
    oxygen{n} = sci_oxy3835_wphase_oxygen{n};
  elseif any(~isnan(sci_oxy4_oxygen{n}))
    oxygen{n} = sci_oxy4_oxygen{n};
  else
    oxygen{n} = sci_oxy3835_oxygen{n};  
  end
  bphase{n} = sci_oxy3835_wphase_bphase{n};
  dphase{n} = sci_oxy3835_wphase_dphase{n};
  if length(find(~isnan(sci_oxy4_c1rph{n}))) > length(find(~isnan(sci_oxy3835_wphase_bphase{n})))
    diffphase{n} = sci_oxy4_c1rph{n} - sci_oxy4_c2rph{n};
  else
    diffphase{n} = sci_oxy3835_wphase_bphase{n};
  end
  if any(~isnan(sci_oxy3835_oxygen{n}))
    if all(isnan(sci_oxy3835_timestamp{n}))
      optode_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
    else
      optode_datenum{n} = func_webbtime2mattime(sci_oxy3835_timestamp{n});
    end
    oxygen_temperature{n} = sci_oxy3835_temp{n};
  elseif any(~isnan(sci_oxy3835_wphase_oxygen{n}))
    if all(isnan(sci_oxy3835_wphase_timestamp{n}))
      optode_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
    else
      optode_datenum{n} = func_webbtime2mattime(sci_oxy3835_wphase_timestamp{n});
    end
    oxygen_temperature{n} = sci_oxy3835_wphase_temp{n};
  elseif any(~isnan(sci_oxy4_oxygen{n}))
    if all(isnan(sci_oxy4_timestamp{n}))
      optode_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
    else
      optode_datenum{n} = func_webbtime2mattime(sci_oxy4_timestamp{n});
    end
    oxygen_temperature{n} = sci_oxy3835_wphase_temp{n};
  else
    if all(isnan(sci_oxy3835_timestamp{n}))
      optode_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
    else
      optode_datenum{n} = func_webbtime2mattime(sci_oxy3835_timestamp{n}); 
    end
    oxygen_temperature{n} = sci_oxy3835_temp{n}; 
  end

  % handle FLNTU data  ISABELLE EDIT
  chlorophyll{n} = nan*sci_m_present_time{n};
  turbidity{n} = nan*sci_m_present_time{n};
  cdom{n} = nan*sci_m_present_time{n};
  rhodamin{n} = nan*sci_m_present_time{n};
  b700{n} = nan*sci_m_present_time{n};
  b470{n} = nan*sci_m_present_time{n};

  wetlabs_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
  wetlabs2_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
  wetlabs3_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});

  if any(~isnan(sci_flntu_chlor_units{n}))
    chlorophyll{n} = sci_flntu_chlor_units{n};
    turbidity{n} = sci_flntu_turb_units{n};
    cdom{n} = nan*sci_flntu_turb_units{n};
    b470{n} = sci_flntu_bb_units{n};
    b470{n} = sci_flntu_bb_units{n};

    if all(isnan(sci_flntu_timestamp{n}))
      wetlabs_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
    else
      wetlabs_datenum{n} = func_webbtime2mattime(sci_flntu_timestamp{n});
    end
  end
  if any(~isnan(sci_flbbcd_chlor_units{n}))
    chlorophyll{n} = sci_flbbcd_chlor_units{n};
    turbidity{n} = sci_flbbcd_bb_units{n};
    b700{n} = sci_flbbcd_bb_units{n};
    b470{n} = sci_flbbcd_bb_units{n};
    cdom{n} = sci_flbbcd_cdom_units{n};
    if all(isnan(sci_flbbcd_timestamp{n}))
      wetlabs_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
    else
      wetlabs_datenum{n} = func_webbtime2mattime(sci_flbbcd_timestamp{n});
    end
  end
  if any(~isnan(sci_bb2flsV2_chl_scaled{n}))
    chlorophyll{n} = sci_bb2flsV2_chl_scaled{n};
    if all(isnan(sci_bb2flsV2_timestamp{n}))
      wetlabs_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
    else
      wetlabs_datenum{n} = func_webbtime2mattime(sci_bb2flsV2_timestamp{n});
    end
  end
  if any(~isnan(sci_bb2flsV8_chl_scaled{n}))
    chlorophyll{n} = sci_bb2flsV8_chl_scaled{n};
    if all(isnan(sci_bb2flsV8_timestamp{n}))
      wetlabs_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
    else
      wetlabs_datenum{n} = func_webbtime2mattime(sci_bb2flsV8_timestamp{n});
    end
  end
  if any(~isnan(sci_bb2fls_cdom_scaled{n}))
    cdom{n} = sci_bb2fls_cdom_scaled{n};
    if all(isnan(sci_bb2fls_timestamp{n}))
      wetlabs2_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
    else
      wetlabs2_datenum{n} = func_webbtime2mattime(sci_bb2fls_timestamp{n});
    end
  end
  if any(~isnan(sci_fl2urrh_rhod_units{n}))
    rhodamin{n} = sci_fl2urrh_rhod_units{n};
    if all(isnan(sci_fl2urrh_timestamp{n}))
      wetlabs3_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
    else
      wetlabs3_datenum{n} = func_webbtime2mattime(sci_fl2urrh_timestamp{n});
    end
  end
  %   %isabelle add:
    if any(~isnan(sci_bb2flsv8_chl_scaled{n}))
    chlorophyll{n} = sci_bb2flsv8_chl_scaled{n};
    if all(isnan(sci_bb2flsv8_timestamp{n}))
      wetlabs2_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
    else
      wetlabs2_datenum{n} = func_webbtime2mattime(sci_bb2flsv8_timestamp{n});
    end
  end

  %   %isabelle add:
    if any(~isnan(sci_bb2flsv8_b470_scaled{n}))
    b470{n} = sci_bb2flsv8_b470_scaled{n};
    if all(isnan(sci_bb2flsv8_timestamp{n}))
      wetlabs2_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
    else
      wetlabs2_datenum{n} = func_webbtime2mattime(sci_bb2flsv8_timestamp{n});
    end
    end


  %   %isabelle add:
    if any(~isnan(sci_bb2flsv8_b700_scaled{n}))
    b700{n} = sci_bb2flsv8_b700_scaled{n};
    if all(isnan(sci_bb2flsv8_timestamp{n}))
      wetlabs2_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
    else
      wetlabs2_datenum{n} = func_webbtime2mattime(sci_bb2flsv8_timestamp{n});
    end
  end

  % handle SUNA data
  nitrate{n} = sci_suna_nitrate_um{n};
  if all(isnan(sci_suna_timestamp{n}))
    suna_datenum{n} = func_webbtime2mattime(sci_m_present_time{n});
  else
    suna_datenum{n} = func_webbtime2mattime(sci_suna_timestamp{n});
  end

  % Optode, CTD, and FLNTU datenum vectors can have same consecutive numbers
  % we replace the second of these with NaN
  good = find(~isnan(ctd_datenum{n}));
  bad = find(diff(ctd_datenum{n}(good))==0);
  if ~isempty(bad)
    ctd_datenum{n}(good(bad+1)) = nan;
  end
  good = find(~isnan(optode_datenum{n}));
  bad = find(diff(optode_datenum{n}(good))==0);
  if ~isempty(bad)
    optode_datenum{n}(good(bad+1)) = nan;
  end
  good = find(~isnan(wetlabs_datenum{n}));
  bad = find(diff(wetlabs_datenum{n}(good))==0);
  if ~isempty(bad)
    wetlabs_datenum{n}(good(bad+1)) = nan;
  end
  good = find(~isnan(wetlabs2_datenum{n}));
  bad = find(diff(wetlabs2_datenum{n}(good))==0);
  if ~isempty(bad)
    wetlabs2_datenum{n}(good(bad+1)) = nan;
  end
  good = find(~isnan(suna_datenum{n}));
  bad = find(diff(suna_datenum{n}(good))==0);
  if ~isempty(bad)
    suna_datenum{n}(good(bad+1)) = nan;
  end
  good = find(~isnan(wetlabs3_datenum{n}));
  bad = find(diff(wetlabs3_datenum{n}(good))==0);
  if ~isempty(bad)
    wetlabs3_datenum{n}(good(bad+1)) = nan;
  end

  % create the new time vector
  newtime = main_datenum{n}*86400;
  newtime = [ceil(newtime(1)):floor(newtime(end))]/86400;

  % set interpolation method
  metho = 'linear';

  % set variables to be interpolated from the science board timestamps to
  % the new time vector
  % these have to come first as they they store under different variable names
  % while the following ones overwrite
if 0
  inter_vars = {'ctd_pressure'};

  % interpolate these
  for m=1:length(inter_vars)
    eval(['dummy = ',inter_vars{m},'{n};'])
    good = find(~isnan(dummy.*science_datenum{n}));
    if length(good)>1
      dummy = interp1(science_datenum{n}(good),dummy(good),newtime,metho);
    else
      dummy = nan*newtime;
    end
    eval([inter_vars{m},'{n}=dummy;'])
  end
end

  % set variables to be interpolated from the navigational timestamps to
  % the new time vector
  inter_vars = {'slocum_battery_position',...
	'latitude',...
	'longitude',...
	'pressure',...
	'heading',...
	'fin',...
	'pitch',...
	'roll',...
        'm_air_pump',...
        'internal_pressure',...
        'thruster_power',...
	'slocum_ballast_volume',...
        'main_datenum_before_clock_correction'};

  % interpolate these
  for m=1:length(inter_vars)
    eval(['dummy = ',inter_vars{m},'{n};'])
    good = find(~isnan(dummy.*main_datenum{n}));
    if length(good)>1
      dummy = interp1(main_datenum{n}(good),dummy(good),newtime,metho);
    else
      dummy = nan*newtime;
    end
    eval([inter_vars{m},'{n}=dummy;'])
  end

  % set variables to be interpolated from the best pressure timestamps to
  % the new time vector
  inter_vars = {'best_pressure'};

  % interpolate these
  for m=1:length(inter_vars)
    eval(['dummy = ',inter_vars{m},'{n};'])
    good = find(~isnan(dummy.*best_pressure_datenum{n}));
    if length(good)>1
      dummy = interp1(best_pressure_datenum{n}(good),dummy(good),newtime,metho);
    else
      dummy = nan*newtime;
    end
    eval([inter_vars{m},'{n}=dummy;'])
  end

  % set variables to be interpolated from the CTD timestamps to
  % the new time vector
  inter_vars = {'ctd_pressure',...
	'ctd_temperature',...
	'ctd_conductivity'};

  % interpolate these
  for m=1:length(inter_vars)
    eval(['dummy = ',inter_vars{m},'{n};'])
    good = find(~isnan(dummy.*ctd_datenum{n}));
    if length(good)>1
      dummy = interp1(ctd_datenum{n}(good),dummy(good),newtime,metho);
    else
      dummy = nan*newtime;
    end
    eval([inter_vars{m},'{n}=dummy;'])
  end
    bad = find(isnan(ctd_pressure{n}+ctd_temperature{n}+ctd_conductivity{n}));
    if ~isempty(bad)
      ctd_pressure{n}(bad) = nan;
      ctd_temperature{n}(bad) = nan;
      ctd_conductivity{n}(bad) = nan;
    end

  % set variables to be interpolated from the Optode timestamps to
  % the new time vector
  inter_vars = {'oxygen',...
	'oxygen_temperature',...
	'dphase',...
	'diffphase',...
	'bphase'};

  % interpolate these
  for m=1:length(inter_vars)
    eval(['dummy = ',inter_vars{m},'{n};'])
    good = find(~isnan(dummy.*optode_datenum{n}));
    if length(good)>1
      dummy = interp1(optode_datenum{n}(good),dummy(good),newtime,metho);
    else
      dummy = nan*newtime;
    end
    eval([inter_vars{m},'{n}=dummy;'])
  end

  % set variables to be interpolated from the FLNTU/FLBBCD timestamps to
  % the new time vector  % ISABELLE EDIT
  inter_vars = {'chlorophyll',...
	'turbidity',...
        'cdom',...
        'b700',...
        'b470'};

  % interpolate these
  for m=1:length(inter_vars)
    eval(['dummy = ',inter_vars{m},'{n};'])
    good = find(~isnan(dummy.*wetlabs_datenum{n}));
    if length(good)>1
      dummy = interp1(wetlabs_datenum{n}(good),dummy(good),newtime,metho);
    else
      dummy = nan*newtime;
    end
    eval([inter_vars{m},'{n}=dummy;'])
  end

  % set variables to be interpolated from the FLNTU/FLBBCD timestamps to
  % the new time vector
  inter_vars = {'rhodamin'};

  % interpolate these
  for m=1:length(inter_vars)
    eval(['dummy = ',inter_vars{m},'{n};'])
    good = find(~isnan(dummy.*wetlabs3_datenum{n}));
    if length(good)>1
      dummy = interp1(wetlabs3_datenum{n}(good),dummy(good),newtime,metho);
    else
      dummy = nan*newtime;
    end
    eval([inter_vars{m},'{n}=dummy;'])
  end

  % set variables to be interpolated from the SUNA timestamps to
  % the new time vector
  inter_vars = {'nitrate'};

  % interpolate these
  for m=1:length(inter_vars)
    eval(['dummy = ',inter_vars{m},'{n};'])
    good = find(~isnan(dummy.*suna_datenum{n}));
    if length(good)>1
      dummy = interp1(suna_datenum{n}(good),dummy(good),newtime,metho);
    else
      dummy = nan*newtime;
    end
    eval([inter_vars{m},'{n}=dummy;'])
  end

  % last we replace the time vector
  new_datenum{n} = newtime;

end
fprintf(1,'\n',[]);
ctd_datenum = new_datenum;
optode_datenum = new_datenum;
science_datenum = new_datenum;
main_datenum = new_datenum;
wetlabs_datenum = new_datenum;
wetlabs2_datenum = new_datenum;
wetlabs3_datenum = new_datenum;
suna_datenum = new_datenum;
allglider_variable_list{end+1} = 'best_pressure';
allglider_variable_list{end+1} = 'main_datenum';
allglider_variable_list{end+1} = 'main_datenum_before_clock_correction';
allglider_variable_list{end+1} = 'm_air_pump';
allglider_variable_list{end+1} = 'diffphase';

%
% save the interpolated data
%
str = ['save(''',op.deplname,'_1sec'''];

for m=1:length(allglider_variable_list)
  if exist(allglider_variable_list{m},'var')
    str = [str,',''',allglider_variable_list{m},''''];
  end
end
str = [str,');'];
eval(str)


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 09  at  ',datestr(now)])
disp(' ')
diary off

