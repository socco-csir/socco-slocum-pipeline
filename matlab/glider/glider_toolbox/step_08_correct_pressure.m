function [] = step_08_ccorrect_pressure()
% function [] = step_08_ccorrect_pressure()
% 
% GEOMAR SVN $Id: step_08_correct_pressure.m 557 2019-01-25 10:53:26Z gkrahmann@geomar.de $
%
% try to correct the different shortcomings of the two pressure sensors
%
% version 5.1.0  last change 25.01.2023

% G.Krahmann, GEOMAR  Dec 2012

% serious rewrite                           GK, 26.03.2013  1-->2
% figure cleanup                            GK, 24.02.2017  2-->3
% add processing diary, remove v6 mat comp  GK, 30.01.2018  3-->4
% save plot and extend diary output         GK, 15.01.2019  4-->5
% added labels to plot                      GK, 25.01.2023  5-->5.1.0

% Unfortunately this is a bit involved and we have to do some operations twice.
% First we oversample the two pressure time series onto the same times with 1 second
% resolution. The interpolation for this is linear. 
% Then we create a filtered version of the nav pressure series.
% Then we extract surface offsets for both sensors. This is a bit tricky as we have to
% figure out, what part of the data is really at the surface.
% Then we look for jumps in the nav pressure data that occur at about 5m depth.
% Then we create a jump-corrected series of the nav pressure data.
% Then we reinterpolate and refilter this corrected series.
% Then we apply the offsets determined earlier.
% Then we determine a scaling factor between the two pressure series.

%
% add line to diary
%
diary diary_processing.txt
disp(['Starting step 08  at  ',datestr(now)])
disp('correct pressure records so that nav and science pressures match')
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
% load list of variables
%
slocum_variable_list = func_load_variable_list;


%
% load yos data and create single vectors for 3 variables
% we do this here from yos data since we might have applied corrections after the 1vector data
%
load([op.deplname,'_yos'],'dos_id','sci_water_pressure','m_argos_on');
dos_id = [dos_id{:}];
sci_water_pressure = [sci_water_pressure{:}];
m_argos_on = [m_argos_on{:}];
ids = unique(dos_id);


%
% We will first determine and correct the deck offsets of the science pressure sensor.
% For that we look for all the science pressure data when the Argos transmitter was on
% and determine the median value.
%
sci_pressure_offset = repmat(nan,1,length(ids));
for n=1:length(ids)
  good = find(dos_id==ids(n));
  argos_on = func_slocum_fill_vector(m_argos_on(good));
  ind = find(~isnan(argos_on) & sci_water_pressure(good)<1);
  if length(ind)>20
    sci_pressure_offset(n) = nmin(sci_water_pressure(good(ind)));
  end
end

figure(10)
plot(sci_pressure_offset,'x')
xlabel('Number of times when the glider was at the surface')
ylabel('Pressure of the science sensor at the surface [dbar]')


%
% fill gaps in the pressure offset series
%
for n=1:length(sci_pressure_offset)
  if isnan(sci_pressure_offset((n))) & n==1
    sci_pressure_offset(1) = nmedian(sci_pressure_offset);
  elseif isnan(sci_pressure_offset((n)))
    sci_pressure_offset(n) = sci_pressure_offset(n-1);
  end
end

if op.no_plot~=1
  figure(1)
  clf
  plot(sci_pressure_offset)
  title('Science pressure offset from surface data (to be subtracted from data)')
  ylabel('Surface offset [bar !]')
  xlabel('Data file index')
  print('-djpeg',['plots',filesep,'step_08_science_surface_pressure_correction.jpg'])
end

%
% load yo data
%
load([op.deplname,'_yos'])


%
% apply pressure offset
%
for n=1:length(sci_water_pressure)
  ind = find(dos_id{n}(1)==ids);
  sci_water_pressure2{n} = sci_water_pressure{n} - sci_pressure_offset(ind);
  sci_pressure_offset2(n) = sci_pressure_offset(ind);
end


% 
% Look for science pressures close to 0 and get the offset of the nav sensor
% for the same time steps.
% In theory the nav pressure should have been reset to 0 before the dive. 
% But sometimes there still appear to be offsets.
%
nav_pressure_offset = nan*sci_pressure_offset2;
for n=1:length(m_pressure)
  good_sci = find(~isnan(sci_water_pressure2{n} + sci_ctd41cp_timestamp{n}));
  if length(good_sci)>20
    good_nav = find(~isnan(m_pressure{n} + m_present_time{n}));
    if ~isempty(good_nav)
      n_pressure = interp1(m_present_time{n}(good_nav),m_pressure{n}(good_nav),sci_ctd41cp_timestamp{n}(good_sci));
    else
      n_pressure = m_pressure{n};
    end
    good = find(sci_water_pressure2{n}(good_sci)<0.2);
    if length(good)>20
      nav_pressure_offset(n) = nmedian(n_pressure(good));
    end
  end
end


%
% fill gaps in the nav pressure offset series
%
for n=1:length(nav_pressure_offset)
  if isnan(nav_pressure_offset((n))) & n==1
    nav_pressure_offset(1) = nmedian(nav_pressure_offset);
  elseif isnan(nav_pressure_offset((n)))
    nav_pressure_offset(n) = nav_pressure_offset(n-1);
  end
end


%
% apply pressure offset
%
clear n_pressure
for n=1:length(nav_pressure_offset)
  n_pressure{n} = m_pressure{n} - nav_pressure_offset(n);
end



%
% Look for spurious pressure jumps in nav pressure record.
% For that we go back to the original data.
%
jump_corrections = repmat(nan,[length(m_present_time),4]);
for n=1:length(m_present_time)
  if all(isnan(n_pressure{n}))
    ind = nan;
    corre = nan;
    down_or_up = nan;
  else
    [ind,corre,down_or_up] = func_find_pressure_jumps( n_pressure{n},...
      m_present_time{n},sci_water_pressure2{n},sci_ctd41cp_timestamp{n});
  end
  jump_corrections(n,:) = [ind,corre,down_or_up,dos_id{n}(1)];
end
% determine whether there is a significant amount of consistent pressure jumps
ids = unique(jump_corrections(:,4));
jump_corrections(:,2) = nans(jump_corrections(:,2),nan,0,'==');
for n=1:length(ids)
  ind = find(jump_corrections(:,4)==ids(n));
  if ~isempty(ind)
    ids_corr(n) = nmean(abs(jump_corrections(ind,2)));
    ids_corr_up(n) = nmean(abs(jump_corrections(ind(find(jump_corrections(ind,2)<0)),2)));
    ids_corr_down(n) = nmean(abs(jump_corrections(ind(find(jump_corrections(ind,2)>0)),2)));
  end
end
diary diary_processing.txt
if length(find(~isnan(ids_corr_down)))>0.5*length(ids_corr_down)
  disp('Found jumps in nav pressure data.')
  disp(['Typical jump down is : ',num2str(nmean(ids_corr_down)),' +/- ',num2str(nstd(ids_corr_down)),' bar'])
  disp(['Typical jump up is   : ',num2str(nmean(ids_corr_up)),' +/- ',num2str(nstd(ids_corr_up)),' bar'])
  disp(['Trying to correct them by subtracting ',num2str(nmean(ids_corr)),' bar from the inflicted pressure data.'])
  pressure_jump_correction = nmean(ids_corr);
else
  disp('Found no jumps in nav pressure data.')
  pressure_jump_correction = 0;
end
diary off

% Apply a correction to the pressure jumps.
% If we know when a jump occurs, then we apply it to the data before or after the jump.
% If we do not know when a jump occurs, then we apply it to the whole yo as we have to assume that
% it was a deep yo that never reached the pressure where the jumps occur.
for n=1:length(m_present_time)
  if jump_corrections(n,1)==1
    n_pressure{n}(:) = n_pressure{n}(:) + pressure_jump_correction;
  else
    ind = jump_corrections(n,1);
    down_or_up = jump_corrections(n,3);
    if down_or_up==1
      n_pressure{n}(1:ind) = n_pressure{n}(1:ind) - pressure_jump_correction;
    elseif down_or_up==-1
      n_pressure{n}(ind+1:end) = n_pressure{n}(ind+1:end) - pressure_jump_correction;
    end
  end
end


%
% interpolate nav pressure data onto the time steps of the science pressure
% and determine the offset and scaling factor of the nav pressure data
%
res = repmat(nan,[length(ids),1]);
sci_p = [sci_water_pressure{:}];
sci_t = [sci_ctd41cp_timestamp{:}];
nav_t = [m_present_time{:}];
nav_p = [n_pressure{:}];
dos_id_1 = [dos_id{:}];
for n=1:length(ids)
  good_sci = find(~isnan(sci_p + sci_t) & dos_id_1==ids(n));
  if length(good_sci)>20
    good_nav = find(~isnan(nav_p + nav_t) & dos_id_1==ids(n));
    if ~isempty(good_nav)
      n_pressure2 = interp1(nav_t(good_nav),nav_p(good_nav),sci_t(good_sci));
    else
      n_pressure2 = [];
    end
    good_sci2 = find(~isnan(n_pressure2));
    %dummy = ones(1,length(good_sci2));
    %res(n,:) = sci_p(good_sci(good_sci2))/[n_pressure2(good_sci2);dummy];
    if ~isempty(good_sci2)
      res(n) = sci_p(good_sci(good_sci2))/[n_pressure2(good_sci2)];
    else
      res(n) = nan;
    end
if 0 
figure(1)
clf
plot(sci_t(good_sci),sci_p(good_sci),'+',nav_t(good_nav),nav_p(good_nav),'x')
pause
end
  end
end
scaling_factor = nmedian(res(:,1));
diary diary_processing.txt
disp(['scaling nav pressure data by : ',num2str(scaling_factor)])
diary off

%
% apply the scaling factor to the jump-corrected nav pressure series
%
for n=1:length(m_pressure)
  n_pressure{n} = scaling_factor * n_pressure{n};
  %n_pressure{n} = scaling_factor * (n_pressure{n} - res(n,2));
  n_water_pressure{n} = sci_water_pressure2{n};
end


%
% compose 'best pressure vector' from CTD and  nav data
%
for n=1:length(n_pressure)
  if length(find(~isnan(n_water_pressure{n})))>20
    best_pressure{n} = n_water_pressure{n};
    best_pressure_timestamp{n} = sci_ctd41cp_timestamp{n};
  else
    best_pressure{n} = n_pressure{n};
    best_pressure_timestamp{n} = m_present_time{n};
  end
end


%
% load yo data again, add new data and resave
%
if 1
dummy = load([op.deplname,'_yos']);
dummy.n_pressure  = n_pressure;
dummy.n_water_pressure = n_water_pressure;
dummy.scaling_factor = scaling_factor;
dummy.nav_pressure_offset = nav_pressure_offset;
dummy.sci_pressure_offset = sci_pressure_offset;
dummy.best_pressure = best_pressure;
dummy.best_pressure_timestamp = best_pressure_timestamp;
dummy.pressure_jump_correction = pressure_jump_correction;
save([op.deplname,'_yos'],'-struct','dummy');
end


%
% plots to check the results
%
if 0
figure(4)
for n=1:length(m_present_time)
  disp(n)
  clf
  plot(m_present_time{n},m_pressure{n},'x',m_present_time{n},n_pressure{n},'+',...
    sci_ctd41cp_timestamp{n},sci_water_pressure{n},'o',m_present_time{n},m_argos_on{n},'*')
  pause
end
end


%
% add line to diary
%
diary diary_processing.txt
disp(['Ending step 08  at  ',datestr(now)])
disp(' ')
diary off

