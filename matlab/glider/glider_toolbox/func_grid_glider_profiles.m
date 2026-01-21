function [gridded] = func_grid_glider_profiles(data,mode)
% function [gridded] = func_grid_glider_profiles(data,mode)
% 
% GEOMAR SVN $Id: func_grid_glider_profiles.m 720 2020-09-10 09:32:43Z gkrahmann@geomar.de $
%
% grid 1 sec interpolated glider profile data to 1dbar steps
%
% input  :  data                - structure containing 1sec interpolated data
%           mode      [1]       -
%
% version 4  last change 06.03.2018

% G.Krahmann, GEOMAR,  Aug 2012

% bug fixes, new pressure basis                          GK, 03.04.2013  1-->2
% new mode 4 with fewer variables                        GK, 25.04.2014  2-->3
% new mode 5 to grid nitrate only                        GK, 06.03.2018  3-->4


%
% look for and load processing parameters
%
if exist('./processing_parameters.m') | exist('.\processing_parameters.m')
  processing_parameters
else
  error(['could not find file ''processing_parameters.m'' in current folder ',pwd])
end


%
% load the variable list
%
[slocum_variable_list,allglider_variable_list] = func_load_variable_list;


%
% which variables to grid
% mode 1 is for gridding all variables
% mode 2 is for gridding only a subset necessary for the salinity optimization evaluation
% mode 3 is for gridding only a subset necessary for the oxygen calibration evaluation
%
if mode==2
  allglider_variable_list = {'ctd_temperature','ctd_salinity','pressure_basis'};
elseif mode==3
  allglider_variable_list = {'oxygen','oxygen_calculated',...
    'oxygen_calculated_undelayed','oxygen_calculated_undelayed_filtered',...
    'oxygen_undelayed_filtered',...
    'oxygen_undelayed','pressure_basis'};
elseif mode==4
  allglider_variable_list = {'oxygen_calculated_undelayed','oxygen_undelayed',...
    'pressure_basis'};
elseif mode==5
  allglider_variable_list = {'nitrate','pressure_basis'};
else
  allglider_variable_list{end+1} = 'ctd_salinity_dpdt';
  allglider_variable_list{end+1} = 'ctd_salinity_model';
  allglider_variable_list{end+1} = 'ctd_salinity';
  allglider_variable_list{end+1} = 'pressure_basis';
end


%
% pick pressure variable
%
fn = fieldnames(data);
if any(strcmp('best_pressure',fn))
  disp('using best_pressure as pressure variable')
  data.pressure_basis = data.best_pressure;
elseif any(strcmp('pressure_filtered',fn))
  disp('using pressure_filtered as pressure variable')
  data.pressure_basis = data.pressure_filtered;
elseif any(strcmp('pressure',fn))
  disp('using pressure as pressure variable')
  data.pressure_basis = data.pressure;
else
  error('found neither pressure_filtered nor pressure as pressure variable')
end
fn{end+1} = 'pressure_basis';


%
% prepare regular pressure grid
%
pref = [0:ceil(nmax([data.pressure_basis{:}]))];


%
% loop over possible variables, and look for existing ones
%
count = 1;
for n=1:length(allglider_variable_list)
  if any(strcmp(allglider_variable_list{n},fn))
    variable_list{count} = allglider_variable_list{n};
    count = count+1;
  end
end
  

% 
% figure out the grid-index to which each of the glider's pressure levels belongs
% this will speed up the calculation later
%
for n=1:length(data.pressure_basis)
    pressure_index{n} = nans(round(data.pressure_basis{n}/diff(pref(1:2)))+1,nan,length(pref),'>');
end


%
% prepare output structures
%
nprofiles = length(pressure_index);
gridded.dummy = '1';
gridded = setfield(gridded,variable_list{1},repmat(nan,length(pref),nprofiles));
for n=2:length(variable_list)
  gridded = setfield(gridded,{1},variable_list{n},repmat(nan,length(pref),nprofiles));
end
gridded = rmfield(gridded,'dummy');


%
% loop over variables and create box averages
%
for n=1:length(variable_list)
  disp(variable_list{n})
  gridded_dummy = repmat(nan,length(pref),nprofiles);
  for m=1:nprofiles
    eval(['dummy = data.',variable_list{n},'{m};'])
    if any(~isnan(dummy))
      for k=1:length(pref)
        good = find(k==pressure_index{m});
        if ~isempty(good)
          gridded_dummy(k,m) = nmean(dummy(good));
        end
      end
    end
  end
  gridded = setfield(gridded,variable_list{n},gridded_dummy);
end


%
% loop over variables again and interpolate onto exact pressure grid
%
p = gridded.pressure_basis;
for n=1:length(variable_list)
  disp(variable_list{n})
  dummy = getfield(gridded,variable_list{n});
  dummy2 = nan*dummy;
  for m=1:size(dummy,2)
    good = find(~isnan(dummy(:,m).*p(:,m)));
    if length(good)>1
      dummy2(:,m) = interp1(p(good,m),dummy(good,m),pref,'linear',nan);
    end
  end
  gridded = setfield(gridded,variable_list{n},dummy2);
end

gridded.pressure = pref';

if mode==1
  m_water_vx = load([op.deplname,'_yos'],'m_water_vx');
  m_water_vy = load([op.deplname,'_yos'],'m_water_vy');
  for n=1:length(m_water_vx.m_water_vx)
    ind = find(~isnan(m_water_vx.m_water_vx{n}));
    if length(ind)>=1
      gridded.m_water_vx(n) = m_water_vx.m_water_vx{n}(ind(end));
      gridded.m_water_vy(n) = m_water_vy.m_water_vy{n}(ind(end));
    else
      gridded.m_water_vx(n) = nan;
      gridded.m_water_vy(n) = nan;
    end
  end

  save([op.deplname,'_gridded'],'-struct','gridded')
end
